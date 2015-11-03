#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cctype> // for tolower
#include "generator/Helpers.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/Ncart.hpp"

namespace {
        QAMSet contq_;  // set of contracted integral AM
        OptionsMap options_;
        size_t memory_;
        size_t nelements_;
        QAM finalam_;

        int simdlen_;
        std::map<std::string, std::string> intrinsicmap_;

        std::set<std::string> cpuflags_;
        std::vector<std::string> includes_;
        std::map<std::string, std::string> constants_;
}



namespace WriterInfo {

void Init(const OptionsMap & options, const QAM & finalam,
          const std::string & cpuflags)
{
    options_ = options;
    finalam_ = finalam;
    ReadCPUFlags(cpuflags);
}



void SetContQ(const QAMSet & topquartets)
{
    contq_.clear();

    // add hrr top level stuff
    for(const auto & it : topquartets)
        contq_.insert(it);

    // calculate the memory
    nelements_ = 0;
    for(const auto & it : contq_)
    {
        if(it != finalam_)
            nelements_ += NCART(it);
    }

    // add memory for AB_{xyz} and CD_{xyz}
    if(HasBraHRR())
        nelements_ += 3;
    if(HasKetHRR())
        nelements_ += 3;

    memory_ = nelements_ * sizeof(double); 
}



size_t MemoryReq(void) 
{
    return memory_;
}



int L(void) 
{
    return finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3];
}



int GetOption(int option) 
{
    if(options_.count(option) > 0)
        return options_.at(option);
    else
    {
        std::stringstream ss;
        ss << "No value for option: " << option;
        throw std::runtime_error(ss.str());
    }
}



bool IsContArray(const QAM & am) 
{
    return contq_.count(am);
}



bool IsFinalAM(const QAM & am) 
{
    return am == finalam_;
}


QAM FinalAM(void) 
{
    return finalam_;
}


void WriteIncludes(std::ostream & os) 
{
    for(const auto & it : includes_)
        os << "#include " << it << "\n";
}


void AddNamedConstant(const std::string & name, const std::string & val)
{
    constants_[name] = val;
}

void AddIntConstant(int i)
{
    std::stringstream ss, ssval;
    ss << "const_" << i;
    ssval << i;
    AddNamedConstant(ss.str(), ssval.str());
}


void WriteConstants(std::ostream & os)
{
    if(Intrinsics())
    {
        os << indent1 << "//Create constant simd vectors\n";
        for(const auto & it : constants_)
            os << indent1 << NewConstDoubleSet1(it.first, it.second) << ";\n";  
    }
}


std::string ArrVarName(const QAM & am, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]];
    return ss.str();
}

std::string ArrVarName(int am1, int am2, const std::string & ketstr, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << amchar[am1] << "_" << amchar[am2] << "_" << ketstr;
    return ss.str();
}

std::string ArrVarName(const std::string & brastr, int am3, int am4, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << brastr << "_" << amchar[am3] << "_" << amchar[am4];
    return ss.str();
}


std::string HRRVarName(const QAM & am)
{
    return ArrVarName(am, "HRR");
}


std::string HRRVarName(int am1, int am2, const std::string & ketstr)
{
    return ArrVarName(am1, am2, ketstr, "HRR");
}

std::string HRRVarName(const std::string & brastr, int am3, int am4)
{
    return ArrVarName(brastr, am3, am4, "HRR");
}


std::string PrimVarName(const QAM & am)
{
    return ArrVarName(am, "PRIM");
}


std::string PrimPtrName(const QAM & am)
{
    return ArrVarName(am, "PRIM_PTR");
}
        

void DeclareContwork(std::ostream & os) 
{

    if(memory_ > 0)
    {
        os << "    // Workspace for contracted integrals\n";

        if(memory_ > (size_t)GetOption(OPTION_STACKMEM))
            os << "    double * const contwork = ALLOC(SIMINT_NSHELL_SIMD * " << memory_ << ");\n";
        else
            os << "    double contwork[SIMINT_NSHELL_SIMD * " << nelements_ << "] SIMINT_ALIGN_ARRAY_DBL;\n";
        os << "\n";

        os << "    // partition workspace into shells\n";

        size_t ptidx = 0;
        for(const auto & it : contq_)
        {
            if(it != finalam_)
            {
                os << "    double * const " << ArrVarName(it) << " = contwork + (SIMINT_NSHELL_SIMD * " << ptidx << ");\n";
                ptidx += NCART(it);
            }
        
        }
        os << "\n";

        // AB_ and CD_ arrays
        if(HasBraHRR())
        {
            os << "    // Holds AB_{xyz} in a flattened fashion for later\n";
            os << "    double * const restrict AB_x = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict AB_y = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict AB_z = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "\n";
        }

        if(HasKetHRR())
        {
            os << "    // Holds CD_{xyz} in a flattened fashion for later\n";
            os << "    double * const restrict CD_x = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict CD_y = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict CD_z = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "\n";
        }
    }
}


void ZeroContWork(std::ostream & os) 
{
    if(memory_ > 0)
        os << "            memset(contwork, 0, SIMINT_NSHELL_SIMD * " << memory_ << ");\n";
}


void FreeContwork(std::ostream & os) 
{
    if((memory_ > 0) && (memory_ > (size_t)GetOption(OPTION_STACKMEM)))
    {
        os << "    // Free contracted work space\n";
        os << "    FREE(contwork);\n";
        os << "\n";
    }
}

bool HasBraVRR(void) 
{
    return ((finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3]) > 0);
}

bool HasKetVRR(void)
{
    return false; // for now
}

bool HasVRR(void) 
{
    return ((finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3]) > 0);
}

bool HasVRR_I(void)
{
    return HasVRR() && finalam_[0] >= finalam_[1];
}

bool HasVRR_J(void)
{
    return HasVRR() && finalam_[1] > finalam_[0];
}

bool HasVRR_K(void)
{
    return HasVRR() && finalam_[2] >= finalam_[3];
}

bool HasVRR_L(void)
{
    return HasVRR() && finalam_[3] > finalam_[2];
}

bool HasET(void) 
{
    return (!GetOption(OPTION_NOET)) && (finalam_[2]+finalam_[3] > 0);
}

bool HasHRR(void) 
{
    return (HasBraHRR() || HasKetHRR());
}

bool HasBraHRR(void) 
{
    return ( (finalam_[0] > 0) && (finalam_[1] > 0) );
}

bool HasBraHRR_I(void) // going from J -> I
{
    return HasBraHRR() && (finalam_[1] > finalam_[0]);
}

bool HasBraHRR_J(void) // going from I -> J
{
    return HasBraHRR() && (finalam_[0] >= finalam_[1]);
}

bool HasKetHRR(void) 
{
    return ( (finalam_[2] > 0) && (finalam_[3] > 0) );
}

bool HasKetHRR_K(void) // going from L -> K
{
    return HasKetHRR() && (finalam_[3] > finalam_[2]); 
}

bool HasKetHRR_L(void) // going from K -> L 
{
    return HasKetHRR() && (finalam_[2] >= finalam_[3]); 
}

void ReadCPUFlags(const std::string & flags)
{
    std::set<std::string> cpuflags_tmp;

    // separate by spaces rather than commas
    std::string spacesep = flags;
    for(auto & it : spacesep)
        if(it == ',')
            it = ' ';

    std::stringstream ss(spacesep);

    while(ss)
    {
        std::string s;
        ss >> s;

        // convert to lower case
        for(auto & it : s)
            it = (char)tolower(it);

        if(s.size())
            cpuflags_tmp.insert(s);
    }    

    /*
    std::cout << "Read " << cpuflags_tmp.size() << " flags\n";
    for(const auto & it : cpuflags_tmp)
        std::cout << "  \"" << it << "\"\n";
    std::cout << "\n";
    */

    // add equivalent ones where underscore is replaced by a period
    cpuflags_ = cpuflags_tmp;

    for(auto it : cpuflags_tmp)
    {
        for(auto & it2 : it)
        {
            if(it2 == '_')
                it2 = '.';
        }

        // its a set, so inserting duplicates is ok
        cpuflags_.insert(it);
    }

    
    std::cout << "Processed cpuflags into " << cpuflags_.size() << " flags\n";
    //for(const auto & it : cpuflags_)
    //    std::cout << "  \"" << it << "\"\n";
    std::cout << "\n";


    // defaults
    simdlen_ = 1;
    intrinsicmap_["dbl_type"] = "double";
    intrinsicmap_["cdbl_type"] = "const double";
    intrinsicmap_["dbl_set"] = "";
    intrinsicmap_["dbl_set1"] = "";
    intrinsicmap_["dbl_load"] = "";
    intrinsicmap_["dbl_store"] = "";
    intrinsicmap_["union_type"] = "";
    intrinsicmap_["sqrt"] = "sqrt";
    intrinsicmap_["pow"] = "pow";
    intrinsicmap_["exp"] = "exp";
    

    // determine simdlen, types, etc
    if(HasCPUFlag("kncni"))
    {
        simdlen_ = 8;  // 4 packed doubles
        includes_.push_back("\"vectorization/vectorization.h\"");

        if(Intrinsics())
        {
            includes_.push_back("\"vectorization/intrinsics_kncni.h\"");

            intrinsicmap_["dbl_type"] = "__m512d";
            intrinsicmap_["cdbl_type"] = "const __m512d";
            intrinsicmap_["dbl_set"] = "MM512_SET_PD";  // macro in intrinsics_kncni.h
            intrinsicmap_["dbl_set1"] = "MM512_SET1_PD";  // macro in intrinsics_kncni.h
            intrinsicmap_["dbl_load"] = "_mm512_load_pd";
            intrinsicmap_["dbl_store"] = "_mm512_store_pd";
            intrinsicmap_["union_type"] = "union double8";
            intrinsicmap_["sqrt"] = "MM512_SQRT_PD"; // macro in intrinsics_kncni.h
            intrinsicmap_["pow"] = "_mm512_pow_pd";
            intrinsicmap_["exp"] = "_mm512_exp_pd";

            intrinsicmap_["fmadd"] = "_mm512_fmadd_pd";
            intrinsicmap_["fmsub"] = "_mm512_fmsub_pd";
        }
    }
    else if(HasCPUFlag("avx"))
    {
        simdlen_ = 4;  // 4 packed doubles
        includes_.push_back("\"vectorization/vectorization.h\"");

        if(Intrinsics())
        {
            includes_.push_back("\"vectorization/intrinsics_avx.h\"");

            intrinsicmap_["dbl_type"] = "__m256d";
            intrinsicmap_["cdbl_type"] = "const __m256d";
            intrinsicmap_["dbl_set"] = "_mm256_set_pd";
            intrinsicmap_["dbl_set1"] = "_mm256_set1_pd";
            intrinsicmap_["dbl_load"] = "_mm256_load_pd";
            intrinsicmap_["dbl_store"] = "_mm256_store_pd";
            intrinsicmap_["union_type"] = "union double4";
            intrinsicmap_["sqrt"] = "_mm256_sqrt_pd";
            intrinsicmap_["pow"] = "_mm256_pow_pd";
            intrinsicmap_["exp"] = "_mm256_exp_pd";
        }

        if(HasCPUFlag("fma"))
        {
            intrinsicmap_["fmadd"] = "_mm256_fmadd_pd";
            intrinsicmap_["fmsub"] = "_mm256_fmsub_pd";
        }

    }
    else if(HasCPUFlag("sse2"))
    {
        simdlen_ = 2;
        includes_.push_back("\"vectorization/vectorization.h\"");

        if(Intrinsics())
        {
            includes_.push_back("\"vectorization/intrinsics_sse.h\"");

            intrinsicmap_["dbl_type"] = "__m128d";
            intrinsicmap_["cdbl_type"] = "const __m128d";
            intrinsicmap_["dbl_set"] = "_mm_set_pd";
            intrinsicmap_["dbl_set1"] = "_mm_set1_pd";
            intrinsicmap_["dbl_load"] = "_mm_load_pd";
            intrinsicmap_["dbl_store"] = "_mm_store_pd";
            intrinsicmap_["union_type"] = "union double2";
            intrinsicmap_["sqrt"] = "_mm_sqrt_pd";
            intrinsicmap_["pow"] = "_mm_pow_pd";
            intrinsicmap_["exp"] = "_mm_exp_pd";
        }
    }
    else
    {
        // keep the ones from Init()
        // but force intrinsics off
        options_[OPTION_INTRINSICS] = 0;
    }


}

bool HasCPUFlag(const std::string & flag) 
{
    return (cpuflags_.count(flag) > 0);
}
       

bool HasFMA(void)
{
    return HasCPUFlag("fma");
}
 
bool Intrinsics(void) 
{
    return options_[OPTION_INTRINSICS];
}


bool Scalar(void)
{
    return options_[OPTION_SCALAR];
}


int SimdLen(void) 
{
    return simdlen_;
}

int ByteAlign(void) 
{
    return 8*simdlen_;
}

std::string DoubleType(void) 
{
    return intrinsicmap_.at("dbl_type");
}

std::string ConstDoubleType(void) 
{
    return intrinsicmap_.at("cdbl_type");
}

std::string DoubleSet(const std::vector<std::string> & dbls)
{
    // todo - exception if dbls.size() != SimdLen?
    if(Intrinsics())
    {
        std::stringstream ss;
        ss << intrinsicmap_.at("dbl_set") <<  "(" << dbls[0];
        for(size_t i = 1; i < dbls.size(); ++i)
            ss << ", " << dbls[i];
        ss << ");\n";
        return ss.str();
    }
    else
        return dbls[0];
}

std::string DoubleSet1(const std::string & dbl) 
{
    if(Intrinsics())
        return intrinsicmap_.at("dbl_set1") + "(" + dbl + ")";
    else
        return dbl;
}

std::string NamedConstant(const std::string & cname)
{
    if(Intrinsics())
        return cname;
    else
        return DoubleSet1(constants_.at(cname));
}

std::string IntConstant(int i)
{
    std::stringstream ss;
    ss << "const_" << i;
    return NamedConstant(ss.str());
}

std::string DoubleLoad(const std::string & ptr, const std::string & idx) 
{
    if(Intrinsics())
    {
        std::string ptrstr = ptr;
        if(idx.size())
            ptrstr += " + " + idx;

        return intrinsicmap_.at("dbl_load") + "(" + ptrstr + ")";
    }
    else
        return ptr + "[" + idx + "]";
}

std::string DoubleStore(const std::string & var, const std::string & ptr, const std::string & idx) 
{
    if(Intrinsics())
    {
        std::string ptrstr = ptr;
        if(idx.size())
            ptrstr += " + " + idx;

        return intrinsicmap_.at("dbl_store") + "(" + ptrstr + ", " + var + ")";
    }
    else
        return ptr + "[" + idx + "] = " + var;
}

std::string NewDoubleSet1(const std::string & var, const std::string & val) 
{
    std::stringstream ss;
    ss << DoubleType() << " " << var << " = " << DoubleSet1(val);
    return ss.str();
}

std::string NewDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx) 
{
    std::stringstream ss;
    ss << DoubleType() << " " << var << " = " << DoubleLoad(ptr, idx);
    return ss.str();
}

std::string NewConstDoubleSet1(const std::string & var, const std::string & val) 
{
    return std::string("const ") + NewDoubleSet1(var, val);
}

std::string NewConstDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx) 
{
    return std::string("const ") + NewDoubleLoad(var, ptr, idx);
}

std::string UnionType(void)
{
    return intrinsicmap_.at("union_type");
}

std::string ConstUnionType(void)
{
    return std::string("const ") + UnionType();
}

std::string FMAdd(const std::string & a, const std::string & b, const std::string & c)
{
    return intrinsicmap_.at("fmadd") + "(" + a + ", " + b + ", " + c + ")";
}
std::string FMSub(const std::string & a, const std::string & b, const std::string & c)
{
    return intrinsicmap_.at("fmsub") + "(" + a + ", " + b + ", " + c + ")";
}
        
std::string Sqrt(const std::string & val) 
{
    return intrinsicmap_.at("sqrt") + "(" + val + ")";
}

std::string RSqrt(const std::string & val) 
{
    return std::string("1.0 / ") + Sqrt(val);
}
       
std::string Power(const std::string & base, const std::string & exp) 
{
    return intrinsicmap_.at("pow") + "(" + base + ", " + exp + ")";
}

std::string Exp(const std::string & exp) 
{
    return intrinsicmap_.at("exp") + "(" + exp + ")";
}


void WriteAccumulation(std::ostream & os)
{

    os << "\n\n";
    os << indent5 << "////////////////////////////////////\n";
    os << indent5 << "// Accumulate contracted integrals\n";
    os << indent5 << "////////////////////////////////////\n";
    if(Intrinsics())
    {
        os << indent5 << "if(hasoffset == 0)\n";
        os << indent5 << "{\n";

        for(const auto qam : contq_)
        {
            int ncart = NCART(qam);

            if(HasCPUFlag("avx"))
            {
                if(ncart > 3)
                {
                    os << indent6 << "for(n = 0; n < " << (ncart/4)*4 << "; n += 4)\n";
                    os << indent6 << "{\n";
                    os << indent7 << "__m256d t1 = _mm256_hadd_pd(" << WriterInfo::PrimVarName(qam) << "[n], " << WriterInfo::PrimVarName(qam) << "[n+1]);\n";
                    os << indent7 << "__m256d t2 = _mm256_hadd_pd(" << WriterInfo::PrimVarName(qam) << "[n+2], " << WriterInfo::PrimVarName(qam) << "[n+3]);\n";
                    os << indent7 << "__m256d t3 = _mm256_set_m128d(  _mm256_extractf128_pd(t2, 0) + _mm256_extractf128_pd(t2, 1), _mm256_extractf128_pd(t1, 0) + _mm256_extractf128_pd(t1, 1));\n";
                    os << indent7 << "_mm256_storeu_pd(" << WriterInfo::PrimPtrName(qam) << " + n, _mm256_loadu_pd(" << WriterInfo::PrimPtrName(qam) << " + n) + t3);\n";
                    os << indent6 << "}\n";
                }

                if((ncart%4) > 1)
                {
                    int n = (ncart/4)*4;

                    std::stringstream tmpname, ssename;
                    tmpname << "tmp_" << WriterInfo::PrimVarName(qam);
                    ssename << "sse_" << WriterInfo::PrimVarName(qam);

                    os << indent6 << "__m256d " << tmpname.str() << " = _mm256_hadd_pd(" << WriterInfo::PrimVarName(qam) << "[" << n << "], " 
                                  << WriterInfo::PrimVarName(qam) << "[" << n+1 << "]);\n";

                    os << indent6 << "__m128d " << ssename.str() << " = _mm256_extractf128_pd("
                                  << tmpname.str() << ", 0) + _mm256_extractf128_pd(" << tmpname.str() << ", 1);\n";
                    os << indent6 << "_mm_storeu_pd(" << WriterInfo::PrimPtrName(qam) << " + " << n << ", _mm_loadu_pd(" << WriterInfo::PrimPtrName(qam)
                                  << " + " << n << ") + " << ssename.str() << ");\n";
                }
                if((ncart%2) > 0)
                {
                    int n = ncart-1;
                    std::stringstream vecname;
                    vecname << "vec_" << WriterInfo::PrimVarName(qam);

                    os << indent6 << "union double4 " << vecname.str()
                                 << " = (union double4)" << WriterInfo::PrimVarName(qam) << "[" << n << "];\n";    
                    os << indent6 << WriterInfo::PrimPtrName(qam) << "[" << n << "] += "
                                  << vecname.str() << ".d[0] + " << vecname.str() << ".d[1] + "
                                  << vecname.str() << ".d[2] + " << vecname.str() << ".d[3];\n";
                }
            }
            else    // assume everyone has at least sse3
            {
                if(ncart > 1)
                {
                    os << indent6 << "for(n = 0; n < " << (ncart/2)*2 << "; n+=2)\n";  // integer division on purpose
                    os << indent6 << "{\n";
                    os << indent7 << "__m128d t1 = _mm_hadd_pd(" << WriterInfo::PrimVarName(qam) << "[n], " << WriterInfo::PrimVarName(qam) << "[n+1]);\n";
                    os << indent7 << "_mm_storeu_pd(" << WriterInfo::PrimPtrName(qam) << " + n, _mm_loadu_pd(" << WriterInfo::PrimPtrName(qam) << " + n) + t1);\n";
                    os << indent6 << "}\n";
                }
                if((ncart % 2) > 0)
                {
                    int n = (ncart/2)*2;
                    std::stringstream vecname;
                    vecname << "vec_" << WriterInfo::PrimVarName(qam);

                    os << indent6 << "union double2 " << vecname.str()
                                 << " = (" << WriterInfo::UnionType() << ")" << WriterInfo::PrimVarName(qam) << "[" << n << "];\n";    
                    os << indent6 << WriterInfo::PrimPtrName(qam) << "[" << n << "] += " << vecname.str() << ".d[0] + " << vecname.str() << ".d[1];\n";
                }
            }

            /*
            os << indent6 << "for(np = 0; np < " << ncart << "; ++np)\n";
            os << indent6 << "{\n";
            os << indent7 << ConstUnionType() << " tmp = (" << UnionType() << ")" << PrimVarName(qam) << "[np];\n";
            os << indent7 << "for(n = 0; n < SIMINT_SIMD_LEN; ++n)\n";
            os << indent8 << WriterInfo::PrimPtrName(qam) << "[np] += tmp.d[n];\n";
            os << indent6 << "}\n";
            */
        }
        os << indent5 << "}\n";
        os << indent5 << "else\n";
        os << indent5 << "{\n";

        for(const auto qam : contq_)
        {
            int ncart = NCART(qam);
            os << indent6 << "for(np = 0; np < " << ncart << "; ++np)\n";
            os << indent6 << "{\n";
            os << indent7 << WriterInfo::ConstUnionType() << " tmp = (" << WriterInfo::UnionType() << ")" << WriterInfo::PrimVarName(qam) << "[np];\n";
            os << indent7 << WriterInfo::PrimPtrName(qam) << "[np] += tmp.d[0];   // first offset is always zero\n";
            os << indent7 << "for(n = 1; n < SIMINT_SIMD_LEN; ++n)\n";
            os << indent8 << WriterInfo::PrimPtrName(qam) << "[shelloffsets[n]*" << ncart << "+np] += tmp.d[n];\n";
            os << indent6 << "}\n";
            os << indent6 << WriterInfo::PrimPtrName(qam) << " += shelloffsets[SIMINT_SIMD_LEN-1]*" << ncart << ";\n";
        }

        os << indent5 << "}\n";
    }
    else
    {
        for(const auto qam : contq_)
        {
            int ncart = NCART(qam);
            os << indent6 << "for(n = 0; n < " << ncart << "; n++)\n";
                os << indent7 << WriterInfo::PrimPtrName(qam) << "[n] += " << WriterInfo::PrimVarName(qam) << "[n];\n";
        }
    }
}




void WriteShellOffsets(std::ostream & os)
{
    os << indent5 << "// calculate the shell offsets\n";
    os << indent5 << "// these are the offset from the shell pointed to by cd\n";
    os << indent5 << "// for each element\n";
    os << indent5 << "int shelloffsets[SIMINT_SIMD_LEN] = {0};\n";
    os << indent5 << "int hasoffset = 0;\n";
    os << indent5 << "if((iprimcd + SIMINT_SIMD_LEN) >= nprim_icd)\n";
    os << indent5 << "{\n";

    os << indent6 << "// Handle if the first element of the vector is a new shell\n";
    os << indent6 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os << indent6 << "{\n";
    os << indent7 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    
    for(const auto qam : contq_)
        os << indent7 << WriterInfo::PrimPtrName(qam) << " += " << NCART(qam) << ";\n";

    os << indent6 << "}\n";
    os << indent6 << "iprimcd++;\n";

    os << indent6 << "for(n = 1; n < SIMINT_SIMD_LEN; ++n)\n";
    os << indent6 << "{\n";
    os << indent7 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os << indent7 << "{\n";
    os << indent8 << "hasoffset = 1;\n";
    os << indent8 << "shelloffsets[n] = shelloffsets[n-1] + 1;\n";
    os << indent8 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    os << indent7 << "}\n";
    os << indent7 << "else\n";
    os << indent8 << "shelloffsets[n] = shelloffsets[n-1];\n";
    os << indent7 << "iprimcd++;\n";
    os << indent6 << "}\n";
    os << indent5 << "}\n";
    os << indent5 << "else\n";
    os << indent6 << "iprimcd += SIMINT_SIMD_LEN;\n\n";
}




} // close namespace WriterInfo
