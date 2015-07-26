#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cctype> // for tolower
#include "generator/Helpers.hpp"
#include "generator/WriterInfo.hpp"

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
          const std::string & cpuinfofile)
{
    options_ = options;
    finalam_ = finalam;
    ReadCPUFlags(cpuinfofile);
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
            nelements_ += NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]);
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
            os << indent1 << NewConstDoubleSet(it.first, it.second) << ";\n";  
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
            os << "    double contwork[SIMINT_NSHELL_SIMD * " << nelements_ << "] SIMINT_ALIGN_ARRAY;\n";
        os << "\n";

        os << "    // partition workspace into shells\n";

        size_t ptidx = 0;
        for(const auto & it : contq_)
        {
            if(it != finalam_)
            {
                os << "    double * const " << ArrVarName(it) << " = contwork + (SIMINT_NSHELL_SIMD * " << ptidx << ");\n";
                ptidx += NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]);
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

bool HasVRR(void) 
{
    return ((finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3]) > 0);
}

bool HasET(void) 
{
    return (finalam_[2]+finalam_[3] > 0);
}

bool HasHRR(void) 
{
    return (HasBraHRR() || HasKetHRR());
}

bool HasBraHRR(void) 
{
    return (finalam_[1] > 0);
}

bool HasKetHRR(void) 
{
    return (finalam_[3] > 0);
}


void ReadCPUFlags(const std::string & file)
{
    std::set<std::string> cpuflags_tmp;

    std::ifstream cpufile(file.c_str());
    if(!cpufile.is_open())
        throw std::runtime_error(std::string("Error - cannot open cpu flags file: ") + file);

    while(cpufile)
    {
        std::string s;
        cpufile >> s;

        // convert to lower case
        for(auto & it : s)
            it = tolower(it);

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
            intrinsicmap_["dbl_set"] = "MM512_SET1_PD";  // macro in intrinsics_kncni.h
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
            intrinsicmap_["dbl_set"] = "_mm256_set1_pd";
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
            intrinsicmap_["dbl_set"] = "_mm_set1_pd";
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

std::string DoubleSet(const std::string & dbl) 
{
    if(Intrinsics())
        return intrinsicmap_.at("dbl_set") + "(" + dbl + ")";
    else
        return dbl;
}

std::string NamedConstant(const std::string & cname)
{
    if(Intrinsics())
        return cname;
    else
        return DoubleSet(constants_.at(cname));
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

std::string NewDoubleSet(const std::string & var, const std::string & val) 
{
    std::stringstream ss;
    ss << DoubleType() << " " << var << " = " << DoubleSet(val);
    return ss.str();
}

std::string NewDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx) 
{
    std::stringstream ss;
    ss << DoubleType() << " " << var << " = " << DoubleLoad(ptr, idx);
    return ss.str();
}

std::string NewConstDoubleSet(const std::string & var, const std::string & val) 
{
    return std::string("const ") + NewDoubleSet(var, val);
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


void WriteAccumulation(std::ostream & os, QAM qam, int ncart)
{
    if(Intrinsics())
    {
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
                              << vecname.str() << ".d[0] + " << vecname.str() << ".d[1] +"
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
    }
    else
    {
        os << indent6 << "for(n = 0; n < " << ncart << "; n++)\n";
            os << indent7 << WriterInfo::PrimPtrName(qam) << "[n] += " << WriterInfo::PrimVarName(qam) << "[n];\n";
    }
}

} // close namespace WriterInfo
