#include <stdexcept>
#include <fstream>
#include <iostream>
#include <cctype> // for tolower
#include "generator/Helpers.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/Ncart.hpp"
#include "generator/Naming.hpp"

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


void WriteIncludes(std::ostream & os) 
{
    for(const auto & it : includes_)
        os << "#include " << it << "\n";
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


void WriteAccumulation(std::ostream & os)
{

    os << "\n\n";
    os << indent5 << "////////////////////////////////////\n";
    os << indent5 << "// Accumulate contracted integrals\n";
    os << indent5 << "////////////////////////////////////\n";
    if(Intrinsics() && !Scalar())
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
                    os << indent7 << "__m256d t1 = _mm256_hadd_pd(" << PrimVarName(qam) << "[n], " << PrimVarName(qam) << "[n+1]);\n";
                    os << indent7 << "__m256d t2 = _mm256_hadd_pd(" << PrimVarName(qam) << "[n+2], " << PrimVarName(qam) << "[n+3]);\n";
                    os << indent7 << "__m256d t3 = _mm256_set_m128d(  _mm256_extractf128_pd(t2, 0) + _mm256_extractf128_pd(t2, 1), _mm256_extractf128_pd(t1, 0) + _mm256_extractf128_pd(t1, 1));\n";
                    os << indent7 << "_mm256_storeu_pd(" << PrimPtrName(qam) << " + n, _mm256_loadu_pd(" << PrimPtrName(qam) << " + n) + t3);\n";
                    os << indent6 << "}\n";
                }

                if((ncart%4) > 1)
                {
                    int n = (ncart/4)*4;

                    std::stringstream tmpname, ssename;
                    tmpname << "tmp_" << PrimVarName(qam);
                    ssename << "sse_" << PrimVarName(qam);

                    os << indent6 << "__m256d " << tmpname.str() << " = _mm256_hadd_pd(" << PrimVarName(qam) << "[" << n << "], " 
                                  << PrimVarName(qam) << "[" << n+1 << "]);\n";

                    os << indent6 << "__m128d " << ssename.str() << " = _mm256_extractf128_pd("
                                  << tmpname.str() << ", 0) + _mm256_extractf128_pd(" << tmpname.str() << ", 1);\n";
                    os << indent6 << "_mm_storeu_pd(" << PrimPtrName(qam) << " + " << n << ", _mm_loadu_pd(" << PrimPtrName(qam)
                                  << " + " << n << ") + " << ssename.str() << ");\n";
                }
                if((ncart%2) > 0)
                {
                    int n = ncart-1;
                    std::stringstream vecname;
                    vecname << "vec_" << PrimVarName(qam);

                    os << indent6 << "union double4 " << vecname.str()
                                 << " = (union double4)" << PrimVarName(qam) << "[" << n << "];\n";    
                    os << indent6 << PrimPtrName(qam) << "[" << n << "] += "
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
                    os << indent7 << "__m128d t1 = _mm_hadd_pd(" << PrimVarName(qam) << "[n], " << PrimVarName(qam) << "[n+1]);\n";
                    os << indent7 << "_mm_storeu_pd(" << PrimPtrName(qam) << " + n, _mm_loadu_pd(" << PrimPtrName(qam) << " + n) + t1);\n";
                    os << indent6 << "}\n";
                }
                if((ncart % 2) > 0)
                {
                    int n = (ncart/2)*2;
                    std::stringstream vecname;
                    vecname << "vec_" << PrimVarName(qam);

                    os << indent6 << "union double2 " << vecname.str()
                                 << " = (" << WriterInfo::UnionType() << ")" << PrimVarName(qam) << "[" << n << "];\n";    
                    os << indent6 << PrimPtrName(qam) << "[" << n << "] += " << vecname.str() << ".d[0] + " << vecname.str() << ".d[1];\n";
                }
            }

            /*
            os << indent6 << "for(np = 0; np < " << ncart << "; ++np)\n";
            os << indent6 << "{\n";
            os << indent7 << ConstUnionType() << " tmp = (" << UnionType() << ")" << PrimVarName(qam) << "[np];\n";
            os << indent7 << "for(n = 0; n < SIMINT_SIMD_LEN; ++n)\n";
            os << indent8 << PrimPtrName(qam) << "[np] += tmp.d[n];\n";
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
            os << indent7 << WriterInfo::ConstUnionType() << " tmp = (" << WriterInfo::UnionType() << ")" << PrimVarName(qam) << "[np];\n";
            os << indent7 << PrimPtrName(qam) << "[np] += tmp.d[0];   // first offset is always zero\n";
            os << indent7 << "for(n = 1; n < SIMINT_SIMD_LEN; ++n)\n";
            os << indent8 << PrimPtrName(qam) << "[shelloffsets[n]*" << ncart << "+np] += tmp.d[n];\n";
            os << indent6 << "}\n";
            os << indent6 << PrimPtrName(qam) << " += shelloffsets[SIMINT_SIMD_LEN-1]*" << ncart << ";\n";
        }

        os << indent5 << "}\n";
    }
    else
    {
        for(const auto qam : contq_)
        {
            int ncart = NCART(qam);
            os << indent6 << "for(n = 0; n < " << ncart << "; n++)\n";
                os << indent7 << PrimPtrName(qam) << "[n] += " << PrimVarName(qam) << "[n];\n";
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
        os << indent7 << PrimPtrName(qam) << " += " << NCART(qam) << ";\n";

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


void WriteShellOffsets_Scalar(std::ostream & os)
{
    os << indent5 << "// Move pointers if this is the end of a shell\n";
    os << indent5 << "// Handle if the first element of the vector is a new shell\n";
    os << indent5 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os << indent5 << "{\n";
    os << indent6 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    
    for(const auto qam : contq_)
        os << indent7 << PrimPtrName(qam) << " += " << NCART(qam) << ";\n";

    os << indent5 << "}\n";
    os << indent5 << "iprimcd++;\n";
}



} // close namespace WriterInfo
