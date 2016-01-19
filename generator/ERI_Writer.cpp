#include "generator/ERI_Writer.hpp"
#include "generator/Boys.hpp"
#include "generator/VRR_Writer.hpp"
#include "generator/ET_Writer.hpp"
#include "generator/HRR_Writer.hpp"
#include "generator/ERIGeneratorInfo.hpp"
#include "generator/Classes.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"


void DeclareContwork(std::ostream & os, const ERIGeneratorInfo & info)
{
    if(info.ContMemoryReq() == 0)
        return;

    size_t contmem = info.ContMemoryReq();

    os << indent1 << "//Workspace for contracted integrals\n";
    if(info.UseHeap())
        os << indent1 << "double * const constwork = ALLOC(SIMINT_NSHELL_SIMD * " << contmem << ");\n\n";
    else
        os << indent1 << "double contwork[SIMINT_NSHELL_SIMD * " << contmem << "] SIMINT_ALIGN_ARRAY_DBL;\n\n";

    os << indent1 << "// partition workspace into shells\n";
    size_t ptidx = 0;

    for(const auto & it : info.GetContQ())
    {
        if(!info.IsFinalAM(it))
        {
            os << indent1 << "double * const " << ArrVarName(it) << " = contwork + (SIMINT_NSHELL_SIMD * " << ptidx << ");\n";
            ptidx += NCART(it);
        }
    }

    os << "\n";

}


void ZeroContWork(std::ostream & os, const ERIGeneratorInfo & info)
{
    size_t contmem = info.ContMemoryReq();
    if(contmem > 0)
        os << indent3 << "memset(contwork, 0, SIMINT_NSHELL_SIMD * " << contmem << ");\n";
    
}


void FreeContwork(std::ostream & os, const ERIGeneratorInfo & info)
{
    size_t contmem = info.ContMemoryReq();

    if(contmem > 0 && info.UseHeap())
    {
        os << indent1 << "// Free contracted workspace\n";
        os << indent1 << "FREE(contwork);\n\n";
    }       
}

void WriteShellOffsets(std::ostream & os, const ERIGeneratorInfo & info)
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
    
    for(const auto qam : info.GetContQ())
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


void WriteShellOffsets_Scalar(std::ostream & os, const ERIGeneratorInfo & info)
{
    os << indent5 << "// Move pointers if this is the end of a shell\n";
    os << indent5 << "// Handle if the first element of the vector is a new shell\n";
    os << indent5 << "if(iprimcd >= nprim_icd && ((icd+1) < nshellbatch))\n";
    os << indent5 << "{\n";
    os << indent6 << "nprim_icd += Q.nprim12[cd + (++icd)];\n";
    
    for(const auto qam : info.GetContQ())
        os << indent7 << PrimPtrName(qam) << " += " << NCART(qam) << ";\n";

    os << indent5 << "}\n";
    os << indent5 << "iprimcd++;\n";
}


void WriteAccumulation(std::ostream & os, const ERIGeneratorInfo & info)
{
    const VectorInfo & vinfo = info.GetVectorInfo();

    os << "\n\n";
    os << indent5 << "////////////////////////////////////\n";
    os << indent5 << "// Accumulate contracted integrals\n";
    os << indent5 << "////////////////////////////////////\n";
    if(info.Vectorized())
    {
        os << indent5 << "if(hasoffset == 0)\n";
        os << indent5 << "{\n";

        for(const auto qam : info.GetContQ())
        {
            int ncart = NCART(qam);

            if(info.HasCPUFlag("avx"))
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

                    std::string tmpname = StringBuilder("tmp_", PrimVarName(qam));
                    std::string ssename = StringBuilder("sse_", PrimVarName(qam));

                    os << indent6 << "__m256d " << tmpname << " = _mm256_hadd_pd(" << PrimVarName(qam) << "[" << n << "], " 
                                  << PrimVarName(qam) << "[" << n+1 << "]);\n";

                    os << indent6 << "__m128d " << ssename << " = _mm256_extractf128_pd("
                                  << tmpname << ", 0) + _mm256_extractf128_pd(" << tmpname << ", 1);\n";
                    os << indent6 << "_mm_storeu_pd(" << PrimPtrName(qam) << " + " << n << ", _mm_loadu_pd(" << PrimPtrName(qam)
                                  << " + " << n << ") + " << ssename << ");\n";
                }
                if((ncart%2) > 0)
                {
                    int n = ncart-1;
                    std::string vecname = StringBuilder("vec_", PrimVarName(qam));

                    os << indent6 << "union double4 " << vecname
                                 << " = (union double4)" << PrimVarName(qam) << "[" << n << "];\n";    
                    os << indent6 << PrimPtrName(qam) << "[" << n << "] += "
                                  << vecname << ".d[0] + " << vecname << ".d[1] + "
                                  << vecname << ".d[2] + " << vecname << ".d[3];\n";
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
                    std::string vecname = StringBuilder("vec_", PrimVarName(qam));

                    os << indent6 << "union double2 " << vecname
                                 << " = (" << vinfo.UnionType() << ")" << PrimVarName(qam) << "[" << n << "];\n";    
                    os << indent6 << PrimPtrName(qam) << "[" << n << "] += " << vecname << ".d[0] + " << vecname << ".d[1];\n";
                }
            }
        }
        os << indent5 << "}\n";
        os << indent5 << "else\n";
        os << indent5 << "{\n";

        for(const auto qam : info.GetContQ())
        {
            int ncart = NCART(qam);
            os << indent6 << "for(np = 0; np < " << ncart << "; ++np)\n";
            os << indent6 << "{\n";
            os << indent7 << vinfo.ConstUnionType() << " tmp = (" << vinfo.UnionType() << ")" << PrimVarName(qam) << "[np];\n";
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
        for(const auto qam : info.GetContQ())
        {
            int ncart = NCART(qam);
            os << indent6 << "for(n = 0; n < " << ncart << "; n++)\n";
                os << indent7 << PrimPtrName(qam) << "[n] += " << PrimVarName(qam) << "[n];\n";
        }
    }
}


void WriteFile_Permute(std::ostream & os, ERIGeneratorInfo & info)
{
    QAM am = info.FinalAM();
    QAM tocall = am;

    bool swap_ab = false;
    bool swap_cd = false;

    // TODO - more thoroughly check?
    if(am[0] == 0)
    {
        swap_ab = true;
        std::swap(tocall[0], tocall[1]);
    }
    if(am[2] == 0)
    {
        swap_cd = true;
        std::swap(tocall[2], tocall[3]);
    }

    std::string funcline = StringBuilder("int eri_", "int eri_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string indent(funcline.length(), ' ');

    // start output to the file
    // we only need this one include
    os << "#include \"eri/eri.h\"\n";
    os << "\n";

    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair P,\n";
    os << indent << "struct multishell_pair Q,\n";
    os << indent << "double * const restrict " << ArrVarName(am) << ")\n";
    os << "{\n";
    os << indent1 << "// Can be accomplished by swapping some variables\n";
    os << indent1 << "// and calling another function\n";
    os << indent1 << "// Note that the struct was passed by copy\n";
    os << "\n";
    os << indent1 << "double * tmp;\n";
    if(swap_ab)
    {
        os << indent1 << "tmp = P.PA_x;   P.PA_x = P.PB_x;   P.PB_x = tmp;\n";
        os << indent1 << "tmp = P.PA_y;   P.PA_y = P.PB_y;   P.PB_y = tmp;\n";
        os << indent1 << "tmp = P.PA_z;   P.PA_z = P.PB_z;   P.PB_z = tmp;\n";
    }
    if(swap_cd)
    {
        os << indent1 << "tmp = Q.PA_x;   Q.PA_x = Q.PB_x;   Q.PB_x = tmp;\n";
        os << indent1 << "tmp = Q.PA_y;   Q.PA_y = Q.PB_y;   Q.PB_y = tmp;\n";
        os << indent1 << "tmp = Q.PA_z;   Q.PA_z = Q.PB_z;   Q.PB_z = tmp;\n";
    }

    os << indent1 << "return eri_" << amchar[tocall[0]] << "_" 
                                   << amchar[tocall[1]] << "_" 
                                   << amchar[tocall[2]] << "_" 
                                   << amchar[tocall[3]] << "(P, Q, " << ArrVarName(am) << ");\n"; 


    os << "}\n";
    os << "\n";
}


void WriteFile(std::ostream & os,
               ERIGeneratorInfo & info,
               const BoysGen & bg,
               const VRR_Writer & vrr_writer,
               const ET_Writer & et_writer,
               const HRR_Writer & hrr_writer)
{
    const VectorInfo & vinfo = info.GetVectorInfo();
    const QAM am = info.FinalAM();
    const int ncart = NCART(am);

    


    // some helper bools
    bool hashrr = hrr_writer.HasHRR();
    bool hasbrahrr = hrr_writer.HasBraHRR();
    bool haskethrr = hrr_writer.HasKetHRR();
    bool inline_hrr = (hashrr && hrr_writer.IsInline());

    bool hasbravrr = vrr_writer.HasBraVRR();
    bool hasketvrr = vrr_writer.HasKetVRR();

    bool haset = et_writer.HasET();
    bool hasbraet = et_writer.HasBraET(); 
    bool hasketet = et_writer.HasKetET(); 

    bool hasoneoverp = (hasbravrr || hasbraet);
    bool hasoneoverq = (hasketvrr || hasketet);
    bool hasoneover2p = (hasbraet || (hasbravrr && (am[0]+am[1]) > 1)); 
    bool hasoneover2q = (hasketet || (hasketvrr && (am[2]+am[3]) > 1)); 
    bool hasoneover2pq = haset && (am[0] + am[1] > 0) && (am[2] + am[3] > 0);


    std::string dbltype = vinfo.DoubleType();
    std::string cdbltype = vinfo.ConstDoubleType();

    // we need a constant one for 1/x
    info.AddIntegerConstant(1);

    // add includes
    info.AddInclude("<string.h>");
    info.AddInclude("<math.h>");
    info.AddInclude("\"eri/eri.h\"");
    bg.AddIncludes(info);

    // add constants from the various steps
    bg.AddConstants(info);
    et_writer.AddConstants(info);
    vrr_writer.AddConstants(info);
    hrr_writer.AddConstants(info);

    // need these factors sometimes
    if(hasoneover2p || hasoneover2q || hasoneover2pq)
        info.AddNamedConstant("one_half", "0.5");



    ///////////////////////////////////////
    // Beginning of file writing
    ///////////////////////////////////////

    // Write out all the includes
    for(const auto & it : info.GetIncludes())
        os << "#include " << it << "\n";


    //////////////////////////////
    // Function name & signature
    //////////////////////////////
    std::string funcline = StringBuilder("int eri_", "int eri_", amchar[am[0]], "_", amchar[am[1]], "_" , amchar[am[2]], "_", amchar[am[3]], "(");
    std::string indent(funcline.length(), ' ');


    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair const P,\n";
    os << indent << "struct multishell_pair const Q,\n";
    os << indent << "double * const restrict " << ArrVarName(am) << ")\n";
    os << "{\n";
    os << "\n";


    ///////////////////////////////////
    // NOW IN THE ACTUAL ERI FUNCTION
    ///////////////////////////////////

    // If there is no HRR, integrals are accumulated from inside the primitive loop
    // directly into the final integral array that was passed into this function, so it must be zeroed first
    if(!hashrr)
        os << indent1 << "memset(" << ArrVarName(am) << ", 0, P.nshell12 * Q.nshell12 * " << ncart << " * sizeof(double));\n";
    os << "\n";


    // abcd = index within simd loop, 
    os << indent1 << "int ab, cd, cdbatch, abcd;\n";
    os << indent1 << "int istart, jstart;\n";
    os << indent1 << "int iprimcd, nprim_icd, icd;\n";
    os << indent1 << "int i, j, n;\n";

    // real_abcd is the absolute actual abcd in terms of all the shells that we are doing
    // (only needed if we do HRR)
    if(hashrr)
        os << indent1 << "int real_abcd;\n";


    // Needed for determining offsets
    // But that's only if we are vectorizing
    if(info.Vectorized())
        os << indent1 << "int np;\n";
    

    // Needed only if we are doing inline HRR
    if(inline_hrr)
    {
        if(hasbrahrr)
            os << indent1 << "int iket;\n";
        if(haskethrr)
            os << indent1 << "int ibra;\n";
    }

    os << "\n";

    // Declare the temporary space 
    // Only needed if we are doing HRR
    if(hashrr)
        DeclareContwork(os, info);



    // Write out all the constants 
    os << indent1 << "//Create constants\n";
    for(const auto & it : info.GetConstants())
        os << indent1 << vinfo.NewConstDoubleSet1(it.first, it.second) << ";\n";

    


    os << "\n\n";
    os << indent1 << "////////////////////////////////////////\n";
    os << indent1 << "// Loop over shells and primitives\n";
    os << indent1 << "////////////////////////////////////////\n";
    os << "\n";
    if(hashrr)
        os << indent1 << "real_abcd = 0;\n";
    else
        os << indent1 << "abcd = 0;\n";

    os << indent1 << "istart = 0;\n";
    os << indent1 << "for(ab = 0; ab < P.nshell12; ++ab)\n";
    os << indent1 << "{\n";

    os << indent2 << "const int iend = istart + P.nprim12[ab];\n";
    os << "\n";

    os << indent2 << "cd = 0;\n";
    os << indent2 << "jstart = 0;\n";
    os << "\n";

    os << indent2 << "for(cdbatch = 0; cdbatch < Q.nbatch; ++cdbatch)\n";
    os << indent2 << "{\n";
    os << indent3 << "const int nshellbatch = ((cd + SIMINT_NSHELL_SIMD) > Q.nshell12) ? Q.nshell12 - cd : SIMINT_NSHELL_SIMD;\n";

    os << indent3 << "const int jend = jstart + Q.nbatchprim[cdbatch];\n";


    if(hashrr)
    {
        ZeroContWork(os, info);
        os << indent3 << "abcd = 0;\n";
        os << "\n";
    }

    os << indent3 << "for(i = istart; i < iend; ++i)\n";
    os << indent3 << "{\n";
    os << "\n";

    os << indent4 << "icd = 0;\n";
    os << indent4 << "iprimcd = 0;\n";
    os << indent4 << "nprim_icd = Q.nprim12[cd];\n";

    if(vrr_writer.HasVRR())
        vrr_writer.DeclarePrimPointers(os, info);
    if(et_writer.HasET())
        et_writer.DeclarePrimPointers(os, info);
    os << "\n";

    os << indent4 << "// Load these one per loop over i\n";
    os << indent4 << vinfo.NewConstDoubleSet1("P_alpha", "P.alpha[i]") << ";\n";
    os << indent4 << vinfo.NewConstDoubleSet1("P_prefac", "P.prefac[i]") << ";\n";
    os << indent4 << vinfo.NewConstDoubleSet1("P_x", "P.x[i]") << ";\n";
    os << indent4 << vinfo.NewConstDoubleSet1("P_y", "P.y[i]") << ";\n";
    os << indent4 << vinfo.NewConstDoubleSet1("P_z", "P.z[i]") << ";\n";

    if(hasbravrr)
    {
        if(vrr_writer.HasVRR_I())
        {
            os << indent4 << vinfo.NewConstDoubleSet1("P_PA_x", "P.PA_x[i]") << ";\n";
            os << indent4 << vinfo.NewConstDoubleSet1("P_PA_y", "P.PA_y[i]") << ";\n";
            os << indent4 << vinfo.NewConstDoubleSet1("P_PA_z", "P.PA_z[i]") << ";\n";
        }
        else
        {
            os << indent4 << vinfo.NewConstDoubleSet1("P_PB_x", "P.PB_x[i]") << ";\n";
            os << indent4 << vinfo.NewConstDoubleSet1("P_PB_y", "P.PB_y[i]") << ";\n";
            os << indent4 << vinfo.NewConstDoubleSet1("P_PB_z", "P.PB_z[i]") << ";\n";
        }
    }

    if(hasketet)
    {
        os << indent4 << vinfo.NewConstDoubleSet1("P_bAB_x", "P.bAB_x[i]") << ";\n";
        os << indent4 << vinfo.NewConstDoubleSet1("P_bAB_y", "P.bAB_y[i]") << ";\n";
        os << indent4 << vinfo.NewConstDoubleSet1("P_bAB_z", "P.bAB_z[i]") << ";\n";
    }

    os << "\n";


    os << indent4 << "for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)\n";
    os << indent4 << "{\n";

    if(info.Vectorized())
        WriteShellOffsets(os, info);
    else
        WriteShellOffsets_Scalar(os, info);

    os << "\n";

    if(vrr_writer.HasVRR())
        vrr_writer.DeclarePrimArrays(os, info);
    if(et_writer.HasET())
        et_writer.DeclarePrimArrays(os, info);

    os << indent5 << vinfo.NewConstDoubleLoad("Q_alpha", "Q.alpha", "j") << ";\n";
    os << indent5 << cdbltype << " PQalpha_mul = P_alpha * Q_alpha;\n";
    os << indent5 << cdbltype << " PQalpha_sum = P_alpha + Q_alpha;\n";
    os << indent5 << cdbltype << " one_over_PQalpha_sum = " << "const_1 / PQalpha_sum;\n";
    os << "\n";
    os << "\n";
    os << indent5 << "/* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os << indent5 << cdbltype << " PQ_x = P_x - " << vinfo.DoubleLoad("Q.x", "j") << ";\n";
    os << indent5 << cdbltype << " PQ_y = P_y - " << vinfo.DoubleLoad("Q.y", "j") << ";\n";
    os << indent5 << cdbltype << " PQ_z = P_z - " << vinfo.DoubleLoad("Q.z", "j") << ";\n";


    os << indent5 << cdbltype << " R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;\n";
    os << "\n";
    os << indent5 << cdbltype << " alpha = PQalpha_mul * one_over_PQalpha_sum;   // alpha from MEST\n";

    if(hasoneoverp)
        os << indent5 << cdbltype << " one_over_p = " << "const_1 / P_alpha;\n";

    if(hasoneoverq)
        os << indent5 << cdbltype << " one_over_q = " << "const_1 / Q_alpha;\n";

    if(hasoneover2p)    
        os << indent5 << cdbltype << " one_over_2p = " << "one_half * one_over_p;  // gets multiplied in VRR\n";

    if(hasoneover2q)    
        os << indent5 << cdbltype << " one_over_2q = " << "one_half * one_over_q;  // gets multiplied in VRR\n";

    if(hasoneover2pq)
        os << indent5 << cdbltype << " one_over_2pq = " << "one_half * one_over_PQalpha_sum;\n";

    if(hasketvrr)
    {
        if(vrr_writer.HasVRR_K())
        {
            os << indent5 << vinfo.NewConstDoubleLoad("Q_PA_x", "Q.PA_x", "j") << ";\n";
            os << indent5 << vinfo.NewConstDoubleLoad("Q_PA_y", "Q.PA_y", "j") << ";\n";
            os << indent5 << vinfo.NewConstDoubleLoad("Q_PA_z", "Q.PA_z", "j") << ";\n";
        }
        else
        {
            os << indent5 << vinfo.NewConstDoubleLoad("Q_PB_x", "Q.PB_x", "j") << ";\n";
            os << indent5 << vinfo.NewConstDoubleLoad("Q_PB_y", "Q.PB_y", "j") << ";\n";
            os << indent5 << vinfo.NewConstDoubleLoad("Q_PB_z", "Q.PB_z", "j") << ";\n";
        }
    }

    if(hasbravrr)
    {
        os << "\n";
        os << indent5 << "// NOTE: Minus sign!\n";
        os << indent5 << cdbltype << " a_over_p =  -alpha * one_over_p;     // a/p from MEST\n";
        os << indent5 << cdbltype << " aop_PQ_x = a_over_p * PQ_x;\n"; 
        os << indent5 << cdbltype << " aop_PQ_y = a_over_p * PQ_y;\n"; 
        os << indent5 << cdbltype << " aop_PQ_z = a_over_p * PQ_z;\n"; 
    }

    if(hasketvrr)
    {
        os << "\n";
        os << indent5 << "// NOTE: Minus sign\n";
        os << indent5 << "// NOTE2: Plus sign taken care of on aoq_PQ_x!\n";
        os << indent5 << cdbltype << " a_over_q =  -alpha * one_over_q;     // a/q from MEST\n";
        os << indent5 << cdbltype << " aoq_PQ_x = -a_over_q * PQ_x;\n"; 
        os << indent5 << cdbltype << " aoq_PQ_y = -a_over_q * PQ_y;\n"; 
        os << indent5 << cdbltype << " aoq_PQ_z = -a_over_q * PQ_z;\n"; 

    }

    if(hasketet || hasbraet)
    {
        os << indent5 << vinfo.NewConstDoubleLoad("Q_bAB_x", "Q.bAB_x", "j") << ";\n";
        os << indent5 << vinfo.NewConstDoubleLoad("Q_bAB_y", "Q.bAB_y", "j") << ";\n";
        os << indent5 << vinfo.NewConstDoubleLoad("Q_bAB_z", "Q.bAB_z", "j") << ";\n";
    }

    if(hasketet)
    {
        os << "\n";
        os << indent5 << cdbltype << " p_over_q = P_alpha * one_over_q;\n";
        os << indent5 << cdbltype << " etfac_k[3] = {\n";
        os << indent6 << "-(P_bAB_x + Q_bAB_x) * one_over_q,\n";
        os << indent6 << "-(P_bAB_y + Q_bAB_y) * one_over_q,\n";
        os << indent6 << "-(P_bAB_z + Q_bAB_z) * one_over_q,\n";
        os << indent6 << "};\n";
        os << "\n";
    }

    if(hasbraet)
    {
        os << "\n";
        os << indent5 << cdbltype << " q_over_p = Q_alpha * one_over_p;\n";
        os << indent5 << cdbltype << " etfac_b[3] = {\n";
        os << indent6 << "-(P_bAB_x + Q_bAB_x) * one_over_p,\n";
        os << indent6 << "-(P_bAB_y + Q_bAB_y) * one_over_p,\n";
        os << indent6 << "-(P_bAB_z + Q_bAB_z) * one_over_p,\n";
        os << indent6 << "};\n";
        os << "\n";
    }

    os << "\n";
    os << "\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << indent5 << "// Boys function section\n";
    os << indent5 << "// Maximum v value: " << info.L() << "\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << indent5 << "// The paremeter to the boys function\n";
    os << indent5 << cdbltype << " F_x = R2 * alpha;\n";
    os << "\n";
    os << "\n";

    bg.WriteBoys(os, info);

    vrr_writer.WriteVRR(os, info);

    et_writer.WriteET(os, info);

    WriteAccumulation(os, info);

        
    os << "\n";
    os << indent4 << "}  // close loop over j\n";
    os << indent3 << "}  // close loop over i\n";

    os << indent3 << "\n";
    os << indent3 << "//Advance to the next batch\n";
    os << indent3 << "jstart = SIMINT_SIMD_ROUND(jend);\n";
    if(!hashrr)
        os << indent3 << "abcd += nshellbatch;\n";
    os << indent3 << "\n";


    hrr_writer.WriteHRR(os, info);

    os << "\n";

    os << indent3 << "cd += nshellbatch;\n";

    os << indent2 << "}   // close loop cdbatch\n";

    os << "\n";
    os << indent2 << "istart = iend;\n";

    os << indent2 << "// if this is the end of a batch in the bra part, skip the padding\n";
    os << indent2 << "if( ((ab+1) % SIMINT_NSHELL_SIMD) == 0)\n";
    os << indent3 << "istart = SIMINT_SIMD_ROUND(istart);\n";
    os << "\n";

    os << indent1 << "}  // close loop over ab\n";
    os << "\n";
    os << "\n";

    os << "\n";

    FreeContwork(os, info);

    os << indent1 << "return P.nshell12 * Q.nshell12;\n";
    os << "}\n";
    os << "\n";
}




