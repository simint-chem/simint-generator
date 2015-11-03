#include <sstream>
#include <iostream>

#include "generator/Classes.hpp"
#include "generator/Helpers.hpp"

#include "generator/FileWriter.hpp"
#include "generator/Boys.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/VRR_Writer.hpp"
#include "generator/ET_Writer.hpp"
#include "generator/HRR_Writer.hpp"
#include "generator/Ncart.hpp"

void WriteFile_Permute(std::ostream & os)
{
    QAM am = WriterInfo::FinalAM();
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

    std::stringstream ss;
    ss << "int eri_"
       << amchar[am[0]] << "_" << amchar[am[1]] << "_"
       << amchar[am[2]] << "_" << amchar[am[3]] << "(";

    std::string funcline = ss.str();
    std::string indent(funcline.length(), ' ');

    // start output to the file
    os << "#include \"eri/eri.h\"\n";
    os << "\n";

    // rest of includes are not needed

    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair P,\n";
    os << indent << "struct multishell_pair Q,\n";
    os << indent << "double * const restrict " << WriterInfo::ArrVarName(am) << ")\n";
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
                                   << amchar[tocall[3]] << "(P, Q, " << WriterInfo::ArrVarName(am) << ");\n"; 


    os << "}\n";
    os << "\n";
}


void WriteFile(std::ostream & os,
               const BoysGen & bg,
               const VRR_Writer & vrr_writer,
               const ET_Writer & et_writer,
               const HRR_Writer & hrr_writer)
{
    const QAM am = WriterInfo::FinalAM();
    int ncart = NCART(am);

    // some helper bools
    bool hashrr = WriterInfo::HasHRR();
    bool hasbrahrr = WriterInfo::HasBraHRR();
    bool haskethrr = WriterInfo::HasKetHRR();
    bool inline_hrr = (hashrr && WriterInfo::GetOption(OPTION_INLINEHRR) != 0);

    bool hasbravrr = vrr_writer.HasBraVRR();
    bool hasketvrr = vrr_writer.HasKetVRR();

    bool hasbraet = et_writer.HasBraET(); 
    bool hasketet = et_writer.HasKetET(); 

    bool hasoneoverp = (hasbravrr || hasbraet);
    bool hasoneoverq = (hasketvrr || hasketet);
    bool hasoneover2p = (hasbraet || (hasbravrr && (am[0]+am[1]) > 1)); 
    bool hasoneover2q = (hasketet || (hasketvrr && (am[2]+am[3]) > 1)); 
    bool hasoneover2pq = WriterInfo::GetOption(OPTION_NOET) && 
                         (am[0] + am[1] > 0) && (am[2] + am[3] > 0);


    // load this once here
    std::string dbltype = WriterInfo::DoubleType();
    std::string cdbltype = WriterInfo::ConstDoubleType();

    // we need a constant one for 1/x
    WriterInfo::AddIntConstant(1);

    std::stringstream ss;
    ss << "int eri_"
       << amchar[am[0]] << "_" << amchar[am[1]] << "_"
       << amchar[am[2]] << "_" << amchar[am[3]] << "(";

    std::string funcline = ss.str();
    std::string indent(funcline.length(), ' ');

    // start output to the file
    os << "#include <string.h>\n";
    os << "#include <math.h>\n";
    os << "\n";

    os << "#include \"eri/eri.h\"\n";
    os << "\n";
    WriterInfo::WriteIncludes(os);
    os << "\n";

    bg.WriteIncludes(os);

    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair const P,\n";
    os << indent << "struct multishell_pair const Q,\n";
    os << indent << "double * const restrict " << WriterInfo::ArrVarName(am) << ")\n";
    os << "{\n";
    os << "\n";

    // if we are manually using intrinsics, we don't need these assume lines
    // TODO: We won't need them for intrinsic calculations either, but HRR is still
    //       auto vectorized
    if(!WriterInfo::Scalar())
    {
        os << indent1 << "ASSUME_ALIGN_DBL(P.x);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.y);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.z);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.PA_x);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.PA_y);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.PA_z);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.PB_x);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.PB_y);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.PB_z);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.alpha);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(P.prefac);\n";
        os << "\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.x);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.y);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.z);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.PA_x);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.PA_y);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.PA_z);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.PB_x);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.PB_y);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.PB_z);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.alpha);\n";
        os << indent1 << "ASSUME_ALIGN_DBL(Q.prefac);\n";

        os << "\n";
        os << indent1 << "ASSUME_ALIGN_DBL(" << WriterInfo::ArrVarName(am) << ");\n";
        os << "\n";
        os << "\n";
    }

    // if there is no HRR, integrals are accumulated from inside the primitive loop
    // into the final integral array, so it must be zeroed first
    if(!hashrr)
        os << indent1 << "memset(" << WriterInfo::ArrVarName(am) << ", 0, P.nshell12 * Q.nshell12 * " << ncart << " * sizeof(double));\n";
    
    os << "\n";

    // abcd =  index within simd loop, real_abcd is the absolute
    // full abcd in terms of all the shells
    os << indent1 << "int ab, cd, cdbatch, abcd;\n";
    os << indent1 << "int istart, jstart;\n";

    os << indent1 << "int iprimcd, nprim_icd, np, icd;\n";

    if(hashrr)
        os << indent1 << "int real_abcd;\n";

    os << indent1 << "int i, j;\n";
    os << indent1 << "int n;\n";

    if(inline_hrr)
    {
        if(hasbrahrr)
            os << indent1 << "int iket;\n";
        if(haskethrr)
            os << indent1 << "int ibra;\n";
    }

    os << "\n";

    if(hashrr)
        WriterInfo::DeclareContwork(os);

    // need these factors sometimes
    if(hasoneover2p || hasoneover2q || hasoneover2pq)
        WriterInfo::AddNamedConstant("one_half", "0.5");

    bg.AddConstants();
    et_writer.AddConstants();
    vrr_writer.AddConstants();
    hrr_writer.AddConstants();
    WriterInfo::WriteConstants(os);

    


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
        WriterInfo::ZeroContWork(os);

        os << indent3 << "for(icd = 0; icd < nshellbatch; ++icd)\n";
        os << indent3<< "{\n"; 


        if(hasbrahrr)
        {
            os << indent4 << "AB_x[icd] = P.AB_x[ab];\n";
            os << indent4 << "AB_y[icd] = P.AB_y[ab];\n";
            os << indent4 << "AB_z[icd] = P.AB_z[ab];\n";
            os << "\n";
        }
        if(haskethrr)
        {
            os << indent4 << "CD_x[icd] = Q.AB_x[cd+icd];\n";
            os << indent4 << "CD_y[icd] = Q.AB_y[cd+icd];\n";
            os << indent4 << "CD_z[icd] = Q.AB_z[cd+icd];\n";
            os << "\n";
        }
        os << indent3<< "}\n"; 
    }


    os << "\n";

 

    if(hashrr)
        os << indent3 << "abcd = 0;\n";
    os << indent3 << "for(i = istart; i < iend; ++i)\n";
    os << indent3 << "{\n";
    os << "\n";

    os << indent4 << "icd = 0;\n";
    os << indent4 << "iprimcd = 0;\n";
    os << indent4 << "nprim_icd = Q.nprim12[cd];\n";

    vrr_writer.DeclarePrimPointers(os);
    os << "\n";
    et_writer.DeclarePrimPointers(os);
    os << "\n";

    os << indent4 << "// Load these one per loop over i\n";
    os << indent4 << WriterInfo::NewConstDoubleSet1("P_alpha", "P.alpha[i]") << ";\n";
    os << indent4 << WriterInfo::NewConstDoubleSet1("P_prefac", "P.prefac[i]") << ";\n";
    os << indent4 << WriterInfo::NewConstDoubleSet1("P_x", "P.x[i]") << ";\n";
    os << indent4 << WriterInfo::NewConstDoubleSet1("P_y", "P.y[i]") << ";\n";
    os << indent4 << WriterInfo::NewConstDoubleSet1("P_z", "P.z[i]") << ";\n";

    if(hasbravrr)
    {
        if(WriterInfo::HasVRR_I())
        {
            os << indent4 << WriterInfo::NewConstDoubleSet1("P_PA_x", "P.PA_x[i]") << ";\n";
            os << indent4 << WriterInfo::NewConstDoubleSet1("P_PA_y", "P.PA_y[i]") << ";\n";
            os << indent4 << WriterInfo::NewConstDoubleSet1("P_PA_z", "P.PA_z[i]") << ";\n";
        }
        else
        {
            os << indent4 << WriterInfo::NewConstDoubleSet1("P_PB_x", "P.PB_x[i]") << ";\n";
            os << indent4 << WriterInfo::NewConstDoubleSet1("P_PB_y", "P.PB_y[i]") << ";\n";
            os << indent4 << WriterInfo::NewConstDoubleSet1("P_PB_z", "P.PB_z[i]") << ";\n";
        }
    }

    if(hasketet)
    {
        os << indent4 << WriterInfo::NewConstDoubleSet1("P_bAB_x", "P.bAB_x[i]") << ";\n";
        os << indent4 << WriterInfo::NewConstDoubleSet1("P_bAB_y", "P.bAB_y[i]") << ";\n";
        os << indent4 << WriterInfo::NewConstDoubleSet1("P_bAB_z", "P.bAB_z[i]") << ";\n";
    }

    os << "\n";


    os << indent4 << "for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)\n";
    os << indent4 << "{\n";

    WriterInfo::WriteShellOffsets(os);

    os << "\n";

    vrr_writer.DeclarePrimArrays(os);
    et_writer.DeclarePrimArrays(os);

    os << indent5 << WriterInfo::NewConstDoubleLoad("Q_alpha", "Q.alpha", "j") << ";\n";
    os << indent5 << cdbltype << " PQalpha_mul = P_alpha * Q_alpha;\n";
    os << indent5 << cdbltype << " PQalpha_sum = P_alpha + Q_alpha;\n";
    os << indent5 << cdbltype << " one_over_PQalpha_sum = " << WriterInfo::NamedConstant("const_1") << " / PQalpha_sum;\n";
    os << "\n";
    os << "\n";
    os << indent5 << "/* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os << indent5 << cdbltype << " PQ_x = P_x - " << WriterInfo::DoubleLoad("Q.x", "j") << ";\n";
    os << indent5 << cdbltype << " PQ_y = P_y - " << WriterInfo::DoubleLoad("Q.y", "j") << ";\n";
    os << indent5 << cdbltype << " PQ_z = P_z - " << WriterInfo::DoubleLoad("Q.z", "j") << ";\n";


    os << indent5 << cdbltype << " R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;\n";
    os << "\n";
    os << indent5 << cdbltype << " alpha = PQalpha_mul * one_over_PQalpha_sum;   // alpha from MEST\n";

    if(hasoneoverp)
        os << indent5 << cdbltype << " one_over_p = " << WriterInfo::NamedConstant("const_1") << " / P_alpha;\n";

    if(hasoneoverq)
        os << indent5 << cdbltype << " one_over_q = " << WriterInfo::NamedConstant("const_1") << " / Q_alpha;\n";

    if(hasoneover2p)    
        os << indent5 << cdbltype << " one_over_2p = " << WriterInfo::NamedConstant("one_half") << " * one_over_p;  // gets multiplied in VRR\n";

    if(hasoneover2q)    
        os << indent5 << cdbltype << " one_over_2q = " << WriterInfo::NamedConstant("one_half") << " * one_over_q;  // gets multiplied in VRR\n";

    if(hasoneover2pq)
        os << indent5 << cdbltype << " one_over_2pq = " << WriterInfo::NamedConstant("one_half") << " * one_over_PQalpha_sum;\n";

    if(hasketvrr)
    {
        if(WriterInfo::HasVRR_K())
        {
            os << indent5 << WriterInfo::NewConstDoubleLoad("Q_PA_x", "Q.PA_x", "j") << ";\n";
            os << indent5 << WriterInfo::NewConstDoubleLoad("Q_PA_y", "Q.PA_y", "j") << ";\n";
            os << indent5 << WriterInfo::NewConstDoubleLoad("Q_PA_z", "Q.PA_z", "j") << ";\n";
        }
        else
        {
            os << indent5 << WriterInfo::NewConstDoubleLoad("Q_PB_x", "Q.PB_x", "j") << ";\n";
            os << indent5 << WriterInfo::NewConstDoubleLoad("Q_PB_y", "Q.PB_y", "j") << ";\n";
            os << indent5 << WriterInfo::NewConstDoubleLoad("Q_PB_z", "Q.PB_z", "j") << ";\n";
        }
    }

    if(hasbravrr)
    {
        os << indent5 << cdbltype << " a_over_p =  alpha * one_over_p;     // a/p from MEST\n";
        os << indent5 << cdbltype << " aop_PQ_x = a_over_p * PQ_x;\n"; 
        os << indent5 << cdbltype << " aop_PQ_y = a_over_p * PQ_y;\n"; 
        os << indent5 << cdbltype << " aop_PQ_z = a_over_p * PQ_z;\n"; 
    }

    if(hasketvrr)
    {
        os << indent5 << cdbltype << " a_over_q =  alpha * one_over_q;     // a/q from MEST\n";
        os << indent5 << cdbltype << " aoq_PQ_x = a_over_q * PQ_x;\n"; 
        os << indent5 << cdbltype << " aoq_PQ_y = a_over_q * PQ_y;\n"; 
        os << indent5 << cdbltype << " aoq_PQ_z = a_over_q * PQ_z;\n"; 

    }

    if(hasketet || hasbraet)
    {
        os << indent5 << WriterInfo::NewConstDoubleLoad("Q_bAB_x", "Q.bAB_x", "j") << ";\n";
        os << indent5 << WriterInfo::NewConstDoubleLoad("Q_bAB_y", "Q.bAB_y", "j") << ";\n";
        os << indent5 << WriterInfo::NewConstDoubleLoad("Q_bAB_z", "Q.bAB_z", "j") << ";\n";
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
    os << indent5 << "// Maximum v value: " << WriterInfo::L() << "\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << indent5 << "// The paremeter to the boys function\n";
    os << indent5 << cdbltype << " F_x = R2 * alpha;\n";
    os << "\n";
    os << "\n";

    bg.WriteBoys(os);

    vrr_writer.WriteVRR(os);

    et_writer.WriteET(os);

    WriterInfo::WriteAccumulation(os);

        
    os << "\n";
    os << indent4 << "}  // close loop over j\n";
    os << indent3 << "}  // close loop over i\n";

    os << indent3 << "\n";
    os << indent3 << "//Advance to the next batch\n";
    os << indent3 << "jstart = SIMINT_SIMD_ROUND(jend);\n";
    os << indent3 << "cd += nshellbatch;\n";
    if(!hashrr)
        os << indent3 << "abcd += nshellbatch;\n";
    os << indent3 << "\n";


    hrr_writer.WriteHRR(os);

    os << "\n";


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

    WriterInfo::FreeContwork(os);

    os << indent1 << "return P.nshell12 * Q.nshell12;\n";
    os << "}\n";
    os << "\n";
}

