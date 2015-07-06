#include <sstream>
#include <iostream>

#include "generator/Classes.hpp"
#include "generator/Helpers.hpp"

#include "generator/VRR_Writer.hpp"
#include "generator/ET_Writer.hpp"
#include "generator/HRR_Writer.hpp"

#include "generator/Boys.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/VRR_Writer.hpp"
#include "generator/ET_Writer.hpp"
#include "generator/HRR_Writer.hpp"


static void WriteFile_NotFlat(std::ostream & os,
                              const BoysGen & bg,
                              const VRR_Writer & vrr_writer,
                              const ET_Writer & et_writer,
                              const HRR_Writer & hrr_writer)
{
    const QAM am = WriterInfo::FinalAM();
    int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);

    // some helper bools
    bool hashrr = WriterInfo::HasHRR();
    bool hasbrahrr = WriterInfo::HasBraHRR();
    bool haskethrr = WriterInfo::HasKetHRR();
    bool inline_hrr = (hashrr && WriterInfo::GetOption(OPTION_INLINEHRR) != 0);

    bool hasvrr = WriterInfo::HasVRR();
    bool haset = WriterInfo::HasET();
    bool hasoneover2p = ((am[0] + am[1] + am[2] + am[3]) > 1);


    // load this once here
    std::string dbltype = WriterInfo::DoubleType();
    std::string cdbltype = WriterInfo::ConstDoubleType();

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

    os << "#include \"constants.h\"\n";
    os << "#include \"shell/shell.h\"\n";
    os << "\n";
    WriterInfo::WriteIncludes(os);
    os << "\n";

    bg.WriteIncludes(os);
    vrr_writer.WriteIncludes(os);
    hrr_writer.WriteIncludes(os);

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
        os << indent1 << "ASSUME_ALIGN(P.x);\n";
        os << indent1 << "ASSUME_ALIGN(P.y);\n";
        os << indent1 << "ASSUME_ALIGN(P.z);\n";
        os << indent1 << "ASSUME_ALIGN(P.PA_x);\n";
        os << indent1 << "ASSUME_ALIGN(P.PA_y);\n";
        os << indent1 << "ASSUME_ALIGN(P.PA_z);\n";
        os << indent1 << "ASSUME_ALIGN(P.bAB_x);\n";
        os << indent1 << "ASSUME_ALIGN(P.bAB_y);\n";
        os << indent1 << "ASSUME_ALIGN(P.bAB_z);\n";
        os << indent1 << "ASSUME_ALIGN(P.alpha);\n";
        os << indent1 << "ASSUME_ALIGN(P.prefac);\n";
        os << "\n";
        os << indent1 << "ASSUME_ALIGN(Q.x);\n";
        os << indent1 << "ASSUME_ALIGN(Q.y);\n";
        os << indent1 << "ASSUME_ALIGN(Q.z);\n";
        os << indent1 << "ASSUME_ALIGN(Q.PA_x);\n";
        os << indent1 << "ASSUME_ALIGN(Q.PA_y);\n";
        os << indent1 << "ASSUME_ALIGN(Q.PA_z);\n";
        os << indent1 << "ASSUME_ALIGN(Q.bAB_x);\n";
        os << indent1 << "ASSUME_ALIGN(Q.bAB_y);\n";
        os << indent1 << "ASSUME_ALIGN(Q.bAB_z);\n";
        os << indent1 << "ASSUME_ALIGN(Q.alpha);\n";
        os << indent1 << "ASSUME_ALIGN(Q.prefac);\n";

        os << "\n";
        os << indent1 << "ASSUME_ALIGN(" << WriterInfo::ArrVarName(am) << ");\n";
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
    os << indent1 << "int ab, cd, abcd;\n";

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

    // load constants into variables
    // (important only really for instrinsics)
    // include the factor of TWO_PI_52

    WriterInfo::AddConstant("const_two_pi_52", "TWO_PI_52");

    // need these factors if there is VRR or ET
    if(hasoneover2p || haset)
        WriterInfo::AddConstant("one_half", "0.5");
    if(hasvrr || haset)
        WriterInfo::AddConstant("const_1", "1.0");

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

    os << indent1 << "for(ab = 0; ab < P.nshell12; ++ab)\n";
    os << indent1 << "{\n";

    os << indent2 << "const int istart = P.primstart[ab];\n";
    os << indent2 << "const int iend = istart + P.nprim12[ab];\n";
    os << "\n";


    // if we are manually using intrinsics, we don't need these assume lines
    if(!WriterInfo::Scalar() && !WriterInfo::Intrinsics())
    {
        os << indent2 << "// this should have been set/aligned in fill_multishell_pair or something else\n";
        os << indent2 << "ASSUME(istart%SIMINT_SIMD_ALIGN_DBL == 0);\n";
        os << "\n";
    }


    if(hashrr)
    {
        os << indent2 << "// holds the main counter over the ket parts\n";
        os << indent2 << "cd = 0;\n";
        os << "\n";
        os << indent2 << "while(cd < Q.nshell12)\n";
        os << indent2 << "{\n";
        os << "\n";
        os << indent3 << "int cdstop = cd + SIMINT_NSHELL_SIMD;\n";
        os << indent3 << "cdstop = (cdstop > Q.nshell12 ? Q.nshell12 : cdstop);\n";
        os << "\n";
        os << indent3 << "const int nshell1234 = cdstop - cd;   // how many we are actually calcualting\n";
        os << "\n";
        
        WriterInfo::ZeroContWork(os);

        os << "\n";
        os << indent3 << "for(abcd = 0; abcd < nshell1234; ++cd, ++abcd)\n";
        os << indent3 << "{\n";
    }
    else
    {
        os << indent3 << "for(cd = 0; cd < Q.nshell12; ++cd, ++abcd)\n";
        os << indent3 << "{\n";
    }

    if(hasbrahrr)
    {
        os << indent4 << "AB_x[abcd] = P.AB_x[ab];\n";
        os << indent4 << "AB_y[abcd] = P.AB_y[ab];\n";
        os << indent4 << "AB_z[abcd] = P.AB_z[ab];\n";
        os << "\n";
    }
    if(haskethrr)
    {
        os << indent4 << "CD_x[abcd] = Q.AB_x[cd];\n";
        os << indent4 << "CD_y[abcd] = Q.AB_y[cd];\n";
        os << indent4 << "CD_z[abcd] = Q.AB_z[cd];\n";
        os << "\n";
    }

    os << indent4 << "const int jstart = Q.primstart[cd];\n";


    if(WriterInfo::Intrinsics())
        os << indent4 << "const int jend = Q.primend[cd];\n"; // we want the end, including the padding
    else
        os << indent4 << "const int jend = jstart + Q.nprim12[cd];\n";


    os << "\n";

    // if we are manually using intrinsics, we don't need these assume lines
    if(!WriterInfo::Scalar() && !WriterInfo::Intrinsics())
    {
        os << indent4 << "// this should have been set/aligned in fill_multishell_pair or something else\n";
        os << indent4 << "ASSUME(jstart%SIMINT_SIMD_ALIGN_DBL == 0);\n";
        os << "\n";
    }
 
    vrr_writer.DeclarePrimPointers(os);
    os << "\n";
    et_writer.DeclarePrimPointers(os);
    os << "\n";


    os << indent4 << "for(i = istart; i < iend; ++i)\n";
    os << indent4 << "{\n";
    os << "\n";
    os << indent5 << "// Load these one per loop over i\n";

    os << indent5 << WriterInfo::NewConstDoubleSet("P_alpha", "P.alpha[i]") << ";\n";
    os << indent5 << WriterInfo::NewConstDoubleSet("P_prefac", "P.prefac[i]") << ";\n";
    os << indent5 << WriterInfo::NewConstDoubleSet("P_x", "P.x[i]") << ";\n";
    os << indent5 << WriterInfo::NewConstDoubleSet("P_y", "P.y[i]") << ";\n";
    os << indent5 << WriterInfo::NewConstDoubleSet("P_z", "P.z[i]") << ";\n";

    if(hasvrr)
    {
        os << indent5 << WriterInfo::NewConstDoubleSet("P_PA_x", "P.PA_x[i]") << ";\n";
        os << indent5 << WriterInfo::NewConstDoubleSet("P_PA_y", "P.PA_y[i]") << ";\n";
        os << indent5 << WriterInfo::NewConstDoubleSet("P_PA_z", "P.PA_z[i]") << ";\n";
    }

    if(haset)
    {
        os << indent5 << WriterInfo::NewConstDoubleSet("P_bAB_x", "P.bAB_x[i]") << ";\n";
        os << indent5 << WriterInfo::NewConstDoubleSet("P_bAB_y", "P.bAB_y[i]") << ";\n";
        os << indent5 << WriterInfo::NewConstDoubleSet("P_bAB_z", "P.bAB_z[i]") << ";\n";
    }

    os << "\n";


    if(WriterInfo::Intrinsics())
    {
        os << indent5 << "#pragma novector\n";
        os << indent5 << "for(j = jstart; j < jend; j += SIMINT_SIMD_LEN)\n";
    }
    else
    {
        // commented out due to intel bug
        if(!WriterInfo::Scalar())
            os << indent5 << "//#pragma omp simd private(n)\n";

        os << indent5 << "for(j = jstart; j < jend; ++j)\n";
    }


    os << indent5 << "{\n";
    os << "\n";

    vrr_writer.DeclarePrimArrays(os);
    et_writer.DeclarePrimArrays(os);

    os << indent6 << WriterInfo::NewConstDoubleLoad("Q_alpha", "Q.alpha", "j") << ";\n";
    os << indent6 << cdbltype << " PQalpha_mul = P_alpha * Q_alpha;\n";
    os << indent6 << cdbltype << " PQalpha_sum = P_alpha + Q_alpha;\n";
    os << "\n";
    os << indent6 << cdbltype << " pfac = " << WriterInfo::DoubleConstant("const_two_pi_52") << " / (PQalpha_mul * " << WriterInfo::Sqrt("PQalpha_sum") << ");\n";
    os << "\n";
    os << indent6 << "/* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os << indent6 << cdbltype << " PQ_x = P_x - " << WriterInfo::DoubleLoad("Q.x", "j") << ";\n";
    os << indent6 << cdbltype << " PQ_y = P_y - " << WriterInfo::DoubleLoad("Q.y", "j") << ";\n";
    os << indent6 << cdbltype << " PQ_z = P_z - " << WriterInfo::DoubleLoad("Q.z", "j") << ";\n";


    os << indent6 << cdbltype << " R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;\n";
    os << "\n";
    os << indent6 << "// collected prefactors\n";
    os << indent6 << cdbltype << " allprefac =  pfac * P_prefac * " << WriterInfo::DoubleLoad("Q.prefac", "j") << ";\n";
    os << "\n";
    os << indent6 << "// various factors\n";
    os << indent6 << cdbltype << " alpha = PQalpha_mul/PQalpha_sum;   // alpha from MEST\n";

    if(hasvrr)
    {
        os << indent6 << "// for VRR\n";
        os << indent6 << cdbltype << " one_over_p = " << WriterInfo::DoubleConstant("const_1") << " / P_alpha;\n";
        os << indent6 << cdbltype << " a_over_p =  alpha * one_over_p;     // a/p from MEST\n";
        if(hasoneover2p)    
            os << indent6 << cdbltype << " one_over_2p = " << WriterInfo::DoubleConstant("one_half") << " * one_over_p;  // gets multiplied by i in VRR\n";

        os << "\n";
        os << indent6 << "// a_over_p * PQ_{xyz}\n";
        os << indent6 << cdbltype << " aop_PQ_x = a_over_p * PQ_x;\n"; 
        os << indent6 << cdbltype << " aop_PQ_y = a_over_p * PQ_y;\n"; 
        os << indent6 << cdbltype << " aop_PQ_z = a_over_p * PQ_z;\n"; 
        os << "\n";
    }

    if(haset)
    {
        os << indent6 << "// for electron transfer\n";
        os << indent6 << cdbltype << " one_over_q = " << WriterInfo::DoubleConstant("const_1") << " / Q_alpha;\n";
        os << indent6 << cdbltype << " one_over_2q = " << WriterInfo::DoubleConstant("one_half") << " * one_over_q;\n";
        os << indent6 << cdbltype << " p_over_q = P_alpha * one_over_q;\n";
        os << "\n";

        os << indent6 << cdbltype << " etfac[3] = {\n";
        os << indent7 << "-(P_bAB_x + " << WriterInfo::DoubleLoad("Q.bAB_x", "j") << ") * one_over_q,\n";
        os << indent7 << "-(P_bAB_y + " << WriterInfo::DoubleLoad("Q.bAB_y", "j") << ") * one_over_q,\n";
        os << indent7 << "-(P_bAB_z + " << WriterInfo::DoubleLoad("Q.bAB_z", "j") << ") * one_over_q,\n";
        os << indent7 << "};\n";
    }

    os << "\n";
    os << "\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << indent6 << "// Boys function section\n";
    os << indent6 << "// Maximum v value: " << WriterInfo::L() << "\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << indent6 << "// The paremeter to the boys function\n";
    os << indent6 << cdbltype << " F_x = R2 * alpha;\n";
    os << "\n";
    os << "\n";

    bg.WriteBoys(os);

    vrr_writer.WriteVRR(os);

    et_writer.WriteETInline(os);
        
    os << "\n";
    os << indent5 << "}  // close loop over j\n";
    os << indent4 << "}  // close loop over i\n";
    os << indent3 << "}\n";  // close loop over abcd or cd

    os << "\n";
    os << "\n";

    hrr_writer.WriteHRR(os);

    os << "\n";


    if(hashrr)
        os << indent2 << "}\n";   // close loop over ab or cd

    os << indent1 << "}  // close loop over ab\n";
    os << "\n";
    os << "\n";

    os << "\n";

    WriterInfo::FreeContwork(os);

    os << indent1 << "return P.nshell12 * Q.nshell12;\n";
    os << "}\n";
    os << "\n";
}



///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// MAIN ENTRY POINT
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
void WriteFile(std::ostream & os,
               const BoysGen & bg,
               const VRR_Writer & vrr_writer,
               const ET_Writer & et_writer,
               const HRR_Writer & hrr_writer)
{


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    // Create the function
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    WriteFile_NotFlat(os, bg, vrr_writer, et_writer, hrr_writer);
}

