#include <sstream>
#include <iostream>

#include "generator/Classes.hpp"
#include "generator/Helpers.hpp"
#include "generator/AlgorithmBase.hpp"

#include "generator/Boys.hpp"
#include "generator/WriterBase.hpp"
#include "generator/VRRWriter.hpp"
#include "generator/ETWriter.hpp"
#include "generator/HRRWriter.hpp"


size_t memory_cont;          // memory required for contracted integral storage (bytes)



static void WriteFile_NotFlat(std::ostream & os,
                              const QAMList & am,
                              const OptionsMap & options,
                              const std::string & prefix,
                              const BoysGen & bg,
                              const WriterBase & base,
                              const VRRWriter & vrr_writer,
                              const ETWriter & et_writer,
                              const HRRWriter & hrr_writer)
{
    int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);

    // some helper bools
    bool hashrr = hrr_writer.HasHRR();
    bool hasbrahrr = hrr_writer.HasBraHRR();
    bool haskethrr = hrr_writer.HasKetHRR();

    bool hasvrr = vrr_writer.HasVRR();
    bool hasvrr_m = hasvrr && (options.at(OPTION_INLINEVRR) > 0);
    bool haset = et_writer.HasET();
    bool hasoneover2p = ((am[0] + am[1] + am[2] + am[3]) > 1);


    std::stringstream ss;
    ss << "int eri_" << prefix << "_"
       << amchar[am[0]] << "_" << amchar[am[1]] << "_"
       << amchar[am[2]] << "_" << amchar[am[3]] << "(";

    std::string funcline = ss.str();
    std::string indent(funcline.length(), ' ');

    // start output to the file
    os << "#include <string.h>\n";
    os << "#include <math.h>\n";
    os << "\n";

    os << "#include \"vectorization.h\"\n";
    os << "#include \"constants.h\"\n";
    os << "#include \"eri/shell.h\"\n";

    bg.WriteIncludes(os);
    vrr_writer.WriteIncludes(os, base);
    hrr_writer.WriteIncludes(os, base);

    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair const P,\n";
    os << indent << "struct multishell_pair const Q,\n";
    os << indent << "double * const restrict " << base.ArrVarName(am) << ")\n";
    os << "{\n";
    os << "\n";
    os << "    ASSUME_ALIGN(P.x);\n";
    os << "    ASSUME_ALIGN(P.y);\n";
    os << "    ASSUME_ALIGN(P.z);\n";
    os << "    ASSUME_ALIGN(P.PA_x);\n";
    os << "    ASSUME_ALIGN(P.PA_y);\n";
    os << "    ASSUME_ALIGN(P.PA_z);\n";
    os << "    ASSUME_ALIGN(P.bAB_x);\n";
    os << "    ASSUME_ALIGN(P.bAB_y);\n";
    os << "    ASSUME_ALIGN(P.bAB_z);\n";
    os << "    ASSUME_ALIGN(P.alpha);\n";
    os << "    ASSUME_ALIGN(P.prefac);\n";
    os << "\n";
    os << "    ASSUME_ALIGN(Q.x);\n";
    os << "    ASSUME_ALIGN(Q.y);\n";
    os << "    ASSUME_ALIGN(Q.z);\n";
    os << "    ASSUME_ALIGN(Q.PA_x);\n";
    os << "    ASSUME_ALIGN(Q.PA_y);\n";
    os << "    ASSUME_ALIGN(Q.PA_z);\n";
    os << "    ASSUME_ALIGN(Q.bAB_x);\n";
    os << "    ASSUME_ALIGN(Q.bAB_y);\n";
    os << "    ASSUME_ALIGN(Q.bAB_z);\n";
    os << "    ASSUME_ALIGN(Q.alpha);\n";
    os << "    ASSUME_ALIGN(Q.prefac);\n";

    os << "\n";
    os << "    ASSUME_ALIGN(" << base.ArrVarName(am) << ");\n";
    os << "\n";
    os << "    const int nshell1234 = P.nshell12 * Q.nshell12;\n";
    os << "\n";

    // if there is no HRR, integrals are accumulated from inside the primitive loop
    // into the final integral array, so it must be zeroed first
    if(!hashrr)
        os << "    memset(" << base.ArrVarName(am) << ", 0, nshell1234*" << ncart << "*sizeof(double));\n";
    
    os << "\n";

    if(hashrr)
    {
        os << "    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later\n";
        os << "    double AB_x[nshell1234];  double CD_x[nshell1234];\n";
        os << "    double AB_y[nshell1234];  double CD_y[nshell1234];\n";
        os << "    double AB_z[nshell1234];  double CD_z[nshell1234];\n";
        os << "\n";
    }

    os << "    int ab, cd, abcd;\n";
    os << "    int i, j;\n";

    if(hasvrr_m)
        os << "    int m;\n";
    if(hasvrr || haset)
        os << "    int n;\n";

    if(base.GetOption(OPTION_INLINEHRR) != 0)
    {
        if(hasbrahrr)
            os << "    int iket;\n";
        if(haskethrr)
            os << "    int ibra;\n";
    }

    os << "\n";

    base.DeclareContwork(os);

    os << "\n\n";
    os << "    ////////////////////////////////////////\n";
    os << "    // Loop over shells and primitives\n";
    os << "    ////////////////////////////////////////\n";
    os << "    for(ab = 0, abcd = 0; ab < P.nshell12; ++ab)\n";
    os << "    {\n";
    os << "        const int abstart = P.primstart[ab];\n";
    os << "        const int abend = P.primend[ab];\n";
    os << "\n";
    os << "        // this should have been set/aligned in fill_multishell_pair or something else\n";
    os << "        ASSUME(abstart%SIMD_ALIGN_DBL == 0);\n";
    os << "\n";
    os << "        for(cd = 0; cd < Q.nshell12; ++cd, ++abcd)\n";
    os << "        {\n";

    vrr_writer.DeclarePointers(os, base);
    et_writer.DeclarePointers(os, base);

    os << "\n";
    os << "            const int cdstart = Q.primstart[cd];\n";
    os << "            const int cdend = Q.primend[cd];\n";
    os << "\n";
    os << "            // this should have been set/aligned in fill_multishell_pair or something else\n";
    os << "            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);\n";
    os << "\n";

    if(hashrr) 
    {
        os << "            // Store for later\n";
        os << "            AB_x[abcd] = P.AB_x[ab];  CD_x[abcd] = Q.AB_x[cd];\n";
        os << "            AB_y[abcd] = P.AB_y[ab];  CD_y[abcd] = Q.AB_y[cd];\n";
        os << "            AB_z[abcd] = P.AB_z[ab];  CD_z[abcd] = Q.AB_z[cd];\n";
        os << "\n";
    }

    os << "            for(i = abstart; i < abend; ++i)\n";
    os << "            {\n";
    os << "\n";
    os << "                // Load these one per loop over i\n";
    os << "                const double P_alpha = P.alpha[i];\n";
    os << "                const double P_prefac = P.prefac[i];\n";
    os << "                const double P_x = P.x[i];\n";
    os << "                const double P_y = P.y[i];\n";
    os << "                const double P_z = P.z[i];\n";

    if(hasvrr)
    {
        os << "                const double P_PA_x = P.PA_x[i];\n";
        os << "                const double P_PA_y = P.PA_y[i];\n";
        os << "                const double P_PA_z = P.PA_z[i];\n";
    }

    if(haset)
    {
        os << "                const double P_bAB_x = P.bAB_x[i];\n";
        os << "                const double P_bAB_y = P.bAB_y[i];\n";
        os << "                const double P_bAB_z = P.bAB_z[i];\n";
    }

    os << "\n";
    os << "                //#pragma simd\n";
    os << "                for(j = cdstart; j < cdend; ++j)\n";
    os << "                {\n";
    os << "\n";

    vrr_writer.DeclareAuxArrays(os, base);
    et_writer.DeclareAuxArrays(os, base);

    os << "                    const double PQalpha_mul = P_alpha * Q.alpha[j];\n";
    os << "                    const double PQalpha_sum = P_alpha + Q.alpha[j];\n";
    os << "\n";
    os << "                    const double pfac = TWO_PI_52 / (PQalpha_mul * sqrt(PQalpha_sum));\n";
    os << "\n";
    os << "                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os << "                    const double PQ_x = P_x - Q.x[j];\n";
    os << "                    const double PQ_y = P_y - Q.y[j];\n";
    os << "                    const double PQ_z = P_z - Q.z[j];\n";
    os << "                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;\n";
    os << "\n";
    os << "                    // collected prefactors\n";
    os << "                    const double allprefac =  pfac * P_prefac * Q.prefac[j];\n";
    os << "\n";
    os << "                    // various factors\n";
    os << "                    const double alpha = PQalpha_mul/PQalpha_sum;   // alpha from MEST\n";

    if(hasvrr)
    {
        os << "                    // for VRR\n";
        os << "                    const double one_over_p = 1.0 / P_alpha;\n";
        os << "                    const double a_over_p =  alpha * one_over_p;     // a/p from MEST\n";
        if(hasoneover2p)
            os << "                    const double one_over_2p = 0.5 * one_over_p;  // gets multiplied by i in VRR\n";

        os << "\n";
        os << "                    // a_over_p * PQ_{xyz}\n";
        os << "                    const double aop_PQ_x = a_over_p * PQ_x;\n"; 
        os << "                    const double aop_PQ_y = a_over_p * PQ_y;\n"; 
        os << "                    const double aop_PQ_z = a_over_p * PQ_z;\n"; 
        os << "\n";
    }

    if(haset)
    {
        os << "                    // for electron transfer\n";
        os << "                    const double one_over_q = 1.0 / Q.alpha[j];\n";
        os << "                    const double one_over_2q = 0.5 * one_over_q;\n";
        os << "                    const double p_over_q = P_alpha * one_over_q;\n";
        os << "\n";

        os << "                    const double etfac[3] = {\n";
        os << "                                             -(P_bAB_x + Q.bAB_x[j]) * one_over_q,\n";
        os << "                                             -(P_bAB_y + Q.bAB_y[j]) * one_over_q,\n";
        os << "                                             -(P_bAB_z + Q.bAB_z[j]) * one_over_q,\n";
        os << "                                            };\n";
    }

    os << "\n";
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Boys function section\n";
    os << "                    // Maximum v value: " << base.L() << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // The paremeter to the boys function\n";
    os << "                    const double F_x = R2 * alpha;\n";
    os << "\n";
    os << "\n";

    bg.WriteBoys(os, base);

    vrr_writer.WriteVRR(os, base);

    et_writer.WriteETInline(os, base);
        
    os << "\n";
    os << "                 }\n";
    os << "            }\n";
    os << "        }\n";
    os << "    }\n";
    os << "\n";
    os << "\n";
    os << "\n";

    hrr_writer.WriteHRR(os, base);

    base.FreeContwork(os);

    os << "    return nshell1234;\n";
    os << "}\n";
    os << "\n";
}



static void WriteFile_Flat(std::ostream & os,
                           const QAMList & am,
                           const OptionsMap & options,
                           const std::string & prefix,
                           const BoysGen & bg,
                           const WriterBase & base,
                           const VRRWriter & vrr_writer,
                           const ETWriter & et_writer,
                           const HRRWriter & hrr_writer)
{
    int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);


    // some helper bools
    bool hashrr = hrr_writer.HasHRR();
    bool hasbrahrr = hrr_writer.HasBraHRR();
    bool haskethrr = hrr_writer.HasKetHRR();

    bool hasvrr = vrr_writer.HasVRR();
    bool hasvrr_m = hasvrr && (options.at(OPTION_INLINEVRR) > 0);
    bool haset = et_writer.HasET();
    bool hasoneover2p = ((am[0] + am[1] + am[2] + am[3]) > 1);


    std::stringstream ss;
    ss << "int eri_" << prefix << "_"
       << amchar[am[0]] << "_" << amchar[am[1]] << "_"
       << amchar[am[2]] << "_" << amchar[am[3]] << "(";

    std::string funcline = ss.str();
    std::string indent(funcline.length(), ' ');

    // start output to the file
    os << "#include <string.h>\n";
    os << "#include <math.h>\n";
    os << "\n";

    os << "#include \"vectorization.h\"\n";
    os << "#include \"constants.h\"\n";
    os << "#include \"eri/shell.h\"\n";

    bg.WriteIncludes(os);
    vrr_writer.WriteIncludes(os, base);
    hrr_writer.WriteIncludes(os, base);

    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair_flat const P,\n";
    os << indent << "struct multishell_pair_flat const Q,\n";
    os << indent << "double * const restrict " << base.ArrVarName(am) << ")\n";
    os << "{\n";
    os << "\n";
    os << "    ASSUME_ALIGN(P.x);\n";
    os << "    ASSUME_ALIGN(P.y);\n";
    os << "    ASSUME_ALIGN(P.z);\n";
    os << "    ASSUME_ALIGN(P.PA_x);\n";
    os << "    ASSUME_ALIGN(P.PA_y);\n";
    os << "    ASSUME_ALIGN(P.PA_z);\n";
    os << "    ASSUME_ALIGN(P.bAB_x);\n";
    os << "    ASSUME_ALIGN(P.bAB_y);\n";
    os << "    ASSUME_ALIGN(P.bAB_z);\n";
    os << "    ASSUME_ALIGN(P.alpha);\n";
    os << "    ASSUME_ALIGN(P.prefac);\n";
    os << "    ASSUME_ALIGN(P.shellidx);\n";
    os << "\n";
    os << "    ASSUME_ALIGN(Q.x);\n";
    os << "    ASSUME_ALIGN(Q.y);\n";
    os << "    ASSUME_ALIGN(Q.z);\n";
    os << "    ASSUME_ALIGN(Q.PA_x);\n";
    os << "    ASSUME_ALIGN(Q.PA_y);\n";
    os << "    ASSUME_ALIGN(Q.PA_z);\n";
    os << "    ASSUME_ALIGN(Q.bAB_x);\n";
    os << "    ASSUME_ALIGN(Q.bAB_y);\n";
    os << "    ASSUME_ALIGN(Q.bAB_z);\n";
    os << "    ASSUME_ALIGN(Q.alpha);\n";
    os << "    ASSUME_ALIGN(Q.prefac);\n";
    os << "    ASSUME_ALIGN(Q.shellidx);\n";
    os << "\n";
    os << "    ASSUME_ALIGN(" << base.ArrVarName(am) << ");\n";
    os << "\n";
    os << "    const int nshell1234 = P.nshell12 * Q.nshell12;\n";
    os << "\n";

    // if there is no HRR, integrals are accumulated from inside the primitive loop
    // into the final integral array, so it must be zeroed first
    if(!hashrr)
        os << "    memset(" << base.ArrVarName(am) << ", 0, nshell1234*" << ncart << "*sizeof(double));\n";
    
    os << "\n";

    if(hashrr)
    {
        os << "    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later\n";
        os << "    double AB_x[nshell1234];  double CD_x[nshell1234];\n";
        os << "    double AB_y[nshell1234];  double CD_y[nshell1234];\n";
        os << "    double AB_z[nshell1234];  double CD_z[nshell1234];\n";
        os << "\n";
    }

    if(hashrr)
        os << "    int abcd;\n";

    os << "    int i, j;\n";

    if(hasvrr_m)
        os << "    int m;\n";
    if(hasvrr || haset)
        os << "    int n;\n";

    if(base.GetOption(OPTION_INLINEHRR) != 0)
    {
        if(hasbrahrr)
            os << "    int iket;\n";
        if(haskethrr)
            os << "    int ibra;\n";
    }

    os << "\n";

    // save these for later
    if(hashrr)
    {
        os << "\n";
        os << "    // Store for later\n";
        os << "    abcd = 0;\n";
        os << "    for(i = 0; i < P.nshell12; i++)\n";
        os << "    for(j = 0; j < Q.nshell12; j++)\n";
        os << "    {\n";
        os << "        AB_x[abcd] = P.AB_x[i];  CD_x[abcd] = Q.AB_x[j];\n";
        os << "        AB_y[abcd] = P.AB_y[i];  CD_y[abcd] = Q.AB_y[j];\n";
        os << "        AB_z[abcd] = P.AB_z[i];  CD_z[abcd] = Q.AB_z[j];\n";
        os << "        abcd++;\n";
        os << "    }\n";
        os << "\n";
    }

    base.DeclareContwork(os);

    os << "\n\n";
    os << "    ////////////////////////////////////////\n";
    os << "    // Loop over shells and primitives\n";
    os << "    ////////////////////////////////////////\n";

    os << "\n";

    os << "            for(i = 0; i < P.nprim; ++i)\n";
    os << "            {\n";
    os << "                // determine part of the shell quartet index\n";
    os << "                const int Pshellidx = P.shellidx[i] * Q.nshell12;\n";
    os << "\n";
    os << "                // Load these one per loop over i\n";
    os << "                const double P_alpha = P.alpha[i];\n";
    os << "                const double P_prefac = P.prefac[i];\n";
    os << "                const double P_x = P.x[i];\n";
    os << "                const double P_y = P.y[i];\n";
    os << "                const double P_z = P.z[i];\n";

    if(hasvrr)
    {
        os << "                const double P_PA_x = P.PA_x[i];\n";
        os << "                const double P_PA_y = P.PA_y[i];\n";
        os << "                const double P_PA_z = P.PA_z[i];\n";
    }

    if(haset)
    {
        os << "                const double P_bAB_x = P.bAB_x[i];\n";
        os << "                const double P_bAB_y = P.bAB_y[i];\n";
        os << "                const double P_bAB_z = P.bAB_z[i];\n";
    }

    os << "\n";
    os << "                //#pragma simd\n";
    os << "                for(j = 0; j < Q.nprim; ++j)\n";
    os << "                {\n";

    os << "                    // Contracted shell quartet index\n";
    os << "                    const int abcd = Pshellidx + Q.shellidx[j];\n";
    os << "\n";

    vrr_writer.DeclarePointers(os, base);
    et_writer.DeclarePointers(os, base);

    os << "\n\n";

    vrr_writer.DeclareAuxArrays(os, base);
    et_writer.DeclareAuxArrays(os, base);

    os << "\n\n";
    os << "                    const double PQalpha_mul = P_alpha * Q.alpha[j];\n";
    os << "                    const double PQalpha_sum = P_alpha + Q.alpha[j];\n";
    os << "\n";
    os << "                    const double pfac = TWO_PI_52 / (PQalpha_mul * sqrt(PQalpha_sum));\n";
    os << "\n";
    os << "                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os << "                    const double PQ_x = P_x - Q.x[j];\n";
    os << "                    const double PQ_y = P_y - Q.y[j];\n";
    os << "                    const double PQ_z = P_z - Q.z[j];\n";
    os << "                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;\n";
    os << "\n";
    os << "                    // collected prefactors\n";
    os << "                    const double allprefac =  pfac * P_prefac * Q.prefac[j];\n";
    os << "\n";
    os << "                    // various factors\n";
    os << "                    const double alpha = PQalpha_mul/PQalpha_sum;   // alpha from MEST\n";

    if(hasvrr)
    {
        os << "                    // for VRR\n";
        os << "                    const double one_over_p = 1.0 / P_alpha;\n";
        os << "                    const double a_over_p =  alpha * one_over_p;     // a/p from MEST\n";
        if(hasoneover2p)
            os << "                    const double one_over_2p = 0.5 * one_over_p;  // gets multiplied by i in VRR\n";

        os << "\n";
        os << "                    // a_over_p * PQ_{xyz}\n";
        os << "                    const double aop_PQ_x = a_over_p * PQ_x;\n"; 
        os << "                    const double aop_PQ_y = a_over_p * PQ_y;\n"; 
        os << "                    const double aop_PQ_z = a_over_p * PQ_z;\n"; 
        os << "\n";
    }

    if(haset)
    {
        os << "                    // for electron transfer\n";
        os << "                    const double one_over_q = 1.0 / Q.alpha[j];\n";
        os << "                    const double one_over_2q = 0.5 * one_over_q;\n";
        os << "                    const double p_over_q = P_alpha * one_over_q;\n";
        os << "\n";

        os << "                    const double etfac[3] = {\n";
        os << "                                             -(P_bAB_x + Q.bAB_x[j]) * one_over_q,\n";
        os << "                                             -(P_bAB_y + Q.bAB_y[j]) * one_over_q,\n";
        os << "                                             -(P_bAB_z + Q.bAB_z[j]) * one_over_q,\n";
        os << "                                            };\n";
    }

    os << "\n";
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Boys function section\n";
    os << "                    // Maximum v value: " << base.L() << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // The paremeter to the boys function\n";
    os << "                    const double F_x = R2 * alpha;\n";
    os << "\n";
    os << "\n";

    bg.WriteBoys(os, base);

    vrr_writer.WriteVRR(os, base);

    et_writer.WriteETInline(os, base);

        
    os << "\n";
    os << "                 }\n";
    os << "            }\n";
    os << "\n";
    os << "\n";
    os << "\n";


    hrr_writer.WriteHRR(os, base);

    base.FreeContwork(os);

    os << "    return nshell1234;\n";
    os << "}\n";
    os << "\n";
}




///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// MAIN ENTRY POINT
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
void WriteFile(std::ostream & os,
               const QAMList & am,
               const std::string & prefix,
               const OptionsMap & options,
               const BoysGen & bg,
               VRR_Algorithm_Base & vrralgo,
               ET_Algorithm_Base & etalgo,
               HRR_Algorithm_Base & hrralgo)
{
    // set of contracted quartets
    QAMListSet contq;

    // add the am list to the contracted info
    contq.insert(am);

    // Base writer information
    WriterBase base(options, am);


    // Working backwards, I need:
    // 1.) HRR Steps
    HRRBraKetStepList hrrsteps = hrralgo.Create_DoubletStepLists(am);
    HRRWriter hrr_writer(hrrsteps, am);

    // set the contracted quartets
    base.SetContQ(hrr_writer.TopQuartets(), hrr_writer.TopKets());


    // 2.) ET steps
    //     with the HRR top level stuff as the initial targets
    QuartetSet etinit = hrr_writer.TopQuartets();
    ETStepList etsl = etalgo.Create_ETStepList(etinit);
    ETWriter et_writer(etsl);


    // 3.) VRR Steps
    // requirements for vrr are the elements of etrm
    ETReqMap vreq = et_writer.ETRMap();


    // and also any elements from top bra/kets in the form ( X 0 | 0 0 )
    for(const auto & it : hrr_writer.TopBras())
    for(const auto & it2 : hrr_writer.TopKets())
    {
        for(const auto & dit : it.second)
        for(const auto & dit2 : it2.second)
        {
            if(dit.right.am() == 0 && dit2.left.am() == 0 && dit2.right.am() == 0)
                vreq[dit.left.am()].insert(dit.left);
        }
    }

    std::pair<VRRMap, VRRReqMap> vrrinfo = vrralgo.CreateAllMaps(vreq);
    VRRWriter vrr_writer(vrrinfo.first, vrrinfo.second);


    ///////////////////////////////////////////////
    // Done with prerequisites
    ///////////////////////////////////////////////


    // print out some info
    std::cout << "MEMORY (per shell quartet): " << base.MemoryReq() << "\n";


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    // Create the function
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    if(options.count(OPTION_FLATPRIM) && options.at(OPTION_FLATPRIM) != 0)
    {
        WriteFile_Flat(os, am, options, prefix, bg,
                       base, vrr_writer, et_writer, hrr_writer);
    }
    else
        WriteFile_NotFlat(os, am, options, prefix, bg,
                          base, vrr_writer, et_writer, hrr_writer);
}

