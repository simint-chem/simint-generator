#include <sstream>
#include <iostream>

#include "generator/Writer.hpp"
#include "generator/Classes.hpp"
#include "generator/Helpers.hpp"
#include "generator/AlgorithmBase.hpp"

#include "generator/Boys.hpp"
#include "generator/WriterBase.hpp"
#include "generator/VRRWriter.hpp"
#include "generator/ETWriter.hpp"
#include "generator/HRRWriter.hpp"


static const char * amchar = "spdfghijklmnoqrtuvwxyzabe";

size_t memory_cont;          // memory required for contracted integral storage (bytes)



static void Writer_Looped_NotFlat(std::ostream & os,
                                 const QAMList & am,
                                 const std::string & nameappend,
                                 const OptionsMap & options,
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
    bool haset = et_writer.HasET();
    bool hasoneover2p = ((am[0] + am[1] + am[2] + am[3]) > 1);


    std::stringstream ss;
    ss << "int eri_" << nameappend << "_"
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

    std::vector<std::string> boysinc = bg.includes();
    for(const auto & it : boysinc)
        os << "#include \"" << it << "\"\n";

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

    if(hasvrr)
        os << "    int m;\n";
    if(hasvrr || haset)
        os << "    int n;\n";

    if(hasbrahrr)
        os << "    int iket;\n";
    if(haskethrr)
        os << "    int ibra;\n";

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

    os << bg.all_code_lines(base.L());

    vrr_writer.WriteVRR(os, base);

    et_writer.WriteET(os, base);
        
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

#if 0
static void Writer_Looped_Flat(std::ostream & os,
                               const QAMList & am,
                               const std::string & nameappend,
                               const OptionsMap & options,
                               const BoysGen & bg,
                               const std::pair<VRRMap, VRRReqMap> & vrrinfo,
                               const ETStepList & etsl,
                               const std::set<QAMList> & etint,
                               const HRRBraKetStepList & hrrsteps,
                               const DoubletSetMap & hrrtopbras,
                               const DoubletSetMap & hrrtopkets)
{
    int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);
    int ncart_bra = NCART(am[0]) * NCART(am[1]);
    const int L = am[0] + am[1] + am[2] + am[3];

    // some helper bools
    bool hasvrr = ( (vrrinfo.second.size() > 1) || (vrrinfo.second.size() == 1 && vrrinfo.second.begin()->first != 0) );
    bool haset = (etsl.size() > 0);
    bool hasoneover2p = ((am[0] + am[1] + am[2] + am[3]) > 1);


    std::stringstream ss;
    ss << "int eri_" << nameappend << "_"
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

    std::vector<std::string> boysinc = bg.includes();
    for(const auto & it : boysinc)
        os << "#include \"" << it << "\"\n";

    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair_flat const P,\n";
    os << indent << "struct multishell_pair_flat const Q,\n";
    os << indent << "double * const restrict " << ArrVarName(am) << ")\n";
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
    os << "    ASSUME_ALIGN(" << ArrVarName(am) << ");\n";
    os << "\n";
    os << "    const int nshell1234 = P.nshell12 * Q.nshell12;\n";
    os << "\n";

    // if there is no HRR, integrals are accumulated from inside the primitive loop
    // into the final integral array, so it must be zeroed first
    if(!hashrr)
        os << "    memset(" << ArrVarName(am) << ", 0, nshell1234*" << ncart << "*sizeof(double));\n";
    
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

    if(hasvrr)
        os << "    int m;\n";
    if(hasvrr || haset)
        os << "    int n;\n";

    if(hrrsteps.first.size() > 0)
        os << "    int iket;\n";
    if(hrrsteps.second.size() > 0)
        os << "    int ibra;\n";

    os << "\n";

    if(hashrr && contq.size() > 0)
    {
        os << "    // Workspace for contracted integrals\n";
        os << "    double * const contwork = malloc(nshell1234 * " << memory_cont << ");\n";
        os << "    memset(contwork, 0, nshell1234 * " << memory_cont << ");\n";
        os << "\n";
        os << "    // partition workspace into shells\n";

        size_t ptidx = 0;
        for(const auto & it : contq)
        {
            if(it != am)
            {
                os << "    double * const " << ArrVarName(it) << " = contwork + (nshell1234 * " << ptidx << ");\n";
                ptidx += NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]);
            }
        }
    }



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

    if(vrrinfo.second.size() > 0)
    {
        os << "                    // set up pointers to the contracted integrals - VRR\n";

        // pointers for accumulation in VRR
        for(const auto & it : vrrinfo.second)
        {
            int vam = it.first;
            if(IsContArray({vam, 0, 0, 0}))
                os << "                    double * const restrict PRIM_" << ArrVarName({vam, 0, 0, 0}) << " = " 
                   << ArrVarName({vam, 0, 0, 0}) << " + (abcd * " << NCART(vam) << ");\n";
        }
    }

    if(etint.size() > 0)
    {
        // pointers for accumulation in ET
        os << "                    // set up pointers to the contracted integrals - Electron Transfer\n";
        for(const auto & it : etint)
        {
            if(IsContArray(it))
                os << "                    double * const restrict PRIM_" << ArrVarName(it) << " = " 
                   << ArrVarName(it) << " + (abcd * " << (NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3])) << ");\n";
        }
    }

    os << "\n";
    os << "                    // Holds the auxiliary integrals ( i 0 | 0 0 )^m in the primitive basis\n";
    os << "                    // with m as the slowest index\n";

    for(const auto & greq : vrrinfo.second)
    {
        os << "                    // AM = " << greq.first << ": Needed from this AM: " << greq.second.size() << "\n";
        os << "                    double " << AuxName(greq.first) << "[" << (L-greq.first+1) << " * " << NCART(greq.first) << "];\n";
        os << "\n";
    }


    os << "\n\n";
    os << "                    // Holds temporary integrals for electron transfer\n";
    for(const auto & it : etint)
    {
        // only if these aren't from vrr
        if(it[1] > 0 || it[2] > 0 || it[3] > 0)
            os << "                    double AUX_" << ArrVarName(it) << "[" << NCART(it[0]) * NCART(it[2]) << "];\n";
    }


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
    os << "                    // Maximum v value: " << L << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // The paremeter to the boys function\n";
    os << "                    const double F_x = R2 * alpha;\n";
    os << "\n";
    os << "\n";

    os << bg.all_code_lines(L);

    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Vertical recurrance\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "\n";

    WriteVRRInfo(os, vrrinfo, L);
    os << "\n";

    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Electron transfer\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "\n";

    WriteETInfo(os, etsl, etint);
        
    os << "\n";
    os << "                 }\n";
    os << "            }\n";
    os << "\n";
    os << "\n";
    os << "\n";


    if(hrrsteps.first.size() > 0)
    {
        os << "    //////////////////////////////////////////////\n";
        os << "    // Contracted integrals: Horizontal recurrance\n";
        os << "    // Bra part\n";
        os << "    // Steps: " << hrrsteps.first.size() << "\n";
        os << "    //////////////////////////////////////////////\n";
        os << "\n";
        os << "    #pragma simd\n";
        os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << "    {\n";

        for(const auto & it : hrrtopkets)
        {
            os << "        // form " << ArrVarName({am[0], am[1], it.first, 0}) << "\n";
            os << "        for(iket = 0; iket < " << it.second.size() << "; ++iket)\n";
            os << "        {\n";
            for(const auto & hit : hrrsteps.first)
            {
                os << std::string(12, ' ') << "// " << hit << "\n";
                os << HRRBraStepString(hit, it.first) << "\n\n";
            }
            os << "        }\n";
            os << "\n";
        }
        os << "\n";
        os << "    }\n";
        os << "\n";
        os << "\n";
    }

    if(hrrsteps.second.size() > 0)
    {
        os << "    //////////////////////////////////////////////\n";
        os << "    // Contracted integrals: Horizontal recurrance\n";
        os << "    // Ket part\n";
        os << "    // Steps: " << hrrsteps.second.size() << "\n";
        os << "    //////////////////////////////////////////////\n";
        os << "\n";

        DAMList braam{am[0], am[1]};

        os << "    #pragma simd\n";
        os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << "    {\n";
        os << "        for(ibra = 0; ibra < " << ncart_bra << "; ++ibra)\n"; 
        os << "        {\n"; 

        for(const auto & hit : hrrsteps.second)
        {
            os << std::string(12, ' ') << "// " << hit << "\n";
            os << HRRKetStepString(hit, braam) << "\n\n";
        }

        os << "        }\n"; 
        os << "    }\n";
    }
    os << "\n";
    os << "\n";

    if(hashrr && contq.size() > 0)
    {
        os << "    // Free contracted work space\n";
        os << "    free(contwork);\n";
        os << "\n";
    }

    os << "    return nshell1234;\n";
    os << "}\n";
    os << "\n";
}
#endif




///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
// MAIN ENTRY POINT
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////
void Writer_Looped(std::ostream & os,
                   const QAMList & am,
                   const std::string & nameappend,
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
    base.SetContQ(hrr_writer);


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
    /*
    std::cout << "HRR Top Bras: " << hrrtopbras.size() << "\n";
    for(const auto & it : hrrtopbras)
        std::cout << " AM: " << it.first << " : Require: " << it.second.size() << " / " << NCART(it.first) << "\n";
    std::cout << "HRR Top Kets: " << hrrtopkets.size() << "\n";
    for(const auto & it : hrrtopkets)
        std::cout << " AM: " << it.first << " : Require:  " << it.second.size() << " / " << NCART(it.first) << "\n";
    std::cout << "\n";

    std::cout << "ET Requirements: " << etrm.size() << "\n";
    for(const auto & greq : etrm)
    {
        std::cout << "AM = " << greq.first << " : Require: " << greq.second.size() << " / " << NCART(greq.first) << "\n";
        for(const auto & it : greq.second)
            std::cout << "    " << it << "\n";
    }
    std::cout << "\n";

    std::cout << "VRR Requirements: " << vrrinfo.second.size() << "\n";
    for(const auto & greq : vrrinfo.second)
        std::cout << "AM = " << greq.first << " : Require: " << greq.second.size() << " / " << NCART(greq.first) << "\n";
    std::cout << "\n";
    */

    std::cout << "MEMORY (unaligned, per shell quartet): " << base.MemoryReq() << "\n";


    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    // Create the function
    //////////////////////////////////////////////////
    //////////////////////////////////////////////////
    /*
    if(options.count(OPTION_FLATPRIM) && options.at(OPTION_FLATPRIM) != 0)
    {
        Writer_Looped_Flat(os, am, nameappend, options, bg,
                           vrrinfo, etsl, etint,
                           hrrsteps, hrrtopbras, hrrtopkets);
    }
    else
    */
        Writer_Looped_NotFlat(os, am, nameappend, options, bg,
                              base, vrr_writer, et_writer, hrr_writer);
}

