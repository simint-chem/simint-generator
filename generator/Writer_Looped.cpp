#include <sstream>
#include <iostream>
#include "generator/Classes.hpp"
#include "generator/Boys.hpp"
#include "generator/Helpers.hpp"
#include "generator/AlgorithmBase.hpp"

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

static const char * amchar = "spdfghijklmnoqrtuvwxyzabe";


// This is split out to allow for future optimization
static void Write_BoysFunction(std::ostream & os,
                               const BoysMap & bm,
                               int maxv)
{
    for(int i = 0; i <= maxv; i++)
        os << bm.at(i)->code_line() << "\n";
}


static std::string HRRVarString(const QAMList & am)
{
    std::stringstream ss;
    ss << "SQ_"  << am[0] << "_" << am[1] << "_" << am[2] << "_" << am[3];;
    return ss.str();
}


static std::string HRRBraStepVariable(const Doublet & d, const DoubletSetMap::value_type kets, bool istarget)
{
    // get the ket flags
    // They should all be the same, so just get the first one
    int ketflags = kets.second.begin()->flags;

    std::stringstream ss;

    if(d.flags & DOUBLET_HRRTOPLEVEL && ketflags & DOUBLET_HRRTOPLEVEL)
    {
        ss << HRRVarString({d.left.am(), d.right.am(), kets.first, 0}) 
           << "[abcd * " << d.ncart() * kets.second.size() 
                         << " + " << d.idx() << " * " << kets.second.size() << " + ni]";
    }

    else if(d.flags & DOUBLET_INITIAL & ketflags & DOUBLET_INITIAL) 
    {
        // we are calculating a final integral here
        // ie integral is in the form of ( a b | c 0 )
        ss << "integrals" << "[abcd * " << d.ncart() * kets.second.size() 
                                        << " + " << d.idx() << " * " << kets.second.size()
                                        << " + ni]";
    }

    else if(d.flags & DOUBLET_INITIAL)
    {
        // this is a 'bra target'
        // determine the shell quartet
        ss << HRRVarString({d.left.am(), d.right.am(), kets.first, 0})
           << "[abcd * " << d.ncart() * kets.second.size() 
           << " + " << d.idx() << " * " << kets.second.size() << " + ni]";
    }
    else
    {
        if(istarget)
            ss << "const double ";
        ss << "BRA_" << d.left.ijk[0] << d.left.ijk[1] << d.left.ijk[2]
           << "_" <<  d.right.ijk[0] << d.right.ijk[1] << d.right.ijk[2];
    }
    return ss.str();
}


static std::string HRRKetStepVariable(const Doublet & d, const DAMList & braam, bool istarget)
{
    std::stringstream ss;

    const int ncart_bra = NCART(braam[0]) * NCART(braam[1]);

    // bra should always be initial

    if(d.flags & DOUBLET_HRRTOPLEVEL)  // initial bra, hrr top level ket = a bra target
    {
        ss << HRRVarString({braam[0], braam[1], d.left.am(), d.right.am()}) 
           << "[abcd * " << ncart_bra * d.ncart() 
           << " + abi * " << ncart_bra << " + " << d.idx() << "]";
    }

    else if(d.flags & DOUBLET_INITIAL)
    {
        ss << "integrals" << "[abcd * " << ncart_bra * d.ncart() 
                                        << " + abi * " << d.ncart() << " + " << d.idx() << "]";
    }

    else
    {
        if(istarget)
            ss << "const double ";
        ss << "KET_" << d.left.ijk[0] << d.left.ijk[1] << d.left.ijk[2]
           << "_" <<  d.right.ijk[0] << d.right.ijk[1] << d.right.ijk[2];
    }
    return ss.str();
}


static std::string HRRBraStepString(const HRRDoubletStep & hrr, const DoubletSetMap::value_type kets)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "AB_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRBraStepVariable(hrr.target, kets, true);

    ss << " = ";
    ss << HRRBraStepVariable(hrr.src1, kets, false);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRBraStepVariable(hrr.src2, kets, false);
    ss << " );";

    return ss.str();
}

static std::string HRRKetStepString(const HRRDoubletStep & hrr, const DAMList & braam)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "CD_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRKetStepVariable(hrr.target, braam, true);

    ss << " = ";
    ss << HRRKetStepVariable(hrr.src1, braam, false);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRKetStepVariable(hrr.src2, braam, false);
    ss << " );";

    return ss.str();
}


void Writer_Looped(std::ostream & os,
                   const QAMList & am,
                   const std::string & nameappend,
                   const BoysMap & bm,
                   VRR_Algorithm_Base & vrralgo,
                   ET_Algorithm_Base & etalgo,
                   HRR_Algorithm_Base & hrralgo)
{
    //int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);
    int ncart_bra = NCART(am[0]) * NCART(am[1]);

    const int L = am[0] + am[1] + am[2] + am[3];
    const int maxv = L+1;


    // I need:
    // 0.) HRR Steps
    HRRBraKetStepList hrrsteps = hrralgo.Create_DoubletStepLists(am);

    // 1.) HRR Top level Doublets, sorted into their AM
    DoubletSetMap hrrtopbras, hrrtopkets;
    for(const auto & it : hrrsteps.first)
    {
        if(it.src1.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopbras[it.src1.am()].insert(it.src1);
        if(it.src2.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopbras[it.src2.am()].insert(it.src2);
    }

    for(const auto & it : hrrsteps.second)
    {
        if(it.src1.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopkets[it.src1.am()].insert(it.src1);
        if(it.src2.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopkets[it.src2.am()].insert(it.src2);
    }

    // 2.) ET steps
    //     with the HRR top level stuff as the initial targets
    QuartetSet etinit;
    for(const auto & it : hrrtopbras)
    for(const auto & it2 : hrrtopkets)
    {
        for(const auto & dit : it.second)
        for(const auto & dit2 : it2.second)
            etinit.insert({dit, dit2, 0, QUARTET_HRRTOPLEVEL});
    }

    ETStepList etsl = etalgo.Create_ETStepList(etinit);
        

    // 3.) VRR Steps
    //     This is #3 because we will need ET steps for #2
    //     For now, this is dependent on L
    GaussianSet gsreq{AllGaussiansForAM(L+1)};
    std::pair<VRRMap, VRRReqMap> vrrinfo = vrralgo.CreateAllMaps(gsreq);

    // these might be ( X s | or | X s )
    if(hrrtopbras.size() == 0)
    {
        GaussianSet gs = AllGaussiansForAM(am[0]);
        for(const auto & it : gs)
            hrrtopbras[it.am()].insert({DoubletType::BRA, it, {0,0,0}, DOUBLET_INITIAL | DOUBLET_HRRTOPLEVEL});
    }
    if(hrrtopkets.size() == 0)
    {
        GaussianSet gs = AllGaussiansForAM(am[2]);
        for(const auto & it : gs)
            hrrtopkets[it.am()].insert({DoubletType::KET, it, {0,0,0}, DOUBLET_INITIAL | DOUBLET_HRRTOPLEVEL});
    }


    std::cout << "VRR Top BrasS: " << hrrtopbras.size() << "\n";
    for(const auto & it : hrrtopbras)
        std::cout << " AM: " << it.first << "   Elements: " << it.second.size() << " / " << ((it.first+1)*(it.first+2)/2) << "\n";
    std::cout << "VRR Top Kets: " << hrrtopkets.size() << "\n";
    for(const auto & it : hrrtopkets)
        std::cout << " AM: " << it.first << "   Elements: " << it.second.size() << " / " << ((it.first+1)*(it.first+2)/2) << "\n";
    std::cout << "\n";

    //////////////////////////////////////////////////
    // set up the function to make it look pretty
    std::stringstream ss;
    ss << "int eri_" << nameappend << "_"
       << amchar[am[0]] << amchar[am[1]]
       << amchar[am[2]] << amchar[am[3]] << "(";

    std::string funcline = ss.str();
    std::string indent(funcline.length(), ' ');

    // start output to the file
    os << "#include <string.h>\n";
    os << "#include <math.h>\n";
    os << "\n";
    os << "#include \"vectorization.h\"\n";
    os << "#include \"constants.h\"\n";
    os << "#include \"eri/shell.h\"\n";
    os << "\n\n";
    os << funcline;
    os << "struct multishell_pair const P,\n";
    os << indent << "struct multishell_pair const Q,\n";
    os << indent << "double * const restrict integrals)\n";
    os << "{\n";
    os << "\n";
    os << "    ASSUME_ALIGN(P.x);\n";
    os << "    ASSUME_ALIGN(P.y);\n";
    os << "    ASSUME_ALIGN(P.z);\n";
    os << "    ASSUME_ALIGN(P.alpha);\n";
    os << "    ASSUME_ALIGN(P.prefac);\n";
    os << "    ASSUME_ALIGN(Q.x);\n";
    os << "    ASSUME_ALIGN(Q.y);\n";
    os << "    ASSUME_ALIGN(Q.z);\n";
    os << "    ASSUME_ALIGN(Q.alpha);\n";
    os << "    ASSUME_ALIGN(Q.prefac);\n";
    os << "    ASSUME_ALIGN(integrals)\n";
    os << "\n";
    os << "    const int nshell1234 = P.nshell12 * Q.nshell12;\n";
    os << "\n";
    os << "    memset(integrals, 0, nshell1234*sizeof(double));\n";
    os << "\n";
    os << "    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later\n";
    os << "    double AB_x[nshell1234];  double CD_x[nshell1234];\n";
    os << "    double AB_y[nshell1234];  double CD_y[nshell1234];\n";
    os << "    double AB_z[nshell1234];  double CD_z[nshell1234];\n";
    os << "\n";
    os << "    int ab, cd, abcd;\n";
    os << "    int i, j;\n";
    os << "\n";
    os << "    // Top level HRR requirements. Contracted integrals are accumulated here\n";

    // HRR top level stuff
    // each topbra/topket combo
    for(const auto & it : hrrtopbras)
        for(const auto & it2 : hrrtopkets)
            os << "    double " << HRRVarString({it.first, 0, it2.first, 0}) << "[" << it.second.size() << " * " << it2.second.size() << " * nshell1234];\n";

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
    os << "            const int cdstart = Q.primstart[cd];\n";
    os << "            const int cdend = Q.primend[cd];\n";
    os << "\n";
    os << "            // this should have been set/aligned in fill_multishell_pair or something else\n";
    os << "            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);\n";
    os << "\n";
    os << "            // Store for later\n";
    os << "            AB_x[abcd] = P.AB_x[ab];  CD_x[abcd] = P.AB_x[cd];\n";
    os << "            AB_y[abcd] = P.AB_y[ab];  CD_y[abcd] = P.AB_y[cd];\n";
    os << "            AB_z[abcd] = P.AB_z[ab];  CD_z[abcd] = P.AB_z[cd];\n";
    os << "\n";
    os << "            for(i = abstart; i < abend; ++i)\n";
    os << "            {\n";
    os << "                for(j = cdstart; j < cdend; ++j)\n";
    os << "                {\n";
    os << "\n";
    os << "                    // Holds the auxiliary integrals ( i 0 | 0 0 )^m in the primitive basis\n";
    os << "                    // with m as the slowest index\n";
    for(const auto & greq : vrrinfo.second)
    {
        // size of requirements for this AM
        //GaussianSet greq = vrrinfo.second.at(i);
        os << "                    // AM = " << greq.first << ": Needed from this AM: " << greq.second.size() << "\n";
        os << "                    double AUX_" << greq.first << "[" << (maxv-greq.first+1) << " * " << greq.second.size() << "];\n";
        os << "\n";
    }
    os << "\n\n";
    os << "                    const double PQalpha_mul = P.alpha[i] * Q.alpha[j];\n";
    os << "                    const double PQalpha_sum = P.alpha[i] + Q.alpha[j];\n";
    os << "\n";
    os << "                    const double pfac = TWO_PI_52 / (PQalpha_mul * sqrt(PQalpha_sum));\n";
    os << "\n";
    os << "                    /* construct R2 = (Px - Qx)**2 + (Py - Qy)**2 + (Pz -Qz)**2 */\n";
    os << "                    const double PQ_x = P.x[i] - Q.x[j];\n";
    os << "                    const double PQ_y = P.y[i] - Q.y[j];\n";
    os << "                    const double PQ_z = P.z[i] - Q.z[j];\n";
    os << "                    const double R2 = PQ_x*PQ_x + PQ_y*PQ_y + PQ_z*PQ_z;\n";
    os << "\n";
    os << "                    // collected prefactors\n";
    os << "                    const double allprefac =  pfac * P.prefac[i] * Q.prefac[j];\n";
    os << "\n";
    os << "                    // various factors\n";
    os << "                    const double alpha = PQalpha_mul/PQalpha_sum;   // alpha from MEST\n";
    os << "                    const double aoverp =  alpha / P.alpha[i];     // a/p from MEST\n";
    os << "                    const double one_over_2p = 1.0/(2.0 * P.alpha[i]);  // gets multiplied by i in VRR\n";
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Boys function section\n";
    os << "                    // Maximum v value: " << maxv << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // The paremeter to the boys function\n";
    os << "                    const double F_x = R2 * alpha;\n";
    os << "\n";
    os << "\n";

    Write_BoysFunction(os, bm, maxv);

    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Vertical recurrance\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "\n";
    os << "                    int idx = 0;\n";
    os << "\n";

    // i loops over am
    for(int i = 1; i <= L+1; i++)
    {
        // requirements for this am
        GaussianSet greq = vrrinfo.second.at(i);
        os << "                    // Forming AUX_" << i << "[" << (maxv-i+1) << " * " << greq.size() << "];\n";
        os << "                    // Needed from this AM:\n";
        for(const auto & it : greq)
            os << "                    //    " << it << "\n";
        os << "                    idx = 0;\n";
        os << "                    for(int m = 0; m < " << (maxv-i+1) << "; m++)  // loop over orders of boys function\n";
        os << "                    {\n";

        for(const auto & it : greq) // should be in order, since it's a set
        {
            // Get the stepping
            XYZStep step = vrrinfo.first.at(it);

            // and then the required gaussians
            Gaussian g1 = it.StepDown(step, 1);
            Gaussian g2 = it.StepDown(step, 2);

            // number of gaussians in the previous two AM
            int ng1 = -1;
            int ng2 = -1;

            int g1_idx = -1;
            int g2_idx = -1;

            // indices. This is a little messy, but whatever
            if(g1)
            { 
                auto g1_it = vrrinfo.second.at(i-1).find(g1); 
                g1_idx = std::distance(vrrinfo.second.at(i-1).begin(), g1_it);
                ng1 = vrrinfo.second.at(i-1).size();
            }
            if(g2)
            {
                auto g2_it = vrrinfo.second.at(i-2).find(g2); 
                g2_idx = std::distance(vrrinfo.second.at(i-2).begin(), g2_it);
                ng2 = vrrinfo.second.at(i-2).size();
            }

            // the value of i in the VRR eqn
            // = value of exponent of g1 in the position of the step
            int vrr_i = g1.ijk[XYZStepToIdx(step)];

            os << "                        //" << it <<  " : STEP: " << step << "\n";
//            os << "                        NEED: " << g1 << "  ,  " << g2 << "\n";
//            os << "                         IDX: " << g1_idx << "  ,  " << g2_idx << "\n";
//            os << "                       OUTOF: " << ng1 << "  ,  " << ng2 << "\n";
            os << "                        AUX_" << i << "[idx++] = P.PA_" << step << "[i] * ";

            if(g1)
                os << "AUX_" << (i-1) << "[m * " << ng1 << " + " << g1_idx 
                   << "] - aoverp * PQ_" << step << " * AUX_" << (i-1) << "[(m+1) * " << ng1 << " + " << g1_idx << "]";
            if(g2)
            {
                os << "\n";
                os << "                                     + " << vrr_i << " * one_over_2p * ( AUX_" << (i-2)
                                                              << "[m * " << ng2 << " +  " << g2_idx << "]"
                                                              << " - aoverp * AUX_" << (i-2) << "[(m+1) * " << ng2 << " + " << g2_idx << "] )";
            }
            os << ";\n\n"; 
        }

        os << "                    }\n";
        os << "\n\n";
    }
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Electron transfer\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                 }\n";
    os << "            }\n";
    os << "        }\n";
    os << "    }\n";
    os << "\n";
    os << "\n";

    os << "    //////////////////////////////////////////////\n";
    os << "    // Contracted integrals: Horizontal recurrance\n";
    os << "    // Bra part\n";
    os << "    // Steps: " << hrrsteps.first.size() << "\n";
    os << "    //////////////////////////////////////////////\n";
    os << "\n";
    if(hrrsteps.first.size() == 0)
        os << "    // Nothing to do.....\n";
    else
    {
        if(hrrsteps.second.size() == 0)
        {
            os << "    // No bra-specific targets. Should all end up in integrals\n";
        }
        else
        {
            os << "    // Bra targets\n";   // these are the top kets, with the final bra
            for(const auto & it : hrrtopkets)
                os << "    double " << HRRVarString({am[0], am[1], it.first, 0}) << "[" << ncart_bra << " * " << it.second.size() << " * nshell1234];\n";
        }

        os << "\n";    
        os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << "    {\n";

        for(const auto & it : hrrtopkets)
        {
            os << "        // form " << HRRVarString({am[0], am[1], it.first, 0}) << "\n";
            os << "        for(int ni = 0; ni < " << it.second.size() << "; ++ni)\n";
            os << "        {\n";
            for(const auto & hit : hrrsteps.first)
            {
                os << std::string(12, ' ') << "// " << hit << "\n";
                os << HRRBraStepString(hit, it) << "\n\n";
            }
            os << "        }\n";
            os << "\n";
        }
        os << "\n";
        os << "    }\n";
    }
    os << "\n";
    os << "\n";

    os << "    //////////////////////////////////////////////\n";
    os << "    // Contracted integrals: Horizontal recurrance\n";
    os << "    // Ket part\n";
    os << "    // Steps: " << hrrsteps.second.size() << "\n";
    os << "    // Forming final integrals\n";
    os << "    //////////////////////////////////////////////\n";
    os << "\n";
    if(hrrsteps.second.size() == 0)
        os << "    //Nothing to do.....\n";
    else
    { 
        DAMList braam{am[0], am[1]};

        os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << "    {\n";
        os << "        for(int abi = 0; abi < " << ncart_bra << "; ++abi)\n"; 
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
    os << "    return nshell1234;\n";
    os << "}\n";
    os << "\n";
}
