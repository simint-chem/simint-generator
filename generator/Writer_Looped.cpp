#include <sstream>
#include <iostream>
#include "generator/Classes.hpp"
#include "generator/Boys.hpp"
#include "generator/Helpers.hpp"
#include "generator/AlgorithmBase.hpp"

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

static const char * amchar = "spdfghijklmnoqrtuvwxyzabe";

// contains all variables stored in arrays, rather than local 'double'
std::set<QAMList> arrinfo;
std::set<QAMList> continfo;  // is stored contracted

static bool IsArray(const QAMList & am)
{
    auto it = arrinfo.find(am);
    return (it != arrinfo.end());
}

static bool IsContArray(const QAMList & am)
{
    return continfo.count(am);
}

static std::string ArrVarName(const QAMList & am)
{
    std::stringstream ss;
    ss << "S_"  << am[0] << "_" << am[1] << "_" << am[2] << "_" << am[3];
    return ss.str();
}


static std::string VarName(const QAMList & am,
                           const std::string & idx,
                           const std::string & contidx,
                           const std::string & singleid,
                           bool istarget)
{
    std::stringstream ss;

    if(IsArray(am))
    {
        ss << ArrVarName(am) << "[";
        if(IsContArray(am))
            ss << contidx;
        else
            ss << idx;
        ss << "]"; 
    }
    else
    {
        if(istarget)
            ss << "const double ";
        ss << "Q_" << am[0] << "_" << am[1] << "_" << am[2] << "_" << am[3] << "__" << singleid;
    }

    return ss.str();
}

static std::string AuxName(int i)
{
    return std::string("AUX_") + ArrVarName({i, 0, 0, 0});
}




static std::string HRRBraStepArrVar(const Doublet & d, int ketam, bool istarget)
{
    int ncart_ket = NCART(ketam); // all am is on the left part of the ket

    // index to the array (always contracted)
    std::stringstream ssidx;
    ssidx << "abcd * " << d.ncart() * ncart_ket << " + " << d.idx() << " * " << ncart_ket << " + ni";

    std::stringstream ss;
    ss << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2]  << "_"
       << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2] << "_"
       << ketam;

    return VarName({d.left.am(), d.right.am(), ketam, 0},
                    ssidx.str(),
                    ssidx.str(),
                    ss.str(),
                    istarget);
    return ss.str();
}


static std::string HRRKetStepArrVar(const Doublet & d, const DAMList & braam, bool istarget)
{
    const int ncart_bra = NCART(braam[0]) * NCART(braam[1]);

    // index to the array (always contracted)
    std::stringstream ssidx;
    ssidx << "abcd * " << ncart_bra * d.ncart() << " + abi * " << ncart_bra << " + " << d.idx(); 

    std::stringstream ss;
    ss << braam[0] << "_" << braam[1] << "_" 
       << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2]  << "_"
       << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2];

    return VarName({braam[0], braam[1], d.left.am(), d.right.am()},
                    ssidx.str(),
                    ssidx.str(),
                    ss.str(),
                    istarget);

    return ss.str();
}


static std::string ETStepVar(const Quartet & q, bool istarget)
{
    // index to the array
    std::stringstream sscontidx;
    sscontidx << "abcd * " << q.ncart() << " + " << q.idx();

    std::stringstream ssidx;
    ssidx << q.idx();

    std::stringstream ss;
    ss << q.bra.left.ijk[0]  << "_" << q.bra.left.ijk[1]  << "_" << q.bra.left.ijk[2]  << "_"
       << q.bra.right.ijk[0]  << "_" << q.bra.right.ijk[1]  << "_" << q.bra.right.ijk[2]  << "_"
       << q.ket.left.ijk[0]  << "_" << q.ket.left.ijk[1]  << "_" << q.ket.left.ijk[2]  << "_"
       << q.ket.right.ijk[0]  << "_" << q.ket.right.ijk[1]  << "_" << q.ket.right.ijk[2];

    // prefix this if it should be a result from VRR
    std::string prefix;

    // always use 'uncontracted' stuff if not target
    if(!istarget || !IsContArray(q.amlist()))
        prefix = "AUX_";

    if(istarget) 
        return prefix + VarName(q.amlist(),
                                ssidx.str(),
                                sscontidx.str(),
                                ss.str(),
                                istarget);
    else
        return prefix + VarName(q.amlist(),
                                ssidx.str(),
                                ssidx.str(),  // force uncontracted
                                ss.str(),
                                istarget);
        
}



static std::string HRRBraStepString(const HRRDoubletStep & hrr, int ketam)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "AB_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRBraStepArrVar(hrr.target, ketam, true);

    ss << " = ";
    ss << HRRBraStepArrVar(hrr.src1, ketam, false);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRBraStepArrVar(hrr.src2, ketam, false);
    ss << " );";

    return ss.str();
}

static std::string HRRKetStepString(const HRRDoubletStep & hrr, const DAMList & braam)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "CD_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRKetStepArrVar(hrr.target, braam, true);

    ss << " = ";
    ss << HRRKetStepArrVar(hrr.src1, braam, true);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRKetStepArrVar(hrr.src2, braam, true);
    ss << " );";

    return ss.str();
}


static std::string ETStepString(const ETStep & et)
{
    int stepidx = XYZStepToIdx(et.xyz);
    int ival = et.target.bra.left.ijk[stepidx];
    int kval = et.target.ket.left.ijk[stepidx]-1;

    std::stringstream ss;
    ss << std::string(20, ' ');
    ss << ETStepVar(et.target, true);

    // accumulate if is contracted integral workspace, otherwise just assign
    if(IsContArray(et.target.amlist()))
        ss << " += ";
    else
        ss << " = ";

    ss << "etfac[" << stepidx << "] * " << ETStepVar(et.src1, false);

    if(et.src2.bra.left && et.src2.ket.left)
        ss << " + " << ival << " * one_over_2q * " << ETStepVar(et.src2, false);
    if(et.src3.bra.left && et.src3.ket.left)
        ss << " + " << kval << " * one_over_2q * " << ETStepVar(et.src3, false);
    if(et.src4.bra.left && et.src4.ket.left)
        ss << " - p_over_q * " << ETStepVar(et.src4, false);
    ss << ";\n";

    return ss.str();
}






static void WriteVRRInfo(std::ostream & os, const std::pair<VRRMap, VRRReqMap> & vrrinfo, int L)
{
    if( vrrinfo.second.size() > 1 || (vrrinfo.second.size() == 1 && vrrinfo.second.begin()->first != 0) )
        os << "                    int idx = 0;\n\n";

    // accumulate S_0_0_0_0 if needed
    if(IsContArray({0, 0, 0, 0}))
    {
        os << "\n";
        os << "                    // Accumulating S_0_0_0_0 in contracted workspace\n";
        os <<"                    " << ArrVarName({0, 0, 0, 0}) << "[abcd] += " << AuxName(0) << "[0];\n";
        os << "\n";
    }

    // iterate over increasing am
    for(const auto & it3 : vrrinfo.second)
    {
        // i is the am
        int i = it3.first;

        // don't do zero - that is handled by the boys function stuff
        if(i == 0)
            continue;

        // greq is the actual requirements for this am
        const GaussianSet & greq = it3.second;

        // requirements for this am
        os << "                    // Forming " << AuxName(i) << "[" << (L-i+1) << " * " << greq.size() << "];\n";
        os << "                    // Needed from this AM:\n";
        for(const auto & it : greq)
            os << "                    //    " << it << "\n";
        os << "                    idx = 0;\n";
        os << "                    for(int m = 0; m < " << (L-i+1) << "; m++)  // loop over orders of boys function\n";
        os << "                    {\n";

        // iterate over the requirements
        // should be in order since it's a set
        for(const auto & it : greq)
        {
            // Get the stepping
            XYZStep step = vrrinfo.first.at(it);

            // and then step to the the required gaussians
            Gaussian g1 = it.StepDown(step, 1);
            Gaussian g2 = it.StepDown(step, 2);

            // get the number of gaussians in the previous two AM
            int ng1 = -1;
            int ng2 = -1;
            int g1_idx = -1;
            int g2_idx = -1;

            // Ugh indices. This is a little messy, but whatever
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
            //os << "                        NEED: " << g1 << "  ,  " << g2 << "\n";
            //os << "                         IDX: " << g1_idx << "  ,  " << g2_idx << "\n";
            //os << "                       OUTOF: " << ng1 << "  ,  " << ng2 << "\n";
            os << "                        " << AuxName(i);

            os << "[idx++] = P.PA_" << step << "[i] * ";

            if(g1)
                os << AuxName(i-1) << "[m * " << ng1 << " + " << g1_idx 
                   << "] - a_over_p * PQ_" << step << " * " << AuxName(i-1) << "[(m+1) * " << ng1 << " + " << g1_idx << "]";
            if(g2)
            {
                os << "\n";
                os << "                                     + " << vrr_i << " * one_over_2p * ( " << AuxName(i-2)
                   << "[m * " << ng2 << " +  " << g2_idx << "]"
                   << " - a_over_p * " << AuxName(i-2) << "[(m+1) * " << ng2 << " + " << g2_idx << "] )";
            }
            os << ";\n\n"; 
        }

        os << "                    }\n";


        // if this target is also a contracted array, accumulate there
        // note that contracted array has storage for all elements of a
        // given AM, not just what we calculated
        // (ie abcd*NCART(i) not abcd*greq.size() or something)
        if(IsContArray({i, 0, 0, 0}))
        {
            os << "\n";
            os << "                    // Accumulating in contracted workspace\n";
            os << "                    idx = 0;\n";
            int idx = 0;
            for(const auto & it : greq)
            {
                os <<"                    " << ArrVarName({i, 0, 0, 0}) << "[abcd * " << NCART(i) << " + " 
                   << it.idx() << "] += " << AuxName(i) << "[" << idx++ << "];\n";
            }
        }

        os << "\n\n";
    }

}





void Writer_Looped(std::ostream & os,
                   const QAMList & am,
                   const std::string & nameappend,
                   const BoysMap & bm,
                   VRR_Algorithm_Base & vrralgo,
                   ET_Algorithm_Base & etalgo,
                   HRR_Algorithm_Base & hrralgo)
{
    // add the am list to the contracted info, etc
    arrinfo.insert(am);
    continfo.insert(am);

    // for electron transfer
    // ETMap[am1][am2] is a list of steps for forming
    typedef std::map<int, std::map<int, ETStepList>> ETMap;


    int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);
    int ncart_bra = NCART(am[0]) * NCART(am[1]);

    const int L = am[0] + am[1] + am[2] + am[3];


    // Working backwards, I need:
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


    // we may need to add ( a s | as a top bra for AM quartets
    // that do not have HRR in the bra part
    // these might be ( X s |  or  | X s )
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

    // 3.) Top level gaussians from ET, sorted by AM
    ETReqMap etrm;
    for(const auto & it : etsl)
    {
        // should only add the first element of the bra
        if(it.src1.flags & QUARTET_ETTOPLEVEL)
            etrm[it.src1.bra.left.am()].insert(it.src1.bra.left);
        if(it.src2.flags & QUARTET_ETTOPLEVEL)
            etrm[it.src2.bra.left.am()].insert(it.src2.bra.left);
        if(it.src3.flags & QUARTET_ETTOPLEVEL)
            etrm[it.src3.bra.left.am()].insert(it.src3.bra.left);
        if(it.src4.flags & QUARTET_ETTOPLEVEL)
            etrm[it.src4.bra.left.am()].insert(it.src4.bra.left);
    }
    

    // 4.) VRR Steps
    // requirements for vrr are the elements of etrm
    ETReqMap vreq = etrm;

    // but if there are none, this of the type ( X 0 | 0 0 )
    if(vreq.size() == 0)
        vreq[L] = AllGaussiansForAM(L);

    std::pair<VRRMap, VRRReqMap> vrrinfo = vrralgo.CreateAllMaps(vreq);


    // add these to contracted array variable set
    for(const auto & it : hrrtopbras)
    for(const auto & it2 : hrrtopkets)
    {
        QAMList qam{it.first, 0, it2.first, 0};
        arrinfo.insert(qam);
        continfo.insert(qam);
    }

    for(const auto & it : hrrtopkets)
    {
        QAMList qam{am[0], am[1], it.first, 0};
        arrinfo.insert(qam);
        continfo.insert(qam);
    }


    // print out some info
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

    // some helper bools
    bool hasvrr = ( (vrrinfo.second.size() > 1) || (vrrinfo.second.size() == 1 && vrrinfo.second.begin()->first != 0) );
    bool haset = (etsl.size() > 0);
    bool hasoneover2p = true;  // TODO


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
    os << "    ASSUME_ALIGN(integrals)\n";
    os << "\n";
    os << "    const int nshell1234 = P.nshell12 * Q.nshell12;\n";
    os << "\n";
    os << "    memset(" << ArrVarName(am) << ", 0, nshell1234*" << ncart << "*sizeof(double));\n";
    os << "\n";
    os << "    // Holds AB_{xyz} and CD_{xyz} in a flattened fashion for later\n";
    os << "    double AB_x[nshell1234];  double CD_x[nshell1234];\n";
    os << "    double AB_y[nshell1234];  double CD_y[nshell1234];\n";
    os << "    double AB_z[nshell1234];  double CD_z[nshell1234];\n";
    os << "\n";
    os << "    int ab, cd, abcd;\n";
    os << "    int i, j;\n";
    os << "\n";
    os << "    // Workspace for contracted integrals\n";

    for(const auto & it : continfo)
        if(it != am)  // skip if this is the final integrals - don't re-declare
        {
            std::stringstream ss;
            ss << "nshell1234 * " << NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]);

            os << "    double " << ArrVarName(it) << "[" << ss.str() << "];\n";
            os << "    memset(" << ArrVarName(it) << ", 0, (" << ss.str() << ") * sizeof(double));\n";
            os << "\n";
        }

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
    os << "            AB_x[abcd] = P.AB_x[ab];  CD_x[abcd] = Q.AB_x[cd];\n";
    os << "            AB_y[abcd] = P.AB_y[ab];  CD_y[abcd] = Q.AB_y[cd];\n";
    os << "            AB_z[abcd] = P.AB_z[ab];  CD_z[abcd] = Q.AB_z[cd];\n";
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
        QAMList qam{greq.first, 0, 0, 0};
        os << "                    // AM = " << greq.first << ": Needed from this AM: " << greq.second.size() << "\n";
        os << "                    double " << AuxName(greq.first) << "[" << (L-greq.first+1) << " * " << greq.second.size() << "];\n";
        os << "\n";

        // add to contracted array list (but not contracted!)
        arrinfo.insert(qam);

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

    if(hasvrr)
    {
        os << "                    // for VRR\n";
        os << "                    const double one_over_p = 1.0 / P.alpha[i];\n";
        os << "                    const double a_over_p =  alpha * one_over_p;     // a/p from MEST\n";
        if(hasoneover2p)
            os << "                    const double one_over_2p = 0.5 * one_over_p;  // gets multiplied by i in VRR\n";
    }

    if(haset)
    {
        os << "                    // for electron transfer\n";
        os << "                    const double one_over_q = 1.0 / Q.alpha[j];\n";
        os << "                    const double one_over_2q = 0.5 * one_over_q;\n";
        os << "                    const double p_over_q = P.alpha[i] * one_over_q;\n";
        os << "\n";
        os << "                    const double etfac[3] = {\n";
        os << "                                             -(P.bAB_x[i] + Q.bAB_x[j]) * one_over_q,\n";
        os << "                                             -(P.bAB_y[i] + Q.bAB_y[j]) * one_over_q,\n";
        os << "                                             -(P.bAB_z[i] + Q.bAB_z[j]) * one_over_q,\n";
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

    for(int i = 0; i <= L; i++)
        os << bm.at(i)->code_line() << "\n";

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

    if(etsl.size() == 0)
        os << "                    //...nothing to do...\n";
    else
    {
        //Write_ElectronTransfer(os, etsl, etrm, etinit);
        // Start at the highest am, and work down
        for(const auto & it : etsl)
        {
            os << std::string(20, ' ') << "// " << it << "\n";
            os << ETStepString(it);
            os << "\n";
        }
    }
        
    os << "\n";
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
        os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << "    {\n";

        for(const auto & it : hrrtopkets)
        {
            os << "        // form " << ArrVarName({am[0], am[1], it.first, 0}) << "\n";
            os << "        for(int ni = 0; ni < " << it.second.size() << "; ++ni)\n";
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
