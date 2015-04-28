#include <sstream>

#include "generator/Classes.hpp"
#include "generator/Boys.hpp"

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

static std::string ShellQuartetVarString(const ShellQuartet & q)
{
    std::stringstream ss;
    ss << "SQ_"  << q.amlist[0] << "_" << q.amlist[1] << "_"
                 << q.amlist[2] << "_" << q.amlist[3] << "_" << q.m;

    return ss.str();
}


/*
static std::string DoubletVarString(const Doublet & d, DoubletType type)
{
    std::stringstream ss;
    if(type == DoubletType::KET)
        ss << "KET";
    else
        ss << "BRA";

    ss << "_" << d.left.ijk[0] << d.left.ijk[1] << d.left.ijk[2]
       << "_" << d.right.ijk[0] << d.right.ijk[1] << d.right.ijk[2];


    return ss.str();
}

static std::string HRRBraetStepString(const HRRDoubletStep & hrr, int ketam)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "AB_";

    std::stringstream ss;
    ss << std::string(12, ' ');
    if(hrr.target.flags & DOUBLET_INITIAL || hrr.target.flags & DOUBLET_HRRTOPLEVEL)
    {
        ShellQuartet tshellq{{hrr.target.left.am(), hrr.target.right.am(), ketam, 0}, 0, 0};  // don't think I need to set flags
        ss << ShellQuartetVarString(tshellq) << "[somemath]"; 
    }
    else
        ss << "const double " << DoubletVarString(hrr.target, DoubletType::BRA);

    ss << " = ";

    if(hrr.src1.flags & DOUBLET_INITIAL || hrr.src1.flags & DOUBLET_HRRTOPLEVEL)
    {
        ShellQuartet tshellq{{hrr.src1.left.am(), hrr.src1.right.am(), ketam, 0}, 0, 0};  // don't think I need to set flags
        ss << ShellQuartetVarString(tshellq) << "[somemath]"; 
    }
    else
        ss << DoubletVarString(hrr.src1, DoubletType::BRA);

    ss << " + (" << xyztype << hrr.xyz << "[abcd]"
       << " * ";

    if(hrr.src2.flags & DOUBLET_INITIAL || hrr.src2.flags & DOUBLET_HRRTOPLEVEL)
    {
        ShellQuartet tshellq{{hrr.src2.left.am(), hrr.src2.right.am(), ketam, 0}, 0, 0};  // don't think I need to set flags
        ss << ShellQuartetVarString(tshellq) << "[somemath]"; 
    }
    else
        ss << DoubletVarString(hrr.src2, DoubletType::BRA);

    ss << ");";

    return ss.str();
}
*/


void Writer_Looped(std::ostream & os,
                   const QAMList & am,
                   const std::string & nameappend,
                   const BoysMap & bm,
                   const HRRBraKetStepList & hrrsteps)
{
    int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);

    // todo - calculate the max v needed for the boys function
    int maxv = 2;


    // I need:
    // 1.) HRR Top level Doublets
    DoubletSet hrrtopbras, hrrtopkets;
    ShellQuartetSet hrrtop;
    for(const auto & it : hrrsteps.first)
    {
        if(it.src1.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopbras.insert(it.src1);
        if(it.src2.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopbras.insert(it.src2);
    }

    for(const auto & it : hrrsteps.second)
    {
        if(it.src1.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopkets.insert(it.src1);
        if(it.src2.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopkets.insert(it.src2);
    }

    // these might be (ss| or |ss)
    if(hrrtopbras.size() == 0)
        hrrtopbras.insert({DoubletType::BRA, {0,0,0}, {0,0,0}, DOUBLET_HRRTOPLEVEL});
    if(hrrtopkets.size() == 0)
        hrrtopkets.insert({DoubletType::KET, {0,0,0}, {0,0,0}, DOUBLET_HRRTOPLEVEL});

    // 2.) HRR top level Shell Quartets
    ShellQuartetSet hrrtopshells;
    for(const auto & it : hrrtopbras)
        for(const auto & it2 : hrrtopkets)
            hrrtopshells.insert({it, it2});


    // 3.) Targets for the HRR bra-generation step
    //     These are the shell quartets corresponding to
    //     hrr top kets, combined with the final integral bra
    ShellQuartetSet hrrbratargets;
    for(const auto & it : hrrtopkets)
        hrrbratargets.insert({{am[0], am[1], it.left.am(), it.right.am()}, 0, 0});
    


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
    for(const auto & it : hrrtopshells)
        os << "    double " << ShellQuartetVarString(it) << "[" << NCART(it.amlist[0]+it.amlist[1])*NCART(it.amlist[2]+it.amlist[3]) << " * nshell1234];\n";

    os << "\n";
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
    os << "            AB_y[abcd] = P.AB_y[ab];  CD_x[abcd] = P.AB_y[cd];\n";
    os << "            AB_z[abcd] = P.AB_z[ab];  CD_x[abcd] = P.AB_z[cd];\n";
    os << "\n";
    os << "            for(i = abstart; i < abend; ++i)\n";
    os << "            {\n";
    os << "                for(j = cdstart; j < cdend; ++j)\n";
    os << "                {\n";
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
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Boys function section\n";
    os << "                    // Maximum v value: " << maxv << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // The paremeter to the boys function\n";
    os << "                    const double F_x = R2 * PQalpha_mul/PQalpha_sum;\n";
    os << "\n";

    Write_BoysFunction(os, bm, maxv);

    os << "\n\n";
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Vertical recurrance\n";
    os << "                    //////////////////////////////////////////////\n";
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
    os << "    // Bra targets\n";
    for(const auto & it : hrrbratargets)
        os << "    double " << ShellQuartetVarString(it) << "[" << NCART(it.amlist[0]+it.amlist[1])*NCART(it.amlist[2]+it.amlist[3]) << " * nshell1234];\n";

    os << "\n";    
    os << "    int startidx = 0;\n";
    os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
    os << "    {\n";

    for(const auto & it : hrrbratargets)
    {
        os << "        // form " << it << "\n";
        os << "        for(int ni = 0; ni < " << NCART(it.amlist[2]) << "; ++ni)\n";
        os << "        {\n";
        //for(const auto & hit : hrrinfo.bralist)
        //    os << HRRBraStepString(hit, it.amlist[2]) << "\n";
        os << "        }\n";
        os << "\n";
    }
    os << "\n";
    os << "        startidx += " << ncart << ";\n";
    os << "    }\n";
    os << "\n";
    os << "\n";
    os << "    return nshell1234;\n";
    os << "}\n";
    os << "\n";

    //os << "                    integrals[nint] += pfac * P.prefac[i] * Q.prefac[j] * Boys_F0_FO(x);\n";
}
