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
    ss << "SQ_"  << q.bra.amlist[0] << "_" << q.bra.amlist[1] << "_"
                 << q.ket.amlist[0] << "_" << q.ket.amlist[1] << "_" << q.m;

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
*/


static std::string HRRBraStepVariable(const Doublet & d, const ShellDoublet & ket, bool istarget)
{
    std::stringstream ss;

    if(d.flags & DOUBLET_HRRTOPLEVEL && ket.flags & DOUBLET_HRRTOPLEVEL)
    {
        ShellQuartet sq(d, ket, 0);
        ss << ShellQuartetVarString(sq) << "[abcd * " << d.ncart() * ket.ncart() 
                                        << " + " << d.idx() << " * " << ket.ncart() << " + ni]";
    }

    /*
    else if(d.flags & DOUBLET_INITIAL & ket.flags & DOUBLET_INITIAL)
    {
        ShellQuartet sq(d, ket, 0);
        ss << ShellQuartetVarString(sq) << "[INITmath]";
    }
    */

    else if(d.flags & DOUBLET_INITIAL)
    {
        // this is a 'bra target'
        // determine the shell quartet
        ShellQuartet sq(d, ket, 0);
        ss << ShellQuartetVarString(sq) << "[abcd * " << d.ncart() * ket.ncart() 
                                        << " + " << d.idx() << " * " << ket.ncart() << " + ni]";
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

static std::string HRRKetStepVariable(const Doublet & d, const ShellDoublet & bra, bool istarget)
{
    std::stringstream ss;


    // bra should always be initial

    if(d.flags & DOUBLET_HRRTOPLEVEL && bra.flags & DOUBLET_HRRTOPLEVEL)
    {
        ShellQuartet sq(bra, d, 0);
        ss << ShellQuartetVarString(sq) << "[abcd * " << d.ncart() * bra.ncart() 
                                        << " + abi * " << bra.ncart() << " + " << d.idx() << "]";
    }

    else if(d.flags & DOUBLET_HRRTOPLEVEL)  // initial bra, hrr top level ket = a bra target
    {
        ShellQuartet sq(bra, d, 0);
        ss << ShellQuartetVarString(sq) << "[abcd * " << d.ncart() * bra.ncart() 
                                        << " + abi * " << bra.ncart() << " + " << d.idx() << "]";
    }

    else if(d.flags & DOUBLET_INITIAL)
    {
        ss << "integrals" << "[abcd * " << d.ncart() * bra.ncart() 
                                        << " + ni * " << d.ncart() << " + " << d.idx() << "]";
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


static std::string HRRBraStepString(const HRRDoubletStep & hrr, const ShellDoublet & ket)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "AB_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRBraStepVariable(hrr.target, ket, true);

    ss << " = ";
    ss << HRRBraStepVariable(hrr.src1, ket, false);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRBraStepVariable(hrr.src2, ket, false);
    ss << " );";

    return ss.str();
}

static std::string HRRKetStepString(const HRRDoubletStep & hrr, const ShellDoublet & bra)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "AB_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRKetStepVariable(hrr.target, bra, true);

    ss << " = ";
    ss << HRRKetStepVariable(hrr.src1, bra, false);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRKetStepVariable(hrr.src2, bra, false);
    ss << " );";

    return ss.str();
}


void Writer_Looped(std::ostream & os,
                   const QAMList & am,
                   const std::string & nameappend,
                   const BoysMap & bm,
                   const HRRBraKetStepList & hrrsteps)
{
    //int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);

    // todo - calculate the max v needed for the boys function
    int maxv = 2;


    // I need:
    // 1.) HRR Top level Doublets
    DoubletSet hrrtopbras, hrrtopkets;
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
        {
            int flag = 0;
            if(it.flags & DOUBLET_INITIAL && it2.flags & DOUBLET_INITIAL)
                flag |= QUARTET_INITIAL;
            if(it.flags & DOUBLET_HRRTOPLEVEL && it2.flags & DOUBLET_HRRTOPLEVEL)
                flag |= QUARTET_HRRTOPLEVEL;
            hrrtopshells.insert({it, it2, 0});
        }


    // 3.) Targets for the HRR bra-generation step
    //     These are the shell quartets corresponding to
    //     hrr top kets, combined with the final integral bra
    ShellQuartetSet hrrbratargets;
    for(const auto & it : hrrtopkets)
        hrrbratargets.insert(
                             {
                                {DoubletType::BRA, {am[0], am[1]}, DOUBLET_INITIAL}, // initial bra
                                it, // this ket
                                0   // m = 0 at this point
                             }
                            );
    


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
        os << "    double " << ShellQuartetVarString(it) << "[" << it.ncart() << " * nshell1234];\n";

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
    os << "            AB_y[abcd] = P.AB_y[ab];  CD_y[abcd] = P.AB_y[cd];\n";
    os << "            AB_z[abcd] = P.AB_z[ab];  CD_z[abcd] = P.AB_z[cd];\n";
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
        os << "    double " << ShellQuartetVarString(it) << "[" << it.ncart() << " * nshell1234];\n";

    os << "\n";    
    os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
    os << "    {\n";

    for(const auto & it : hrrbratargets)
    {
        os << "        // form " << it << "\n";
        os << "        for(int ni = 0; ni < " << NCART(it.ket.amlist[0]) << "; ++ni)\n";
        os << "        {\n";
        for(const auto & hit : hrrsteps.first)
        {
            os << std::string(12, ' ') << "// " << hit << "\n";
            os << HRRBraStepString(hit, it.ket) << "\n\n";
        }
        os << "        }\n";
        os << "\n";
    }
    os << "\n";
    os << "    }\n";
    os << "\n";
    os << "\n";
    os << "    //////////////////////////////////////////////\n";
    os << "    // Contracted integrals: Horizontal recurrance\n";
    os << "    // Ket part\n";
    os << "    // Steps: " << hrrsteps.second.size() << "\n";
    os << "    // Forming final integrals\n";
    os << "    //////////////////////////////////////////////\n";
    os << "\n";
    os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
    os << "    {\n";
    os << "        for(int abi = 0; abi < " << NCART(am[0]) * NCART(am[1]) << "; ++abi)\n"; 
    os << "        {\n"; 

    ShellDoublet initbra{DoubletType::BRA, {am[0], am[1]}, DOUBLET_INITIAL};
    for(const auto & hit : hrrsteps.second)
    {
        os << std::string(12, ' ') << "// " << hit << "\n";
        os << HRRKetStepString(hit, initbra) << "\n\n";
    }
    os << "        }\n"; 
    os << "    }\n";
    os << "\n";
    os << "\n";
    os << "    return nshell1234;\n";
    os << "}\n";
    os << "\n";

    //os << "                    integrals[nint] += pfac * P.prefac[i] * Q.prefac[j] * Boys_F0_FO(x);\n";
}
