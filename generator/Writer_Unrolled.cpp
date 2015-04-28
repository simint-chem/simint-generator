#include <sstream>

#include "generator/Classes.hpp"
#include "generator/Boys.hpp"

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

static const char * amlist = "spdfghijklmnoqrtuvwxyzabe";


// This is split out to allow for future optimization
static void Write_BoysFunction(std::ostream & os,
                               const BoysMap & bm,
                               int maxv)
{
    for(int i = 0; i <= maxv; i++)
        os << bm.at(i)->code_line() << "\n";
}


static std::string QuartetVarString(const Quartet & q)
{
    std::stringstream ss;

    // is this one of the final integrals
    if(q.flags & QUARTET_INITIAL)
        ss << "integrals[startidx + " << q.idx() << "]";
    else
    {
        // if not, it's a temp variable
        ss << "Q_" 
           << q.bra.left.ijk[0] << q.bra.left.ijk[1] << q.bra.left.ijk[2]      << "_"
           << q.bra.right.ijk[0] << q.bra.right.ijk[1] << q.bra.right.ijk[2]   << "_"
           << q.ket.left.ijk[0] << q.ket.left.ijk[1] << q.ket.left.ijk[2]      << "_"
           << q.ket.right.ijk[0] << q.ket.right.ijk[1] << q.ket.right.ijk[2]   << "_"
           << q.m;
    }

    return ss.str();
}


static std::string HRRQuartetStepString(const HRRQuartetStep & hrr)
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = (hrr.steptype == DoubletType::BRA ? "AB_" : "CD_");


    std::stringstream ss;

    // target is a temp variable
    if(!(hrr.target.flags & QUARTET_INITIAL))
        ss << "const double ";  

    ss << QuartetVarString(hrr.target) 
       << " = " << QuartetVarString(hrr.src1)
       << " + (" << xyztype << hrr.xyz << "[abcd]" 
       << " * " << QuartetVarString(hrr.src2) << ");";

    // is it one of the final integrals?
    // if so, add a comment about which one it is
    if(hrr.target.flags & QUARTET_INITIAL)
        ss << "    // " << hrr.target.str();

    return ss.str();
}


void Writer_Unrolled(std::ostream & os,
                     const QAMList & am,
                     const std::string & nameappend,
                     const BoysMap & bm,
                     const HRRQuartetStepList & hrrsteps)
{
    int nam1 = NCART(am[0]);
    int nam2 = NCART(am[1]);
    int nam3 = NCART(am[2]);
    int nam4 = NCART(am[3]);
    int ncart = nam1 * nam2 * nam3 * nam4;

    // todo - calculate the max v needed for the boys function
    int maxv = 2;

    // calculate the top requirements for the HRR steps
    QuartetSet hrrtopreq;
    for(const auto & it : hrrsteps)
    {
        if(it.src1.flags & QUARTET_HRRTOPLEVEL)
            hrrtopreq.insert(it.src1);
        if(it.src2.flags & QUARTET_HRRTOPLEVEL)
            hrrtopreq.insert(it.src2);
    }


    // set up the function to make it look pretty
    std::stringstream ss;
    ss << "int eri_" << nameappend << "_"
       << amlist[am[0]] << amlist[am[1]]
       << amlist[am[2]] << amlist[am[3]] << "(";

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

    for(const auto & it : hrrtopreq)
    {
        os << "    double " << QuartetVarString(it) << "[nshell1234];\n";
        os << "    memset(" << QuartetVarString(it) << ", 0, nshell1234 * sizeof(double));\n";
        os << "\n";
    }

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
    os << "    // Steps: " << hrrsteps.size() << "\n";
    os << "    //////////////////////////////////////////////\n";
    os << "\n";
    os << "    int startidx = 0;\n";
    os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
    os << "    {\n";
    for(const auto & h : hrrsteps)
        os << "        " << HRRQuartetStepString(h) << "\n";    
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
