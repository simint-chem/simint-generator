#include <ostream>
#include <array>
#include <string>
#include <sstream>

#include "generator/Classes.hpp"
#include "generator/Boys.hpp"

const char * amlist = "spdfghijklmnoqrtuvwxyzabe";


// This is split out to allow for future optimization
void Write_BoysFunction(std::ostream & os,
                        const BoysMap & bm,
                        int maxv)
{
    for(int i = 0; i <= maxv; i++)
        os << bm.at(i)->code_line() << "\n";
}




void Write_Generic(std::ostream & os,
                   const std::array<int, 4> & am,
                   const std::string & nameappend,
                   const BoysMap & bm,
                   const VRRInfo & vrrinfo,
                   const ETInfo & etinfo,
                   const HRRInfo & hrrinfo)
{

    // set up the function to make it look pretty
    std::stringstream ss;
    ss << "int eri_" << nameappend << "_"
       << amlist[am[0]] << amlist[am[1]]
       << amlist[am[2]] << amlist[am[3]] << "(";

    std::string funcline = ss.str();
    std::string indent(funcline.length(), ' ');

    // start output to the file
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
    os << "    int ab, cd;\n";
    os << "    int i, j;\n";
    os << "    int nint = 0;\n";
    os << "\n";
    os << "    // Top level HRR requirements. Contracted integrals are accumulated here\n";

    for(auto & it : hrrinfo.topreq)
    {
        os << "    double " << it.code_var() << "[nshell1234];\n";
        os << "    memset(" << it.code_var() << ", 0, nshell1234 * sizeof(double));\n";
    }

    os << "\n";
    os << "    for(ab = 0; ab < P.nshell12; ++ab)\n";
    os << "    {\n";
    os << "        const int abstart = P.primstart[ab];\n";
    os << "        const int abend = P.primend[ab];\n";
    os << "\n";
    os << "        // this should have been set/aligned in fill_multishell_pair or something else\n";
    os << "        ASSUME(abstart%SIMD_ALIGN_DBL == 0);\n";
    os << "\n";
    os << "        for(cd = 0; cd < Q.nshell12; ++cd)\n";
    os << "        {\n";
    os << "            const int cdstart = Q.primstart[cd];\n";
    os << "            const int cdend = Q.primend[cd];\n";
    os << "\n";
    os << "            // this should have been set/aligned in fill_multishell_pair or something else\n";
    os << "            ASSUME(cdstart%SIMD_ALIGN_DBL == 0);\n";
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
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // The paremeter to the boys function\n";
    os << "                    const double F_x = R2 * PQalpha_mul/PQalpha_sum;\n";
    os << "\n";

    Write_BoysFunction(os, bm, vrrinfo.maxv);

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
    os << "\n";
    os << "            ++nint;\n";
    os << "\n";
    os << "        }\n";
    os << "    }\n";
    os << "\n";
    os << "    return nshell1234;\n";
    os << "}\n";
    os << "\n";

    //os << "                    integrals[nint] += pfac * P.prefac[i] * Q.prefac[j] * Boys_F0_FO(x);\n";
}
