#include <cstdio>
#include <cmath>
#include <sstream>

#include "eri/eri.h"
#include "boys/boys.h"
#include "common.hpp"

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

int main(int argc, char ** argv)
{
    // parse command line
    if(argc != 6)
    {
        printf("Give me 5 arguments! I got %d\n", argc-1);
        return 1;
    }

    // files to read
    std::string basedir(argv[1]);
    if(basedir[basedir.size()-1] != '/')
        basedir = basedir + '/';

    std::array<std::string, 4> files{basedir + "1.dat",
                                     basedir + "2.dat",
                                     basedir + "3.dat",
                                     basedir + "4.dat"};

    // initialize valeev stuff
    Valeev_Init();

    // read in the shell info
    std::array<int, 4> am{atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5])};
    VecQuartet gshells(  ReadQuartets(am, files, true) );


    /* Storage of integrals results */
    const int nshell1234 = gshells[0].size() * gshells[1].size() 
                         * gshells[2].size() * gshells[3].size();

    const int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);
    const int nsize = ncart * nshell1234;
    double * res_valeev = (double *)ALLOC(nsize * sizeof(double));


    // calculate
    ValeevIntegrals(gshells, res_valeev);


    // open the output file
    const char * amchar = "spdfghijklmnoqrtuvwxyzabe";
    std::stringstream ss;
    ss << basedir << "ref_" << amchar[am[0]] << "_"
                            << amchar[am[1]] << "_"
                            << amchar[am[2]] << "_"
                            << amchar[am[3]] << ".dat";

    uint32_t nsize32 = nsize; // just to make sure
    std::ofstream out(ss.str().c_str(), std::ofstream::binary);
    out.write(reinterpret_cast<const char *>(&nsize32), sizeof(uint32_t));
    out.write(reinterpret_cast<const char *>(res_valeev), nsize*sizeof(double));
    out.close();

    //printf("( %d %d | %d %d )\n", am[0], am[1], am[2], am[3]);
    //for(int m = 0; m < ncart; m++)
    //    printf("%20.8e\n", res_valeev[m]);

    FreeQuartets(gshells);
    Valeev_Finalize();
    FREE(res_valeev);

    return 0;
}
