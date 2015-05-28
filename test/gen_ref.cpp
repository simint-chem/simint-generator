#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>

#include "eri/eri.h"
#include "boys/boys.h"
#include "test/common.hpp"
#include "test/valeev.hpp"

int main(int argc, char ** argv)
{
    // parse command line
    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
        return 1;
    }

    // files to read
    std::string basedir(argv[1]);

    // initialize valeev stuff
    Valeev_Init();

    // read in the shell info
    // angular momentum and normalization will be handled later
    std::array<int, 4> initam{0, 0, 0, 0};
    VecQuartet gshells_orig(  ReadQuartets(initam, basedir, false) );


    /* Storage of integrals results */
    const int maxncart = pow(NCART(MAXAM), 4);
    const int nshell1234 = gshells_orig[0].size() * gshells_orig[1].size() 
                         * gshells_orig[2].size() * gshells_orig[3].size();

    double * res_valeev = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));



    // loop over all quartets, choosing only valid ones
    for(int i = 0; i <= MAXAM; i++)
    for(int j = 0; j <= MAXAM; j++)
    for(int k = 0; k <= MAXAM; k++)
    for(int l = 0; l <= MAXAM; l++)
    {
        std::array<int, 4> am{i, j, k, l};

        if(!ValidQuartet(am))
            continue;
 
        printf("( %d %d | %d %d ) ... ", am[0], am[1], am[2], am[3]);
        const int ncart = NCART(i) * NCART(j) * NCART(k) * NCART(l);   
        const uint32_t nsize = ncart * nshell1234;  // make sure 32-bit uint

        // copy the shells, set all the am, and renormalize
        VecQuartet gshells(CopyQuartets(gshells_orig));

        for(int m = 0; m < 4; m++)
        {
            for(auto & it : gshells[m])
                it.am = am[m];
            normalize_gaussian_shells(gshells[m].size(), gshells[m].data());
        }


        // calculate
        ValeevIntegrals(gshells, res_valeev, true);


        // open the output file
        const char * amchar = "spdfghijklmnoqrtuvwxyzabe";
        std::stringstream ss;
        ss << basedir << "ref_" << amchar[am[0]] << "_"
                                << amchar[am[1]] << "_"
                                << amchar[am[2]] << "_"
                                << amchar[am[3]] << ".dat";

        std::ofstream out(ss.str().c_str(), std::ofstream::binary);
        out.write(reinterpret_cast<const char *>(&nsize), sizeof(uint32_t));
        out.write(reinterpret_cast<const char *>(res_valeev), nsize*sizeof(double));
        out.close();

        //for(int m = 0; m < ncart; m++)
        //    printf("%20.8e\n", res_valeev[m]);

        FreeQuartets(gshells);
        printf(" done\n");
    }

    FreeQuartets(gshells_orig);

    Valeev_Finalize();

    FREE(res_valeev);

    return 0;
}
