#include <cstdio>
#include <cmath>

#include "eri/eri.h"
#include "boys/boys.h"
#include "common.hpp"

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

#define MAXAM 2

bool ValidQuartet(std::array<int, 4> am)
{
    if(am[0] < am[1])
        return false;
    if(am[2] < am[3])
        return false;
    if( (am[0] + am[1]) < (am[2] + am[3]) )
        return false;
    if(am[0] < am[2])
        return false;
    return true;
}


int main(int argc, char ** argv)
{
    const int maxncart = pow(NCART(MAXAM), 4);

    // parse command line
    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
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
    std::array<int, 4> initam{0, 0, 0, 0}; // will be set later
    VecQuartet gshells(  ReadQuartets(initam, files, true) );


    /* Storage of integrals results */
    const int nshell1234 = gshells[0].size() * gshells[1].size() 
                         * gshells[2].size() * gshells[3].size();

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
 
        int ncart = NCART(i) * NCART(j) * NCART(k) * NCART(l);   


        // set all the am
        for(int m = 0; m < 4; m++)
            for(auto & it : gshells[m])
                it.am = am[m];


        // calculate
        ValeevIntegrals(gshells, res_valeev);

        printf("( %d %d | %d %d )\n", am[0], am[1], am[2], am[3]);
        for(int m = 0; m < ncart; m++)
            printf("    %20.8e\n", res_valeev[m]);
    }

    printf("\n");

    FreeQuartets(gshells);
    Valeev_Finalize();
    FREE(res_valeev);

    return 0;
}
