#include <cstdio>
#include <cmath>

#include "boys/boys.h"
#include "test/common.hpp"
#include "test/valeev.hpp"
#include "test/erd_interface.hpp"


int main(int argc, char ** argv)
{
    // set up the function pointers
    Init_Test();

    // parse command line
    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
        return 1;
    }

    // files to read
    std::string basedir(argv[1]);


    // read in the shell info
    // angular momentum and normalization will be handled later
    std::array<int, 4> initam{0, 0, 0, 0};
    VecQuartet gshells_orig(  ReadQuartets(initam, basedir, false) );

    std::array<int, 4> nshell{gshells_orig[0].size(), gshells_orig[1].size(), 
                              gshells_orig[2].size(), gshells_orig[3].size()};

    // initialize stuff
    // in this case, all shells have the same number of primitives
    // so we can initialize ERD here
    Valeev_Init();
    Boys_Init();
    ERD_Init(MAXAM, gshells_orig[0][0].nprim, 1, 
             MAXAM, gshells_orig[1][0].nprim, 1, 
             MAXAM, gshells_orig[2][0].nprim, 1, 
             MAXAM, gshells_orig[3][0].nprim, 1);


    /* Storage of integrals results */
    const int maxncart = pow(NCART(MAXAM), 4);
    const int nshell1234 = nshell[0] * nshell[1] * nshell[2] * nshell[3];
    double * res_FO              = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    double * res_vref            = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    double * res_vref_flat       = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    double * res_liberd          = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    double * res_valeev          = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));


    printf("\n");
    printf("%17s    %8s  %8s  %8s  %8s    %8s %8s  %8s  %8s\n", "Quartet", "FO MaxErr", "vref MaxErr", "vrefF MaxErr", "erd MaxErr", 
                                                                "FO MaxRelErr", "vref MaxRelErr", "vrefF MaxRelErr", "erd MaxRelErr");

    // loop over all quartets, choosing only valid ones
    for(int i = 0; i <= MAXAM; i++)
    for(int j = 0; j <= MAXAM; j++)
    for(int k = 0; k <= MAXAM; k++)
    for(int l = 0; l <= MAXAM; l++)
    {
        std::array<int, 4> am{i, j, k, l};

        if(!ValidQuartet(am))
            continue;
 
        const int ncart = NCART(i) * NCART(j) * NCART(k) * NCART(l);   

        // copy the shells, set all the am, and normalize
        VecQuartet gshells(CopyQuartets(gshells_orig));
        VecQuartet gshells_erd(CopyQuartets(gshells_orig));

        for(int m = 0; m < 4; m++)
        {
            for(auto & it : gshells[m])
                it.am = am[m];
            for(auto & it : gshells_erd[m])
                it.am = am[m];

            //normalize_gaussian_shells(gshells[m].size(), gshells[m].data());
            normalize_gaussian_shells_erd(gshells_erd[m].size(), gshells_erd[m].data());
        }

        
        // set up shell pairs
        struct multishell_pair P = create_multishell_pair(nshell[0], gshells[0].data(),
                                                          nshell[1], gshells[1].data());
        struct multishell_pair Q = create_multishell_pair(nshell[2], gshells[2].data(),
                                                          nshell[3], gshells[3].data());

        struct multishell_pair_flat Pf = create_multishell_pair_flat(nshell[0], gshells[0].data(),
                                                                     nshell[1], gshells[1].data());
        struct multishell_pair_flat Qf = create_multishell_pair_flat(nshell[2], gshells[2].data(),
                                                                     nshell[3], gshells[3].data());


        // calculate with my code
        //Integral_FO(P, Q, res_FO);
        Integral_vref(P, Q, res_vref);
        Integral_vref_flat(Pf, Qf, res_vref_flat);

        // test with liberd
        ERDIntegrals(gshells_erd, res_liberd);

        // Calculate or read valeev reference integrals
        ValeevIntegrals(gshells, res_valeev, true);
        //ReadValeevIntegrals(basedir, am, res_valeev);

        std::pair<double, double> err_FO = CalcError(res_FO, res_valeev, nshell1234 * ncart);
        std::pair<double, double> err_vref = CalcError(res_vref, res_valeev, nshell1234 * ncart);
        std::pair<double, double> err_vref_flat = CalcError(res_vref_flat, res_valeev, nshell1234 * ncart);
        std::pair<double, double> err_erd = CalcError(res_liberd, res_valeev, nshell1234 * ncart);

        printf("( %2d %2d | %2d %2d )    %8.3e  %8.3e  %8.3e  %8.3e    %8.3e  %8.3e  %8.3e  %8.3e\n", am[0], am[1], am[2], am[3],
                                                      err_FO.first, err_vref.first, err_vref_flat.first, err_erd.first,
                                                      err_FO.second, err_vref.second, err_vref_flat.second, err_erd.second);


        for(int i = 0; i < ncart * nshell1234; i++)
            printf("%22.8e  %22.8e       %22.8e %22.8e\n", res_liberd[i], res_valeev[i], res_liberd[i]/res_valeev[i], res_valeev[i]/res_liberd[i]);
        free_multishell_pair(P);
        free_multishell_pair(Q);
        free_multishell_pair_flat(Pf);
        free_multishell_pair_flat(Qf);
        FreeQuartets(gshells);
        FreeQuartets(gshells_erd);

    }

    printf("\n");


    FreeQuartets(gshells_orig);

    Valeev_Finalize();
    Boys_Finalize();
    ERD_Finalize();

    FREE(res_valeev);
    FREE(res_liberd);
    FREE(res_vref);
    FREE(res_vref_flat);
    FREE(res_FO);

    return 0;
}
