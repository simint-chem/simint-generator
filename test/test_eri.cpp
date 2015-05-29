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

    // basis functions file to read
    std::string basfile(argv[1]);

    // read in the shell info
    ShellMap shellmap = ReadBasis(basfile);
    ShellMap shellmap_erd = CopyShellMap(shellmap);

    // normalize
    for(auto & it : shellmap)
        normalize_gaussian_shells(it.second.size(), it.second.data());
    for(auto & it : shellmap_erd)
        normalize_gaussian_shells_erd(it.second.size(), it.second.data());
        

    // find the max dimensions
    std::array<int, 3> maxparams = FindMapMaxParams(shellmap);
    const int maxam = maxparams[0];
    const int maxnprim = maxparams[1];
    const int maxsize = maxparams[2];

    /* Storage of integrals */
    double * res_FO              = (double *)ALLOC(maxsize * sizeof(double));
    double * res_vref            = (double *)ALLOC(maxsize * sizeof(double));
    double * res_vref_flat       = (double *)ALLOC(maxsize * sizeof(double));
    double * res_liberd          = (double *)ALLOC(maxsize * sizeof(double));
    double * res_valeev          = (double *)ALLOC(maxsize * sizeof(double));

    // initialize stuff
    Valeev_Init();
    Boys_Init();
    ERD_Init(maxam, maxnprim, 1);
             
             
    printf("\n");
    printf("%17s    %10s  %10s  %10s  %10s    %10s  %10s  %10s  %10s\n", "Quartet", "MaxErr", "", "", "", 
                                                                "MaxRelErr", "", "", "");
    printf("%17s    %10s  %10s  %10s  %10s    %10s  %10s  %10s  %10s\n", "", "FO", "vref", "vrefF", "erd", 
                                                                "FO", "vref", "vrefF", "erd");

    // loop over all quartets, choosing only valid ones
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        std::array<int, 4> am{i, j, k, l};

        if(!ValidQuartet(am))
            continue;

        const int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);   

        const AlignedGaussianVec & it_i = shellmap[i];
        const AlignedGaussianVec & it_j = shellmap[j];
        const AlignedGaussianVec & it_k = shellmap[k];
        const AlignedGaussianVec & it_l = shellmap[l];

        // number of shells
        const int nshell1 = it_i.size();
        const int nshell2 = it_j.size();
        const int nshell3 = it_k.size();
        const int nshell4 = it_l.size();
        const int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;


        /////////////////////////////////
        // Calculate or read valeev
        // reference integrals
        /////////////////////////////////
        ValeevIntegrals(shellmap[i], shellmap[j],
                        shellmap[k], shellmap[l],
                        res_valeev, false);
        //ReadValeevIntegrals(basedir, am, res_valeev);


        /////////////////////////////////
        // Calculate with liberd
        /////////////////////////////////
        ERDIntegrals(shellmap_erd[i], shellmap_erd[j],
                     shellmap_erd[k], shellmap_erd[l],
                     res_liberd);


        /////////////////////////////////
        // calculate with my code
        /////////////////////////////////

        // set up shell pairs
        struct multishell_pair P = create_multishell_pair(it_i.size(), it_i.data(),
                                                          it_j.size(), it_j.data());
        struct multishell_pair Q = create_multishell_pair(it_k.size(), it_k.data(),
                                                          it_l.size(), it_l.data());

        struct multishell_pair_flat Pf = create_multishell_pair_flat(it_i.size(), it_i.data(),
                                                                     it_j.size(), it_j.data());
        struct multishell_pair_flat Qf = create_multishell_pair_flat(it_k.size(), it_k.data(),
                                                                     it_l.size(), it_l.data());

        // acutally calculate
        Integral_FO(P, Q, res_FO);
        Integral_vref(P, Q, res_vref);
        Integral_vref_flat(Pf, Qf, res_vref_flat);


        // chop
        const int arrlen = nshell1234 * ncart;
        Chop(res_valeev, arrlen);
        Chop(res_FO, arrlen);
        Chop(res_vref, arrlen);
        Chop(res_vref_flat, arrlen);


        // Analyze
        std::pair<double, double> err_FO        = CalcError(res_FO,        res_valeev,  arrlen);
        std::pair<double, double> err_vref      = CalcError(res_vref,      res_valeev,  arrlen);
        std::pair<double, double> err_vref_flat = CalcError(res_vref_flat, res_valeev,  arrlen);
        std::pair<double, double> err_erd       = CalcError(res_liberd,    res_valeev,  arrlen);

        printf("( %2d %2d | %2d %2d )    %10.3e  %10.3e  %10.3e  %10.3e    %10.3e  %10.3e  %10.3e  %10.3e\n", am[0], am[1], am[2], am[3],
                                                      err_FO.first, err_vref.first, err_vref_flat.first, err_erd.first,
                                                      err_FO.second, err_vref.second, err_vref_flat.second, err_erd.second);


        // For debugging
        //for(int i = 0; i < ncart * nshell1234; i++)
        //    printf("%25.15e  %25.15e  %25.15e\n", res_valeev[i], res_liberd[i], res_vref[i]);


        free_multishell_pair(P);
        free_multishell_pair(Q);
        free_multishell_pair_flat(Pf);
        free_multishell_pair_flat(Qf);
    }

    printf("\n");

    FreeShellMap(shellmap);
    FreeShellMap(shellmap_erd);

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
