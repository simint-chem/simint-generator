#include <cstdio>
#include <cmath>

#include "boys/boys.h"
#include "test/common.hpp"
#include "test/valeev.hpp"
#include "test/ERD.hpp"
#include "test/Libint2.hpp"
#include "vectorization/vectorization.h"

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
    double * res              = (double *)ALLOC(maxsize * sizeof(double));
    double * res_liberd          = (double *)ALLOC(maxsize * sizeof(double));
    double * res_libint          = (double *)ALLOC(maxsize * sizeof(double));
    double * res_valeev          = (double *)ALLOC(maxsize * sizeof(double));

    // initialize stuff
    Valeev_Init();
    Boys_Init();
    LIBINT2_PREFIXED_NAME(libint2_static_init)();

    ERD_ERI erd(maxam, maxnprim, 1);
    Libint2_ERI libint(maxam, maxnprim);
             
             
    printf("\n");
    printf("%17s    %10s  %10s  %10s   %10s  %10s  %10s\n", 
                                                         "Quartet", "MaxErr", "", "", 
                                                         "MaxRelErr", "", "");
    printf("%17s    %10s  %10s  %10s   %10s  %10s  %10s\n", "",
                                                         "Me", "erd", "libint", 
                                                         "Me", "erd", "libint");


    // Read the reference integrals
    //RefIntegralReader refint(basfile);

    // loop over all quartets, choosing only valid ones
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        if(!ValidQuartet(i, j, k, l))
            continue;

        const int ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);

        const GaussianVec & it_i = shellmap[i];
        const GaussianVec & it_j = shellmap[j];
        const GaussianVec & it_k = shellmap[k];
        const GaussianVec & it_l = shellmap[l];

        const GaussianVec & erd_it_i = shellmap_erd[i];
        const GaussianVec & erd_it_j = shellmap_erd[j];
        const GaussianVec & erd_it_k = shellmap_erd[k];
        const GaussianVec & erd_it_l = shellmap_erd[l];

        // number of shells
        const int nshell1 = it_i.size();
        const int nshell2 = it_j.size();
        const int nshell3 = it_k.size();
        const int nshell4 = it_l.size();
        const int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;

        const int arrlen = nshell1234 * ncart1234;

        /////////////////////////////////
        // Calculate or read valeev
        // reference integrals
        /////////////////////////////////
        ValeevIntegrals(it_i, it_j, it_k, it_l, res_valeev, false);
        //refint.ReadNext(res_valeev, arrlen);

        //////////////////////
        // set up shell pairs
        //////////////////////
        struct multishell_pair P = create_multishell_pair(it_i.size(), it_i.data(),
                                                          it_j.size(), it_j.data());
        struct multishell_pair Q = create_multishell_pair(it_k.size(), it_k.data(),
                                                          it_l.size(), it_l.data());


        /////////////////////////////////
        // Calculate with liberd and libint2
        /////////////////////////////////
        erd.Integrals(erd_it_i, erd_it_j, erd_it_k, erd_it_l, res_liberd);
        libint.Integrals(P, Q, res_libint);

        /////////////////////////////////
        // calculate with my code
        /////////////////////////////////

        /*
        struct multishell_pair_flat Pf = create_multishell_pair_flat(it_i.size(), it_i.data(),
                                                                     it_j.size(), it_j.data());
        struct multishell_pair_flat Qf = create_multishell_pair_flat(it_k.size(), it_k.data(),
                                                                     it_l.size(), it_l.data());
        */

        // acutally calculate
        Integral(P, Q, res);


        // chop
        //Chop(res_valeev, arrlen);
        //Chop(res_libint, arrlen);
        //Chop(res, arrlen);


        // Analyze
        std::pair<double, double> err           = CalcError(res,         res_valeev,  arrlen);
        std::pair<double, double> err_erd       = CalcError(res_liberd,  res_valeev,  arrlen);
        std::pair<double, double> err_libint    = CalcError(res_libint,  res_valeev,  arrlen);

        printf("( %2d %2d | %2d %2d )    %10.3e  %10.3e  %10.3e    %10.3e  %10.3e  %10.3e\n", i, j, k, l,
                                                      err.first,  err_erd.first,  err_libint.first,
                                                      err.second, err_erd.second, err_libint.second);


        // For debugging
        // For debugging
        int m = 0;
        for(int m1 = 0; m1 < nshell1; m1++)
        for(int m2 = 0; m2 < nshell2; m2++)
        for(int m3 = 0; m3 < nshell3; m3++)
        for(int m4 = 0; m4 < nshell4; m4++)
        for(int c1 = 0; c1 < NCART(i); c1++)
        for(int c2 = 0; c2 < NCART(j); c2++)
        for(int c3 = 0; c3 < NCART(k); c3++)
        for(int c4 = 0; c4 < NCART(l); c4++, ++m)
        {
            double diff = fabs(res_valeev[m] - res[m]);
            double rdiff = fabs(diff / res_valeev[m]);
            double ldiff = fabs(res_valeev[m] - res_libint[m]);
            double lrdiff = fabs(ldiff / res_valeev[m]);

            if( (diff > 1e-14 && rdiff > 1e-8) || (ldiff > 1e-14 && lrdiff > 1e-8) )
                printf("%d %d %d %d : %d %d %d %d  %25.16e  %25.16e  %25.16e   %25.16e  %25.16e  %25.16e  %25.16e\n",
                                                                                          m1, m2, m3, m4, c1, c2, c3, c4,
                                                                                          res_valeev[m], res_libint[m], res[m], 
                                                                                          diff, rdiff, ldiff, lrdiff);
        }


        free_multishell_pair(P);
        free_multishell_pair(Q);
        //free_multishell_pair_flat(Pf);
        //free_multishell_pair_flat(Qf);
    }

    printf("\n");

    FreeShellMap(shellmap);
    FreeShellMap(shellmap_erd);

    Valeev_Finalize();
    Boys_Finalize();

    FREE(res_valeev);
    FREE(res_liberd);
    FREE(res_libint);
    FREE(res);

    return 0;
}
