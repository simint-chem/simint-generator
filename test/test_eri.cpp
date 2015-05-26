#include <cstdio>
#include <cmath>

#include "eri/eri.h"
#include "boys/boys.h"
#include "common.hpp"


typedef int (*erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);
typedef int (*eriflatfunc)(struct multishell_pair_flat const, struct multishell_pair_flat const, double * const restrict);


int eri_notyetimplemented(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict dummy)
{
    printf("****************************\n");
    printf("*** NOT YET IMPLEMENTED! ***\n");
    printf("***  ( %2d %2d | %2d %2d )   ***\n", P.am1, P.am2, Q.am1, Q.am2);
    printf("****************************\n");
    exit(1);
    return 0;
}

int eriflat_notyetimplemented(struct multishell_pair_flat const P,
                              struct multishell_pair_flat const Q,
                              double * const restrict dummy)
{
    printf("****************************\n");
    printf("*** NOT YET IMPLEMENTED! ***\n");
    printf("***  ( %2d %2d | %2d %2d )   ***\n", P.am1, P.am2, Q.am1, Q.am2);
    printf("****************************\n");
    exit(1);
    return 0;
}



int main(int argc, char ** argv)
{
    // set up the function pointers
    erifunc funcs_FO[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];
    erifunc funcs_vref[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];
    eriflatfunc funcs_vref_flat[MAXAM+1][MAXAM+1][MAXAM+1][MAXAM+1];

    for(int i = 0; i <= MAXAM; i++)
    for(int j = 0; j <= MAXAM; j++)
    for(int k = 0; k <= MAXAM; k++)
    for(int l = 0; l <= MAXAM; l++)
    {
        funcs_FO[i][j][k][l] = eri_notyetimplemented; 
        funcs_vref[i][j][k][l] = eri_notyetimplemented; 
        funcs_vref_flat[i][j][k][l] = eriflat_notyetimplemented; 
    }


    funcs_FO[0][0][0][0] = eri_FO_s_s_s_s;
    funcs_FO[1][0][0][0] = eri_FO_p_s_s_s;
    funcs_FO[1][0][1][0] = eri_FO_p_s_p_s;
    funcs_FO[1][1][0][0] = eri_FO_p_p_s_s;
    funcs_FO[1][1][1][0] = eri_FO_p_p_p_s;
    funcs_FO[1][1][1][1] = eri_FO_p_p_p_p;
    funcs_FO[2][0][0][0] = eri_FO_d_s_s_s;
    funcs_FO[2][0][1][0] = eri_FO_d_s_p_s;
    funcs_FO[2][0][1][1] = eri_FO_d_s_p_p;
    funcs_FO[2][0][2][0] = eri_FO_d_s_d_s;
    funcs_FO[2][1][0][0] = eri_FO_d_p_s_s;
    funcs_FO[2][1][1][0] = eri_FO_d_p_p_s;
    funcs_FO[2][1][1][1] = eri_FO_d_p_p_p;
    funcs_FO[2][1][2][0] = eri_FO_d_p_d_s;
    funcs_FO[2][1][2][1] = eri_FO_d_p_d_p;
    funcs_FO[2][2][0][0] = eri_FO_d_d_s_s;
    funcs_FO[2][2][1][0] = eri_FO_d_d_p_s;
    funcs_FO[2][2][1][1] = eri_FO_d_d_p_p;
    funcs_FO[2][2][2][0] = eri_FO_d_d_d_s;
    funcs_FO[2][2][2][1] = eri_FO_d_d_d_p;
    funcs_FO[2][2][2][2] = eri_FO_d_d_d_d;

    funcs_vref[0][0][0][0] = eri_vref_s_s_s_s;
    funcs_vref[1][0][0][0] = eri_vref_p_s_s_s;
    funcs_vref[1][0][1][0] = eri_vref_p_s_p_s;
    funcs_vref[1][1][0][0] = eri_vref_p_p_s_s;
    funcs_vref[1][1][1][0] = eri_vref_p_p_p_s;
    funcs_vref[1][1][1][1] = eri_vref_p_p_p_p;
    funcs_vref[2][0][0][0] = eri_vref_d_s_s_s;
    funcs_vref[2][0][1][0] = eri_vref_d_s_p_s;
    funcs_vref[2][0][1][1] = eri_vref_d_s_p_p;
    funcs_vref[2][0][2][0] = eri_vref_d_s_d_s;
    funcs_vref[2][1][0][0] = eri_vref_d_p_s_s;
    funcs_vref[2][1][1][0] = eri_vref_d_p_p_s;
    funcs_vref[2][1][1][1] = eri_vref_d_p_p_p;
    funcs_vref[2][1][2][0] = eri_vref_d_p_d_s;
    funcs_vref[2][1][2][1] = eri_vref_d_p_d_p;
    funcs_vref[2][2][0][0] = eri_vref_d_d_s_s;
    funcs_vref[2][2][1][0] = eri_vref_d_d_p_s;
    funcs_vref[2][2][1][1] = eri_vref_d_d_p_p;
    funcs_vref[2][2][2][0] = eri_vref_d_d_d_s;
    funcs_vref[2][2][2][1] = eri_vref_d_d_d_p;
    funcs_vref[2][2][2][2] = eri_vref_d_d_d_d;

    funcs_vref_flat[0][0][0][0] = eri_vref_flat_s_s_s_s;
    funcs_vref_flat[1][0][0][0] = eri_vref_flat_p_s_s_s;
    funcs_vref_flat[1][0][1][0] = eri_vref_flat_p_s_p_s;
    funcs_vref_flat[1][1][0][0] = eri_vref_flat_p_p_s_s;
    funcs_vref_flat[1][1][1][0] = eri_vref_flat_p_p_p_s;
    funcs_vref_flat[1][1][1][1] = eri_vref_flat_p_p_p_p;
    funcs_vref_flat[2][0][0][0] = eri_vref_flat_d_s_s_s;
    funcs_vref_flat[2][0][1][0] = eri_vref_flat_d_s_p_s;
    funcs_vref_flat[2][0][1][1] = eri_vref_flat_d_s_p_p;
    funcs_vref_flat[2][0][2][0] = eri_vref_flat_d_s_d_s;
    funcs_vref_flat[2][1][0][0] = eri_vref_flat_d_p_s_s;
    funcs_vref_flat[2][1][1][0] = eri_vref_flat_d_p_p_s;
    funcs_vref_flat[2][1][1][1] = eri_vref_flat_d_p_p_p;
    funcs_vref_flat[2][1][2][0] = eri_vref_flat_d_p_d_s;
    funcs_vref_flat[2][1][2][1] = eri_vref_flat_d_p_d_p;
    funcs_vref_flat[2][2][0][0] = eri_vref_flat_d_d_s_s;
    funcs_vref_flat[2][2][1][0] = eri_vref_flat_d_d_p_s;
    funcs_vref_flat[2][2][1][1] = eri_vref_flat_d_d_p_p;
    funcs_vref_flat[2][2][2][0] = eri_vref_flat_d_d_d_s;
    funcs_vref_flat[2][2][2][1] = eri_vref_flat_d_d_d_p;
    funcs_vref_flat[2][2][2][2] = eri_vref_flat_d_d_d_d;

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
    // doesn't do anything at the moment
    Boys_Init(0, 7); // need F0 + 7 for interpolation

    // read in the shell info
    // angular momentum and normalization will be handled later
    std::array<int, 4> initam{0, 0, 0, 0};
    VecQuartet gshells_orig(  ReadQuartets(initam, files, false) );

    std::array<int, 4> nshell{gshells_orig[0].size(), gshells_orig[1].size(), 
                              gshells_orig[2].size(), gshells_orig[3].size()};



    /* Storage of integrals results */
    const int maxncart = pow(NCART(MAXAM), 4);
    const int nshell1234 = nshell[0] * nshell[1] * nshell[2] * nshell[3];
    double * res_FO              = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    double * res_vref            = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    double * res_vref_flat       = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    //double * res_liberd          = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));
    double * res_valeev          = (double *)ALLOC(maxncart * nshell1234 * sizeof(double));


    printf("\n");
    printf("%17s    %20s  %20s  %20s    %20s  %20s  %20s\n", "Quartet", "FO MaxErr",    "vref MaxErr",    "vrefF MaxErr", 
                                                                        "FO MaxRelErr", "vref MaxRelErr", "vrefF MaxRelErr");

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

        // copy the shells, set all the am, and renormalize
        VecQuartet gshells(CopyQuartets(gshells_orig));

        for(int m = 0; m < 4; m++)
        {
            for(auto & it : gshells[m])
                it.am = am[m];
            normalize_gaussian_shells(gshells[m].size(), gshells[m].data());
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
        funcs_FO[am[0]][am[1]][am[2]][am[3]](P, Q, res_FO);
        funcs_vref[am[0]][am[1]][am[2]][am[3]](P, Q, res_vref);
        funcs_vref_flat[am[0]][am[1]][am[2]][am[3]](Pf, Qf, res_vref_flat);

        // test with valeev & liberd
        //ValeevIntegrals(gshells, res_valeev);
        ReadValeevIntegrals(basedir, am, res_valeev);
        //ERDIntegrals(gshells, res_liberd);


        //printf("( %d %d | %d %d )\n", am[0], am[1], am[2], am[3]);
        //printf("%22s %22s %22s\n", "liberd", "FO", "valeev");

        double maxerr_FO = 0;
        double maxerr_vref = 0;
        double maxerr_vref_flat = 0;

        double maxrelerr_FO = 0;
        double maxrelerr_vref = 0;
        double maxrelerr_vref_flat = 0;

        int idx = 0;
        for(int i = 0; i < nshell1234; i++)
        {
            for(int j = 0; j < ncart; j++)
            {
                //printf("%22.4e  %22.4e  %22.4e\n", res_liberd[idx], res_FO[idx], res_valeev[idx]);

                const double v = res_valeev[idx];
                //double diff_liberd  = fabs(res_liberd[idx]         - v);
                double diff_FO          = fabs(res_FO[idx]     - v);
                double diff_vref        = fabs(res_vref[idx]     - v);
                double diff_vref_flat   = fabs(res_vref_flat[idx]     - v);
                double rel_FO           = fabs(diff_FO / v);
                double rel_vref         = fabs(diff_vref / v);
                double rel_vref_flat    = fabs(diff_vref_flat / v);
                //printf("%22.4e  %22.4e\n", diff_liberd, diff_FO);
                //printf("\n");

                if(diff_FO > maxerr_FO)
                    maxerr_FO = diff_FO;
                if(rel_FO > maxrelerr_FO)
                    maxrelerr_FO = rel_FO;

                if(diff_vref > maxerr_vref)
                    maxerr_vref = diff_vref;
                if(rel_vref > maxrelerr_vref)
                    maxrelerr_vref = rel_vref;

                if(diff_vref_flat > maxerr_vref_flat)
                    maxerr_vref_flat = diff_vref_flat;
                if(rel_vref_flat > maxrelerr_vref_flat)
                    maxrelerr_vref_flat = rel_vref_flat;

                idx++;
            }
            //printf("\n");
        }

        printf("( %2d %2d | %2d %2d )    %20.8e  %20.8e  %20.8e    %20.8e  %20.8e  %20.8e\n", am[0], am[1], am[2], am[3],
                                                      maxerr_FO, maxerr_vref, maxerr_vref_flat,
                                                      maxrelerr_FO, maxrelerr_vref, maxrelerr_vref_flat);

        free_multishell_pair(P);
        free_multishell_pair(Q);
        free_multishell_pair_flat(Pf);
        free_multishell_pair_flat(Qf);
        FreeQuartets(gshells);

    }

    printf("\n");


    FreeQuartets(gshells_orig);

    Valeev_Finalize();
    Boys_Finalize();

    FREE(res_valeev);
    //FREE(res_liberd);
    FREE(res_vref);
    FREE(res_vref_flat);
    FREE(res_FO);

    return 0;
}
