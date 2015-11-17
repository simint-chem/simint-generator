#include <cstdio>
#include <memory>

#include <omp.h>

#include "boys/boys.h"
#include "test/common.hpp"
#include "test/valeev.hpp"
#include "test/ERD.hpp"
#include "test/Libint2.hpp"
#include "vectorization/vectorization.h"

typedef std::array<int, 4> QAM;


static void UpdateMap(std::map<QAM, std::pair<double, double>> & m, std::pair<double, double> p, QAM am)
{
    std::pair<double, double> val = m.at(am);
    val.first = std::max(val.first, p.first);
    val.second = std::max(val.second, p.second);
    m[am] = val;
}


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
    const int maxam = (maxparams[0] > MAXAM ? MAXAM : maxparams[0]);
    const int maxnprim = maxparams[1];
    const int maxsize = maxparams[2];

    // get the number of threads
    int nthread = omp_get_max_threads();


    /* Storage of integrals */
    double * all_res                 = (double *)ALLOC(nthread * maxsize * sizeof(double));
    double * all_res_liberd          = (double *)ALLOC(nthread * maxsize * sizeof(double));
    double * all_res_libint          = (double *)ALLOC(nthread * maxsize * sizeof(double));
    double * all_res_valeev          = (double *)ALLOC(nthread * maxsize * sizeof(double));

    // initialize stuff
    Valeev_Init();
    Boys_Init();
    LIBINT2_PREFIXED_NAME(libint2_static_init)();

    std::vector<std::unique_ptr<ERD_ERI>> erd;
    std::vector<std::unique_ptr<Libint2_ERI>> libint;

    for(int i = 0; i < nthread; i++)
    {
        erd.push_back(std::unique_ptr<ERD_ERI>(new ERD_ERI(maxam, maxnprim, 1)));
        libint.push_back(std::unique_ptr<Libint2_ERI>(new Libint2_ERI(maxam, maxnprim)));
    }
             
             
    printf("\n");
    printf("%17s    %10s  %10s  %10s   %10s  %10s  %10s\n", 
                                                         "Quartet", "MaxErr", "", "", 
                                                         "MaxRelErr", "", "");
    printf("%17s    %10s  %10s  %10s   %10s  %10s  %10s\n", "",
                                                         "Me", "erd", "libint", 
                                                         "Me", "erd", "libint");

    std::map<QAM, std::pair<double, double>> errmap, errmap_erd, errmap_libint;

    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        errmap[{i,j,k,l}] = {0,0};
        errmap_erd[{i,j,k,l}] = {0,0};
        errmap_libint[{i,j,k,l}] = {0,0};
    }

        

    // loop over all quartets, choosing only valid ones
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    {
        // do one shell pair at a time on the bra side
        #pragma omp parallel for
        for(size_t a = 0; a < shellmap[i].size(); a++)
        for(size_t b = 0; b < shellmap[j].size(); b++)
        {
            int ithread = omp_get_thread_num();
            double * res           = all_res + ithread * maxsize;
            double * res_liberd    = all_res_liberd + ithread * maxsize;
            double * res_libint    = all_res_libint + ithread * maxsize;
            double * res_valeev    = all_res_valeev + ithread * maxsize;

            const int nshell1 = 1; 
            const int nshell2 = 1; 
    
            struct multishell_pair P = create_multishell_pair(1, &shellmap[i][a],
                                                              1, &shellmap[j][b]);
       


            for(int k = 0; k <= maxam; k++)
            for(int l = 0; l <= maxam; l++)
            {
                if(!ValidQuartet(i, j, k, l))
                    continue;

                const int ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);

                const int nshell3 = shellmap[k].size();
                const int nshell4 = shellmap[l].size();

                const int nshell1234 = nshell1 * nshell2 * nshell3 * nshell4;
                const int arrlen = nshell1234 * ncart1234;


                /////////////////////////////////
                // Calculate valeev
                // reference integrals
                /////////////////////////////////
                ValeevIntegrals(&shellmap[i][a], nshell1,
                                &shellmap[j][b], nshell2,
                                &shellmap[k][0], nshell3,
                                &shellmap[l][0], nshell4,
                                res_valeev, false);


                struct multishell_pair Q = create_multishell_pair(nshell3, shellmap[k].data(),
                                                                  nshell4, shellmap[l].data());


                /////////////////////////////////
                // Calculate with liberd and libint2
                /////////////////////////////////
                erd.at(ithread)->Integrals(&shellmap_erd[i][a], nshell1,
                              &shellmap_erd[j][b], nshell2,
                              &shellmap_erd[k][0], nshell3,
                              &shellmap_erd[l][0], nshell4,
                              res_liberd);

                //libint.at(ithread)->Integrals(P, Q, res_libint);

                ////////////////////////////
                // Calculate with my code
                ////////////////////////////
                Integral(P, Q, res);



                // Analyze
                std::pair<double, double> err           = CalcError(res,         res_valeev,  arrlen);
                std::pair<double, double> err_erd       = CalcError(res_liberd,  res_valeev,  arrlen);
                //std::pair<double, double> err_libint    = CalcError(res_libint,  res_valeev,  arrlen);

                #pragma omp critical
                {
                    UpdateMap(errmap, err, {i,j,k,l});
                    UpdateMap(errmap_erd, err_erd, {i,j,k,l});
                    //UpdateMap(errmap_libint, err_libint, {i,j,k,l});
                }



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
                    double ediff = fabs(res_valeev[m] - res_liberd[m]);
                    double erdiff = fabs(ediff / res_liberd[m]);

                    //if( (diff > 1e-14 && rdiff > 1e-8) || (ldiff > 1e-14 && lrdiff > 1e-8) )
                    if( (diff > 1e-14 && rdiff > 1e-8) || (ediff > 1e-14 && erdiff > 1e-8) )
                        printf(" [%d/%d] %d %d %d %d : %d %d %d %d  %25.16e  %25.16e  %25.16e   %25.16e  %25.16e  %25.16e  %25.16e\n",
                                                                                                  ithread, nthread,
                                                                                                  m1, m2, m3, m4, c1, c2, c3, c4,
                                                                                                  res_valeev[m], res_liberd[m], res[m], 
                                                                                                  diff, rdiff, ediff, erdiff);
                }


                free_multishell_pair(Q);
                //free_multishell_pair_flat(Pf);
                //free_multishell_pair_flat(Qf);
            }

            free_multishell_pair(P);

        } // end threaded loop over shells

        for(int k = 0; k <= maxam; k++)
        for(int l = 0; l <= maxam; l++)
        {
            if(!ValidQuartet(i,j,k,l))
                continue;

            std::pair<double, double> err        = errmap.at({i,j,k,l});
            std::pair<double, double> err_erd    = errmap_erd.at({i,j,k,l});
            std::pair<double, double> err_libint = errmap_libint.at({i,j,k,l});
            printf("( %2d %2d | %2d %2d )    %10.3e  %10.3e  %10.3e    %10.3e  %10.3e  %10.3e\n", i, j, k, l,
                                                          err.first,  err_erd.first,  err_libint.first,
                                                          err.second, err_erd.second, err_libint.second);
        }

    }

    printf("\n");

    FreeShellMap(shellmap);
    FreeShellMap(shellmap_erd);

    Valeev_Finalize();
    Boys_Finalize();

    FREE(all_res_valeev);
    FREE(all_res_liberd);
    FREE(all_res_libint);
    FREE(all_res);

    return 0;
}
