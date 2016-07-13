#include <cstdio>
#include <atomic>
#include <memory>
#include <cmath>

#include <omp.h>

#include "test/Common.hpp"
#include "test/ValeevRef.hpp"
#include "test/Simint.hpp"

#ifdef TESTS_ENABLE_LIBERD
  #include "test/ERD.hpp"
#endif

#ifdef TESTS_ENABLE_LIBINT2
  #include "test/Libint2.hpp"
#endif


typedef std::array<int, 4> QAM;
typedef std::pair<double, double> ErrorPair;

// Prog -> { QAM -> { Absolute error, Relative Error } }
typedef std::map<std::string, std::map<QAM, ErrorPair>> ErrorMap;

static void UpdateErrorMap(ErrorMap & m, const std::string & prog, QAM am, ErrorPair p)
{
    ErrorPair val = m.at(prog).at(am);
    val.first = std::max(val.first, p.first);
    val.second = std::max(val.second, p.second);
    m[prog][am] = val;
}


int main(int argc, char ** argv)
{
    // set up the function pointers
    Simint_Init();

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


    // copy and normalize for liberd
    // (MUST BE DONE BEFORE NORMALIZING THE ORIGINAL)
    #ifdef TESTS_ENABLE_LIBERD
    ShellMap shellmap_erd = CopyShellMap(shellmap);
    for(auto & it : shellmap_erd)
        simint_normalize_shells_erd(it.second.size(), it.second.data());
    #endif


    // normalize the original
    for(auto & it : shellmap)
        simint_normalize_shells(it.second.size(), it.second.data());


    // find the max dimensions
    std::array<int, 3> maxparams = FindMapMaxParams(shellmap);
    const int maxam = (maxparams[0] > SIMINT_ERI_MAXAM ? SIMINT_ERI_MAXAM : maxparams[0]);
    const int maxnprim = maxparams[1];
    const int maxsize = maxparams[2];

    // get the number of threads
    int nthread = omp_get_max_threads();


    /* contracted workspace */
    double * all_simint_work = (double *)ALLOC(nthread * SIMINT_ERI_MAX_WORKMEM);


    /* Storage of integrals */
    double * all_res_simint = (double *)ALLOC(nthread * maxsize * sizeof(double));
    double * all_res_valeev = (double *)ALLOC(nthread * maxsize * sizeof(double));


    // Map containing the errors
    // QAM -> (absolute error, relative error)
    ErrorMap errmap_simint;


    // initialize stuff
    ValeevRef_Init();

    // Map of errors
    ErrorMap errors;

    // Create an entry in the error map with no errors
    std::map<QAM, ErrorPair> initerror;
    // Set initial errors to zero
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
        initerror[{i,j,k,l}] = {0,0};


    // make one for simint
    errors.emplace("Simint", initerror);

    ///////////////////////////////////////////////
    // initialize and allocate space for liberd
    ///////////////////////////////////////////////
    #ifdef TESTS_ENABLE_LIBERD
    double * all_res_liberd = (double *)ALLOC(nthread * maxsize * sizeof(double));

    std::vector<std::unique_ptr<ERD_ERI>> erd;

    for(int i = 0; i < nthread; i++)
        erd.push_back(std::unique_ptr<ERD_ERI>(new ERD_ERI(maxam, maxnprim, 1)));

    errors.emplace("LibERD", initerror);
    #endif

    ///////////////////////////////////////////////
    // initialize and allocate space for libint2
    ///////////////////////////////////////////////
    #ifdef TESTS_ENABLE_LIBINT2
    double * all_res_libint2 = (double *)ALLOC(nthread * maxsize * sizeof(double));
    LIBINT2_PREFIXED_NAME(libint2_static_init)();

    std::vector<std::unique_ptr<Libint2_ERI>> libint;

    for(int i = 0; i < nthread; i++)
        libint.push_back(std::unique_ptr<Libint2_ERI>(new Libint2_ERI(maxam, maxnprim)));

    errors.emplace("Libint2", initerror);
    #endif


    // Print the header for the final results table
    printf("\n");
    printf("%17s  ", "Quartet");
    for(size_t i = 0; i < errors.size(); i++) printf("  %10s", "MaxErr");
    for(size_t i = 0; i < errors.size(); i++) printf("  %10s", "MaxRelErr");
    printf("\n");

    printf("%17s  ", "");

    for(const auto & it : errors)
        printf("  %10s", it.first.c_str());
    for(const auto & it : errors)
        printf("  %10s", it.first.c_str());
    printf("\n");


    // Number of contracted integrals calculated
    std::atomic<long> ncont(0);


    // loop over all AM quartets, choosing only valid ones
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        if(!ValidQuartet(i, j, k, l))
            continue;


        const int nshell3 = shellmap[k].size();
        const int nshell4 = shellmap[l].size();

        // create a multishell pair for all | k l ) quartets
        struct simint_multi_shellpair Q = simint_create_multi_shellpair(nshell3, shellmap[k].data(),
                                                                        nshell4, shellmap[l].data());

        // do bra one at a time
        #pragma omp parallel for
        for(size_t a = 0; a < shellmap[i].size(); a++)
        for(size_t b = 0; b < shellmap[j].size(); b++)
        {
            int ithread = omp_get_thread_num();

            double * res_simint = all_res_simint + ithread * maxsize;
            double * simint_work = all_simint_work + ithread * SIMINT_ERI_MAX_WORKSIZE;
            double * res_valeev = all_res_valeev + ithread * maxsize;

            #ifdef TESTS_ENABLE_LIBERD
            double * res_liberd = all_res_liberd + ithread * maxsize;
            #endif
            #ifdef TESTS_ENABLE_LIBINT2
            double * res_libint2 = all_res_libint2 + ithread * maxsize;
            #endif

            const int nshell1 = 1;
            const int nshell2 = 1;

            struct simint_multi_shellpair P = simint_create_multi_shellpair(nshell1, &shellmap[i][a],
                                                                            nshell2, &shellmap[j][b]);


            const int ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);

            const int nshell1234 = P.nshell12 * Q.nshell12;
            const int arrlen = nshell1234 * ncart1234;


            /////////////////////////////////
            // Calculate valeev
            // reference integrals
            /////////////////////////////////
            ValeevRef_Integrals(&shellmap[i][a], nshell1,
                                &shellmap[j][b], nshell2,
                                &shellmap[k][0], nshell3,
                                &shellmap[l][0], nshell4,
                                res_valeev, false);



            ////////////////////////////
            // Calculate with my code
            ////////////////////////////
            Simint_Integral(P, Q, simint_work, res_simint);


            /////////////////////////////////
            // Calculate with liberd
            /////////////////////////////////
            #ifdef TESTS_ENABLE_LIBERD
            erd.at(ithread)->Integrals(&shellmap_erd[i][a], nshell1,
                          &shellmap_erd[j][b], nshell2,
                          &shellmap_erd[k][0], nshell3,
                          &shellmap_erd[l][0], nshell4,
                          res_liberd);
            #endif


            /////////////////////////////////
            // Calculate with libint2
            /////////////////////////////////
            #ifdef TESTS_ENABLE_LIBINT2
            libint.at(ithread)->Integrals(P, Q, res_libint2);
            #endif


            /////////////////////////////////
            // Update the error map
            /////////////////////////////////
            #pragma omp critical
            {
                UpdateErrorMap(errors, "Simint", {i, j, k, l}, CalcError(res_simint, res_valeev, arrlen));

                #ifdef TESTS_ENABLE_LIBERD
                UpdateErrorMap(errors, "LibERD", {i, j, k, l}, CalcError(res_liberd, res_valeev, arrlen));
                #endif
                #ifdef TESTS_ENABLE_LIBINT2
                UpdateErrorMap(errors, "Libint2", {i, j, k, l}, CalcError(res_libint2, res_valeev, arrlen));
                #endif
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
                double diff_simint = fabs(res_valeev[m] - res_simint[m]);
                double rdiff_simint = fabs(diff_simint / res_valeev[m]);

                double diff_liberd = 0;
                double rdiff_liberd = 0;
                double diff_libint2 = 0;
                double rdiff_libint2 = 0;

                #ifdef TESTS_ENABLE_LIBERD
                diff_liberd = fabs(res_valeev[m] - res_liberd[m]);
                rdiff_liberd = fabs(diff_liberd / res_liberd[m]);
                #endif

                #ifdef TESTS_ENABLE_LIBINT2
                diff_libint2 = fabs(res_valeev[m] - res_libint2[m]);
                rdiff_libint2 = fabs(diff_libint2 / res_valeev[m]);
                #endif


                if( (diff_simint > 1e-14 && rdiff_simint > 1e-8) ||
                    (diff_liberd > 1e-14 && rdiff_liberd > 1e-8) ||
                    (diff_libint2 > 1e-14 && rdiff_libint2 > 1e-8) )
                {
                    printf(" [%d/%d] %d %d %d %d : %d %d %d %d", ithread, nthread, int(a+m1), int(b+m2), m3, m4, c1, c2, c3, c4);

                    printf("   %25.16e  %25.16e  %25.16e", res_simint[m], diff_simint, rdiff_simint);
                    #ifdef TESTS_ENABLE_LIBERD
                    printf("   %25.16e  %25.16e  %25.16e", res_liberd[m], diff_liberd, rdiff_liberd);
                    #endif
                    #ifdef TESTS_ENABLE_LIBINT2
                    printf("   %25.16e  %25.16e  %25.16e", res_libint2[m], diff_libint2, rdiff_libint2);
                    #endif
                    printf("\n");
                }
            }


            simint_free_multi_shellpair(P);

            ncont += arrlen;

        } // end threaded loop over a,b


        simint_free_multi_shellpair(Q);

        // print out the errors for this quartet
        printf("( %2d %2d | %2d %2d )  ", i, j, k, l);

        // should be the same order as the header, right?
        for(const auto & it : errors)
            printf("  %10.3e", it.second.at({i, j, k, l}).first);
        for(const auto & it : errors)
            printf("  %10.3e", it.second.at({i, j, k, l}).second);
        printf("\n");
    }

    printf("\n");
    printf("Calculated %ld contracted integrals\n", static_cast<long>(ncont));
    printf("\n");

    FreeShellMap(shellmap);

    #ifdef TESTS_ENABLE_LIBERD
    FreeShellMap(shellmap_erd);
    FREE(all_res_liberd);
    #endif

    #ifdef TESTS_ENABLE_LIBINT2
    FREE(all_res_libint2);
    #endif

    ValeevRef_Finalize();
    Simint_Finalize();

    FREE(all_res_simint);
    FREE(all_simint_work);
    FREE(all_res_valeev);

    return 0;
}
