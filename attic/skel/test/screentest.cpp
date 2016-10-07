#include "test/Common.hpp"

#include "simint/eri/eri.h"
#include "simint/simint_init.h"
#include "simint/shell/shell_screen.h"

#define SIMINT_SCREEN 1
#define SIMINT_SCREEN_TOL 1e-13
#define SIMINT_SCREEN_TOL2 (SIMINT_SCREEN_TOL * SIMINT_SCREEN_TOL)

typedef std::pair<int, int> IntPair;


int main(int argc, char ** argv)
{
    // set up the function pointers
    simint_init();

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

    // normalize
    for(auto & it : shellmap)
        simint_normalize_shells(it.second.size(), it.second.data());

    // find the max dimensions
    std::pair<int, int> maxparams = FindMaxParams(shellmap);
    const int maxam = (maxparams.first > SIMINT_ERI_MAXAM ? SIMINT_ERI_MAXAM : maxparams.first);


    // expand the shell map
    GaussianVec shellvec;
    for(auto & it : shellmap)
        shellvec.insert(shellvec.end(), it.second.begin(), it.second.end());

    // maximum number of primitives
    int max_nprim = 0;
    for(const auto & it : shellvec)
        max_nprim = std::max(max_nprim, it.nprim);

    const size_t nshell = shellvec.size();
    const size_t npairtotal = (nshell * (nshell+1))/2;


    // find the maximum screening value for all shell pair
    // and all primitive pair
    double shell_screen_max = 0.0;
    double prim_screen_max = 0.0;

    for(size_t i = 0; i < shellvec.size(); i++)
    for(size_t j = 0; j <= i; j++)
    {
        const simint_shell * A = &shellvec[i];
        const simint_shell * B = &shellvec[j];

        double shell_screen = simint_shellscreen_schwarz_max(A, B); 
        double prim_screen = simint_primscreen_schwarz_max(A, B, NULL); 

        shell_screen_max = std::max(shell_screen_max, shell_screen);
        prim_screen_max = std::max(prim_screen_max, prim_screen);
    }

    printf("\nTotal number of shells: %lu\n", shellvec.size());
    printf("Total number of shell pair: %lu\n", npairtotal);
    printf("Maximum screening value for all contracted shell pair: %12.6e\n", shell_screen_max);
    printf("Maximum screening value for all primitive shell pair : %12.6e\n", prim_screen_max);


    /////////////////////////////////////////////
    // CONTRACTED SHELL PAIRS
    /////////////////////////////////////////////
    // how many shell pair are never significant
    size_t npair_insig = 0;
    size_t npair_sig = 0;
    std::vector<IntPair> sig_shellpair;

    for(size_t i = 0; i < shellvec.size(); i++)
    for(size_t j = 0; j <= i; j++)
    {
        const simint_shell * A = &shellvec[i];
        const simint_shell * B = &shellvec[j];
        const double screen = simint_shellscreen_schwarz_max(A, B);
        const double val4 = screen * shell_screen_max;

        if(val4 < SIMINT_SCREEN_TOL2)
            npair_insig++;
        else
        {
            npair_sig++;
            sig_shellpair.push_back({i, j});
        }
    }

    if( (npair_insig + npair_sig) != npairtotal)
        printf("\n\nWARNING WARNING: Bad logic somewhere\n");

    double frac = ((double)(npair_sig))/((double)(npairtotal));
    printf("\n\n");
    printf("=== CONTRACTED SHELLS ===\n");
    printf("         Total number of shell pair: %lu\n", npairtotal);
    printf("         Number of significant pair: %lu (%6.2f%% of total)\n", sig_shellpair.size(), 100.0*frac);
    printf("Number of always-insignificant pair: %lu (%6.2f%% of total)\n", npair_insig, 100.0*(1.0-frac));
        

    /////////////////////////////////////////////
    // PRIMITIVE SHELL PAIRS
    /////////////////////////////////////////////
    double * primscreen_info = (double *)ALLOC(max_nprim * max_nprim * sizeof(double));
    size_t nprimpair_insig = 0;
    size_t nprimpair_sig = 0;

    for(size_t i = 0; i < shellvec.size(); i++)
    for(size_t j = 0; j <= i; j++)
    {
        const simint_shell * A = &shellvec[i];
        const simint_shell * B = &shellvec[j];
        const int nprimA = A->nprim;
        const int nprimB = B->nprim;

        const size_t nprim12 = nprimA * nprimB;
        const size_t nprimtri = (i == j ? ((nprimA)*(nprimA+1))/2 : nprim12);

        const double screen = simint_primscreen_schwarz_max(A, B, primscreen_info);

        for(size_t p = 0; p < nprimtri; p++)
        {
            const double val4 = primscreen_info[p] * prim_screen_max;

            if(val4 < SIMINT_SCREEN_TOL2)
                nprimpair_insig++;
            else
                nprimpair_sig++;
        }
    }

    FREE(primscreen_info);

    const size_t nprimpair_total = nprimpair_insig + nprimpair_sig;
    frac = ((double)(nprimpair_sig))/((double)(nprimpair_total));
    printf("\n\n");
    printf("=== PRIMITIVE SHELLS ===\n");
    printf("               Number of shell pair: %lu\n", nprimpair_total);
    printf("         Number of significant pair: %lu (%6.2f%% of total)\n", nprimpair_sig, 100.0*frac);
    printf("Number of always-insignificant pair: %lu (%6.2f%% of total)\n", nprimpair_insig, 100.0*(1.0-frac));


    //printf("\n");
    //printf("%17s  %15s  %15s  %15s\n", "Quartet", "Nshell", "Nprim", "Ncont");

        // total number of shells, primitives, etc
        //const size_t nshell1234 = P.nshell12 * Q.nshell12;
        //const size_t Pnprim = P.nprim;
        //const size_t Qnprim = Q.nprim;
        //const size_t nprim1234 = Pnprim * Qnprim;
        //const size_t ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);
        //const size_t ncont1234 = nshell1234 * ncart1234;



        //printf("( %2d %2d | %2d %2d )  %15lu  %15lu  %15lu\n",
        //       i, j, k, l, nshell1234, nprim1234, ncont1234); 


    //for(auto & it : pair_map)
    //    simint_free_multi_shellpair(&it.second);

    FreeShellMap(shellmap);
    simint_finalize();

    printf("\n");
    return 0;
}
