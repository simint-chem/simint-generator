#include <stdexcept>

#include "test/Common.hpp"

#include "simint/simint.h"
#include "simint/shell/shell_screen.h"

#define SIMINT_SCREEN 2
#define SIMINT_SCREEN_TOL 1e-14
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

    // expand the shell map
    GaussianVec shellvec;
    for(auto & it : shellmap)
        shellvec.insert(shellvec.end(), it.second.begin(), it.second.end());

    // create all the shell pair
    std::vector<simint_multi_shellpair> shellpairvec;
    for(size_t i = 0; i < shellvec.size(); i++)
    for(size_t j = 0; j <= i; j++)
    {
        struct simint_multi_shellpair P;
        simint_initialize_multi_shellpair(&P);
        simint_create_multi_shellpair(1, &shellvec[i], 1, &shellvec[j], &P, SIMINT_SCREEN);
        shellpairvec.push_back(P);
    }

    // number of shells
    const size_t nshell = shellvec.size();

    // expected number of shell pair
    const size_t nshellpair_total = (nshell * (nshell+1))/2;

    if(nshellpair_total != shellpairvec.size())
        throw std::logic_error("Logic error in number of shell pair");


    // maximum number of primitives
    int max_nprim = 0;
    for(const auto & it : shellvec)
        max_nprim = std::max(max_nprim, it.nprim);


    // find the maximum screening value for all shell pair
    // and all primitive pair. Also store the values for
    // contracted shell screening
    double shell_screen_max = 0.0;
    double prim_screen_max = 0.0;

    std::vector<double> shellpair_screen;
    shellpair_screen.reserve(nshellpair_total);

    // store the number of primitive pair for each shell pair.
    // Only the padded size is stored in the simint_multi_shellpair
    // structure.
    std::vector<size_t> shellpair_nprim;
    shellpair_nprim.reserve(nshellpair_total);

    for(size_t i = 0; i < shellvec.size(); i++)
    for(size_t j = 0; j <= i; j++)
    {
        const simint_shell * A = &shellvec[i];
        const simint_shell * B = &shellvec[j];
        const int nprimA = A->nprim;
        const int nprimB = B->nprim;
        const size_t nprim12 = nprimA * nprimB;
        const size_t nprimtri = (i == j ? ((nprimA)*(nprimA+1))/2 : nprim12);

        double shell_screen = simint_shellscreen(A, B, SIMINT_SCREEN);
        shell_screen_max = std::max(shell_screen_max, shell_screen);

        shellpair_nprim.push_back(nprimtri);
        shellpair_screen.push_back(shell_screen);
    }

    if( shellpair_screen.size() != nshellpair_total )
        throw std::logic_error("Bad logic in contracted shell pair");
    if( shellpair_nprim.size() != nshellpair_total )
        throw std::logic_error("Bad logic in contracted shell pair");

    // we can get the primitive max from the shell pairs
    for(const auto & it : shellpairvec)
        prim_screen_max = std::max(prim_screen_max, it.screen_max);

    printf("\n");
    printf("Total number of shells: %lu\n", shellvec.size());
    printf("Total number of shell pair: %lu\n", nshellpair_total);
    printf("Maximum screening value for all contracted shell pair: %12.6e\n", shell_screen_max);
    printf("Maximum screening value for all primitive shell pair : %12.6e\n", prim_screen_max);


    /////////////////////////////////////////////
    // CONTRACTED SHELL PAIRS
    /////////////////////////////////////////////
    // how many shell pair are never significant
    // also
    size_t npair_insig = 0;
    size_t npair_sig = 0;

    for(const auto & it : shellpair_screen)
    {
        const double val = it * shell_screen_max;

        if(val < SIMINT_SCREEN_TOL2)
            npair_insig++;
        else
            npair_sig++;
    }

    if( (npair_insig + npair_sig) != nshellpair_total )
        throw std::logic_error("Bad logic in contracted shell pair");

    double frac = ((double)(npair_sig))/((double)(nshellpair_total));
    printf("\n\n");
    printf("=== CONTRACTED SHELL PAIRS ===\n");
    printf("         Total number of shell pair: %lu\n", nshellpair_total);
    printf("         Number of significant pair: %lu (%6.2f%% of total)\n", npair_sig, 100.0*frac);
    printf("Number of always-insignificant pair: %lu (%6.2f%% of total)\n", npair_insig, 100.0*(1.0-frac));


    /////////////////////////////////////////////
    // PRIMITIVE SHELL PAIRS
    /////////////////////////////////////////////
    double * primscreen_info = (double *)SIMINT_ALLOC(max_nprim * max_nprim * sizeof(double));
    size_t nprimpair_insig = 0;
    size_t nprimpair_sig = 0;

    for(size_t ij = 0; ij < shellpairvec.size(); ij++)
    {
        const size_t nprimpair = shellpair_nprim[ij];

        // check to make sure our indexing is ok (see note above)
        // the nprim member of the shell pair struct is padded, and we
        // don't want to include those values
        if(SIMINT_SIMD_ROUND(nprimpair) != static_cast<size_t>(shellpairvec[ij].nprim))
            throw std::logic_error("Bad number of primitives stored in the shell pair");

        for(size_t p = 0; p < nprimpair; p++)
        {
            const double val = shellpairvec[ij].screen[p] * prim_screen_max;
            if(val < SIMINT_SCREEN_TOL2)
                nprimpair_insig++;
            else
                nprimpair_sig++;
        }
    }


    SIMINT_FREE(primscreen_info);

    const size_t nprimpair_total = nprimpair_insig + nprimpair_sig;
    frac = ((double)(nprimpair_sig))/((double)(nprimpair_total));
    printf("\n\n");
    printf("=== PRIMITIVE PAIR ===\n");
    printf("           Number of primitive pair: %lu\n", nprimpair_total);
    printf("         Number of significant pair: %lu (%6.2f%% of total)\n", nprimpair_sig, 100.0*frac);
    printf("Number of always-insignificant pair: %lu (%6.2f%% of total)\n", nprimpair_insig, 100.0*(1.0-frac));



    /////////////////////////////////////////////
    // CONTRACTED SHELL QUARTETS
    /////////////////////////////////////////////
    size_t nshellquartet_insig = 0;
    size_t nshellquartet_sig = 0;

    for(size_t ij = 0; ij < shellpairvec.size(); ij++)
    for(size_t kl = 0; kl <= ij; kl++)
    {
        const double val = shellpair_screen[ij] * shellpair_screen[kl];
        if(val < SIMINT_SCREEN_TOL2)
            nshellquartet_insig++;
        else
            nshellquartet_sig++;
    }

    const size_t nshellquartet_total = nshellquartet_insig + nshellquartet_sig;

    if(nshellquartet_total != (nshellpair_total * (nshellpair_total+1))/2)
        throw std::logic_error("Bad number of contracted shell quartets");

    frac = ((double)(nshellquartet_sig))/((double)(nshellquartet_total));
    printf("\n\n");
    printf("=== CONTRACTED SHELL QUARTETS ===\n");
    printf("               Number of shell quartets: %lu\n", nshellquartet_total);
    printf("         Number of significant quartets: %lu (%6.2f%% of total)\n", nshellquartet_sig, 100.0*frac);
    printf("       Number of insignificant quartets: %lu (%6.2f%% of total)\n", nshellquartet_insig, 100.0*(1.0-frac));


    /////////////////////////////////////////////
    // PRIMITIVE QUARTETS
    /////////////////////////////////////////////
    size_t nprimquartet_insig = 0;
    size_t nprimquartet_sig = 0;

    for(size_t ij = 0; ij < shellpairvec.size(); ij++)
    for(size_t kl = 0; kl <= ij; kl++)
    {
        const size_t nprimpair1 = shellpair_nprim[ij];
        const size_t nprimpair2 = shellpair_nprim[kl];

        for(size_t p = 0; p < nprimpair1; p++)
        {
            size_t qend = (ij == kl ? (p+1) : nprimpair2);
            for(size_t q = 0; q < qend; q++)
            {
                const double val = shellpairvec[ij].screen[p] * shellpairvec[kl].screen[q];
                if(val < SIMINT_SCREEN_TOL2)
                    nprimquartet_insig++;
                else
                    nprimquartet_sig++;
            }
        }
    }

    const size_t nprimquartet_total = nprimquartet_insig + nprimquartet_sig;

    if(nprimquartet_total != (nprimpair_total * (nprimpair_total+1))/2)
        throw std::logic_error("Bad number of contracted prim quartets");

    frac = ((double)(nprimquartet_sig))/((double)(nprimquartet_total));
    printf("\n\n");
    printf("=== PRIMITIVE QUARTETS ===\n");
    printf("           Number of primitive quartets: %lu\n", nprimquartet_total);
    printf("         Number of significant quartets: %lu (%6.2f%% of total)\n", nprimquartet_sig, 100.0*frac);
    printf("       Number of insignificant quartets: %lu (%6.2f%% of total)\n", nprimquartet_insig, 100.0*(1.0-frac));



    FreeShellMap(shellmap);
    for(auto & it : shellpairvec)
        simint_free_multi_shellpair(&it);

    simint_finalize();

    printf("\n");
    return 0;
}
