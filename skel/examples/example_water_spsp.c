#include <stdio.h>

#include "simint/simint.h"

/* Quick and dirty example of using Simint
 *
 * We use hard-coded water/STO-3G as a simple example
 */
int main(int argc, char ** argv)
{
    ////////////////////////////////
    // General workflow as follows
    ////////////////////////////////


    ////////////////////////////////////////
    // 1. Initialize the library
    simint_init();


    ////////////////////////////////////////
    // 2. Create your shells
    struct simint_shell s_shells[4];
    struct simint_shell p_shells[1];

    // initialize the shells
    simint_initialize_shell(&s_shells[0]);
    simint_initialize_shell(&s_shells[1]);
    simint_initialize_shell(&s_shells[2]);
    simint_initialize_shell(&s_shells[3]);
    simint_initialize_shell(&p_shells[0]);


    // allocate the memory
    // arguments to simint_allocate_shell
    //    Number of primitives
    //    Pointer to the shell
    // Note that this doesn't actually set the
    // nprim member of the structure
    simint_allocate_shell(3, &s_shells[0]);
    simint_allocate_shell(3, &s_shells[1]);
    simint_allocate_shell(3, &s_shells[2]);
    simint_allocate_shell(3, &s_shells[3]);
    simint_allocate_shell(3, &p_shells[0]);

    //
    // Coordinates are in atomic units!
    //
    // oxygen
    s_shells[0].am = 0;  s_shells[0].nprim = 1;
    s_shells[0].x = 0.00000; s_shells[0].y = 0.0; s_shells[0].z = -0.02;
    s_shells[0].alpha[0] = 1.;  s_shells[0].coef[0] =  1.;

    s_shells[1].am = 0;  s_shells[1].nprim = 1;
    s_shells[1].x = 0.00000; s_shells[1].y = 0.0; s_shells[1].z = -0.02;
    s_shells[1].alpha[0] =   0.4;  s_shells[1].coef[0] = 1.;

    p_shells[0].am = 1;  p_shells[0].nprim = 1;
    p_shells[0].x = 0.00000; p_shells[0].y = 0.0; p_shells[0].z = -0.02;
    p_shells[0].alpha[0] =   0.4;  p_shells[0].coef[0] =  1.;

    // hydrogen 1
    s_shells[2].am = 0;  s_shells[2].nprim = 1;
    s_shells[2].x = -0.74; s_shells[2].y = 0.0; s_shells[2].z = -0.76;
    s_shells[2].alpha[0] =  0.5;  s_shells[2].coef[0] =  1.;

    // hydrogen 2
    s_shells[3].am = 0;  s_shells[3].nprim = 1;
    s_shells[3].x =  0.74; s_shells[3].y = 0.0; s_shells[3].z = -0.76;
    s_shells[3].alpha[0] = 0.5;  s_shells[3].coef[0] =  1.;


    ////////////////////////////////////////
    // 3. Normalize
    // Arguments to simint_normalize_shells:
    //   number of shells
    //   pointer to array of shells
    simint_normalize_shells(4, s_shells);
    simint_normalize_shells(1, p_shells);


    ////////////////////////////////////////
    // 4. Create your multishell pairs
    // Could obviously be done on the fly if needed
    struct simint_multi_shellpair s1p_pair;
    struct simint_multi_shellpair s2p_pair;
    simint_initialize_multi_shellpair(&s1p_pair);
    simint_initialize_multi_shellpair(&s2p_pair);

    // Last argument - set to one of the following to enable primitive screening
    //       1 - Schwarz screening
    //       2 - Fast Schwarz screening
    simint_create_multi_shellpair(1, &s_shells[3], 1, p_shells, &s2p_pair, 0);
    simint_create_multi_shellpair(1, &s_shells[2], 1, p_shells, &s1p_pair, 0);


    ////////////////////////////////////////
    // 5. Compute your integrals
    // For this example, I'll just allocate
    // space for all (although this overestimates
    // since we are exploiting some symmetry)

    // Allocate some required workspace
    // Arguments to simint_ostei_workmem: Derivative order, max am
    double * work = SIMINT_ALLOC(simint_ostei_workmem(0, 1));

    // Argument is the number of doubles
    double * targets = SIMINT_ALLOC(7*7*7*7 * sizeof(double));
    int ncomputed = 0;
    int ntotal = 0;

    // 3rd argument is the screening tolerance. Set to 0.0 to disable
    ncomputed = simint_compute_eri(&s2p_pair, &s1p_pair, 0, work, targets+ntotal);
    printf("Computed %d contracted s2ps1p integrals\n", ncomputed);
    ntotal += ncomputed;
    ntotal = 9;

    printf("++Computed %d contracted integrals\n", ntotal);

    for(int i = 0; i < ntotal; i++)
        printf(" %d  %20.6f\n", i+1, targets[i]);

    ////////////////////////////////////////
    // 6. Clean up and finalize the library
    //
    simint_free_shell(&s_shells[0]);
    simint_free_shell(&s_shells[1]);
    simint_free_shell(&s_shells[2]);
    simint_free_shell(&s_shells[3]);
    simint_free_shell(&p_shells[0]);
    simint_free_multi_shellpair(&s1p_pair);
    simint_free_multi_shellpair(&s2p_pair);
    SIMINT_FREE(targets);
    SIMINT_FREE(work);
    simint_finalize();

    return 0;
}
