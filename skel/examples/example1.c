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
    s_shells[0].am = 0;  s_shells[0].nprim = 3;
    s_shells[0].x = 0.00000; s_shells[0].y = -0.14322; s_shells[0].z = 0.0;
    s_shells[0].alpha[0] = 130.7093200;  s_shells[0].coef[0] =  0.15432897;
    s_shells[0].alpha[1] =  23.8088610;  s_shells[0].coef[1] =  0.53532814; 
    s_shells[0].alpha[2] =   6.4436083;  s_shells[0].coef[2] =  0.44463454;

    s_shells[1].am = 0;  s_shells[1].nprim = 3;
    s_shells[1].x = 0.00000; s_shells[1].y = -0.14322; s_shells[1].z = 0.0;
    s_shells[1].alpha[0] =   5.0331513;  s_shells[1].coef[0] = -0.09996723;
    s_shells[1].alpha[1] =   1.1695961;  s_shells[1].coef[1] =  0.39951283; 
    s_shells[1].alpha[2] =   0.3803890;  s_shells[1].coef[2] =  0.70011547;

    p_shells[0].am = 1;  p_shells[0].nprim = 3;
    p_shells[0].x = 0.00000; p_shells[0].y = -0.14322; p_shells[0].z = 0.0;
    p_shells[0].alpha[0] =   5.0331513;  p_shells[0].coef[0] =  0.15591627;
    p_shells[0].alpha[1] =   1.1695961;  p_shells[0].coef[1] =  0.60768372; 
    p_shells[0].alpha[2] =   0.3803890;  p_shells[0].coef[2] =  0.39195739;

    // hydrogen 1
    s_shells[2].am = 0;  s_shells[2].nprim = 3;
    s_shells[2].x =  1.63804; s_shells[2].y = 1.13654; s_shells[2].z = 0.0;
    s_shells[2].alpha[0] =  3.42525091;  s_shells[2].coef[0] =  0.15432897;
    s_shells[2].alpha[1] =  0.62391373;  s_shells[2].coef[1] =  0.53532814; 
    s_shells[2].alpha[2] =  0.16885540;  s_shells[2].coef[2] =  0.44463454;

    // hydrogen 2
    s_shells[3].am = 0;  s_shells[3].nprim = 3;
    s_shells[3].x = -1.63804; s_shells[3].y = 1.13654; s_shells[3].z = 0.0;
    s_shells[3].alpha[0] =  3.42525091;  s_shells[3].coef[0] =  0.15432897;
    s_shells[3].alpha[1] =  0.62391373;  s_shells[3].coef[1] =  0.53532814; 
    s_shells[3].alpha[2] =  0.16885540;  s_shells[3].coef[2] =  0.44463454;


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
    struct simint_multi_shellpair ss_pair;
    struct simint_multi_shellpair ps_pair;
    struct simint_multi_shellpair pp_pair;
    simint_initialize_multi_shellpair(&ss_pair);
    simint_initialize_multi_shellpair(&ps_pair);
    simint_initialize_multi_shellpair(&pp_pair);

    // Last argument - set to one of the following to enable primitive screening
    //       1 - Schwarz screening
    //       2 - Fast Schwarz screening
    simint_create_multi_shellpair(4, s_shells, 4, s_shells, &ss_pair, 0);
    simint_create_multi_shellpair(1, p_shells, 4, s_shells, &ps_pair, 0);
    simint_create_multi_shellpair(1, p_shells, 1, p_shells, &pp_pair, 0);


    ////////////////////////////////////////
    // 5. Compute your integrals
    // For this example, I'll just allocate
    // space for all (although this overestimates
    // since we are exploiting some symmetry)

    // Allocate some required workspace
    double * work = SIMINT_ALLOC(SIMINT_OSTEI_MAX_WORKMEM);

    // Argument is the number of doubles
    double * targets = SIMINT_ALLOC(7*7*7*7 * sizeof(double));
    int ncomputed = 0;
    int ntotal = 0;

    // 3rd argument is the screening tolerance. Set to 0.0 to disable
    ncomputed = simint_compute_eri(&ss_pair, &ss_pair, 0.0, work, targets+ntotal);
    printf("Computed %d contracted ssss integrals\n", ncomputed);
    ntotal += ncomputed;

    ncomputed = simint_compute_eri(&ps_pair, &ss_pair, 0.0, work, targets+ntotal);
    ncomputed *= 3;
    printf("Computed %d contracted psss integrals\n", ncomputed);
    ntotal += ncomputed;

    ncomputed = simint_compute_eri(&ps_pair, &ps_pair, 0.0, work, targets+ntotal);
    ncomputed *= 9;
    printf("Computed %d contracted psps integrals\n", ncomputed);
    ntotal += ncomputed;

    ncomputed = simint_compute_eri(&pp_pair, &ss_pair, 0.0, work, targets+ntotal);
    ncomputed *= 9;
    printf("Computed %d contracted ppss integrals\n", ncomputed);
    ntotal += ncomputed;

    ncomputed = simint_compute_eri(&pp_pair, &ps_pair, 0.0, work, targets+ntotal);
    ncomputed *= 27;
    printf("Computed %d contracted ppps integrals\n", ncomputed);
    ntotal += ncomputed;

    ncomputed = simint_compute_eri(&pp_pair, &pp_pair, 0.0, work, targets+ntotal);
    ncomputed *= 81;
    printf("Computed %d contracted pppp integrals\n", ncomputed);
    ntotal += ncomputed;

    printf("++Computed %d contracted integrals\n", ntotal);

    //for(int i = 0; i < ntotal; i++)
    //    printf(" %d  %12.8e\n", i, targets[i]);

    ////////////////////////////////////////
    // 6. Clean up and finalize the library
    //
    simint_free_shell(&s_shells[0]);
    simint_free_shell(&s_shells[1]);
    simint_free_shell(&s_shells[2]);
    simint_free_shell(&s_shells[3]);
    simint_free_shell(&p_shells[0]);
    simint_free_multi_shellpair(&ss_pair);
    simint_free_multi_shellpair(&ps_pair);
    simint_free_multi_shellpair(&pp_pair);
    SIMINT_FREE(targets);
    SIMINT_FREE(work);
    simint_finalize();

    return 0;
}
