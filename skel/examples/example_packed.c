#include <stdio.h>

#include "simint/simint.h"

/* Quick and dirty example of using Simint
 *
 * We use hard-coded water/STO-3G as a simple example
 */
int main(int argc, char ** argv)
{
    //////////////////////////////////////////////////////
    // Simple example with calculation + packed storage
    //////////////////////////////////////////////////////


    ////////////////////////////////////////
    // 1. Initialize the library
    simint_init();


    ////////////////////////////////////////
    // 2. Create your shells
    struct simint_shell shells[5];

    // initialize the shells
    for(int i = 0; i < 5; i++)
        simint_initialize_shell(&shells[i]);


    // allocate the memory
    // arguments to simint_allocate_shell
    //    Number of primitives
    //    Pointer to the shell
    // Note that this doesn't actually set the
    // nprim member of the structure
    simint_allocate_shell(3, &shells[0]);
    simint_allocate_shell(3, &shells[1]);
    simint_allocate_shell(3, &shells[2]);
    simint_allocate_shell(3, &shells[3]);
    simint_allocate_shell(3, &shells[4]);

    //
    // Coordinates are in atomic units!
    //
    // oxygen
    shells[0].am = 0;  shells[0].nprim = 3;
    shells[0].x = 0.00000; shells[0].y = -0.14322; shells[0].z = 0.0;
    shells[0].alpha[0] = 130.7093200;  shells[0].coef[0] =  0.15432897;
    shells[0].alpha[1] =  23.8088610;  shells[0].coef[1] =  0.53532814; 
    shells[0].alpha[2] =   6.4436083;  shells[0].coef[2] =  0.44463454;

    shells[1].am = 0;  shells[1].nprim = 3;
    shells[1].x = 0.00000; shells[1].y = -0.14322; shells[1].z = 0.0;
    shells[1].alpha[0] =   5.0331513;  shells[1].coef[0] = -0.09996723;
    shells[1].alpha[1] =   1.1695961;  shells[1].coef[1] =  0.39951283; 
    shells[1].alpha[2] =   0.3803890;  shells[1].coef[2] =  0.70011547;

    shells[2].am = 1;  shells[2].nprim = 3;
    shells[2].x = 0.00000; shells[2].y = -0.14322; shells[2].z = 0.0;
    shells[2].alpha[0] =   5.0331513;  shells[2].coef[0] =  0.15591627;
    shells[2].alpha[1] =   1.1695961;  shells[2].coef[1] =  0.60768372; 
    shells[2].alpha[2] =   0.3803890;  shells[2].coef[2] =  0.39195739;

    // hydrogen 1
    shells[3].am = 0;  shells[3].nprim = 3;
    shells[3].x =  1.63804; shells[3].y = 1.13654; shells[3].z = 0.0;
    shells[3].alpha[0] =  3.42525091;  shells[3].coef[0] =  0.15432897;
    shells[3].alpha[1] =  0.62391373;  shells[3].coef[1] =  0.53532814; 
    shells[3].alpha[2] =  0.16885540;  shells[3].coef[2] =  0.44463454;

    // hydrogen 2
    shells[4].am = 0;  shells[4].nprim = 3;
    shells[4].x = -1.63804; shells[4].y = 1.13654; shells[4].z = 0.0;
    shells[4].alpha[0] =  3.42525091;  shells[4].coef[0] =  0.15432897;
    shells[4].alpha[1] =  0.62391373;  shells[4].coef[1] =  0.53532814; 
    shells[4].alpha[2] =  0.16885540;  shells[4].coef[2] =  0.44463454;


    ////////////////////////////////////////
    // 3. Normalize
    // Arguments to simint_normalize_shells:
    //   number of shells
    //   pointer to array of shells
    simint_normalize_shells(5, shells);

    // we need the starting indices of the individual basis functions within the
    // shells, as well as how many basis functions are in each shell.

    //                     s  s  p  s  s
    int shell_start[5] = { 0, 1, 2, 5, 6 };
    int shell_nbf[5]   = { 1, 1, 3, 1, 1 };


    ////////////////////////////////////////
    // 4. Allocate some memory
    // Final (packed) storage
    // This example has 7 basis functions. That would
    // be 7^4 = 2401 unpacked integrals, or 406 packed integrals
    double * integrals = SIMINT_ALLOC(406 * sizeof(double));

    // This is where individual calls to simint_compute_eri will place
    // its integrals.
    // This must be big enough to hold an integral computation for the highest
    // angular momentum in the basis set (multiplied by number of batches, however
    // we aren't doing batches here).
    // For water/sto-3g, highest am = p, so we have to hold
    // a (p p | p p) cartesian integral, which is 3*3*3*3 = 81 elements
    double * target = SIMINT_ALLOC(81 * sizeof(double));


    // Shared workspace. This prevents needless allocation/deallocation
    // for each call
    //
    // No need to multiply by sizeof(double) here. It is included
    // in SIMINT_OSTEI_MAX_WORKMEM
    //
    // (for your information, OSTEI = Obara-Saika Two-Electron Integral)
    double * work = SIMINT_ALLOC(SIMINT_OSTEI_MAX_WORKMEM);


    ////////////////////////////////////////
    // 5. Loop over the shells
    // The following is somewhat naive due
    // to only calculating the unique shell
    // quartets, but is just a demonstration
    ////////////////////////////////////////

    // number of shell quartets we calculated
    int ncomputed_shell = 0;

    // number of integrals calcualated
    int ncomputed_integrals = 0;

    // we only have to initialize these once
    struct simint_multi_shellpair left_pair;
    struct simint_multi_shellpair right_pair;
    simint_initialize_multi_shellpair(&left_pair);
    simint_initialize_multi_shellpair(&right_pair);

    for(int i = 0; i < 5; i++)
    for(int j = 0; j <= i; j++)
    {
        int ij = (i*(i+1))/2 + j;

        // form the left shell pair
        simint_create_multi_shellpair(1, &shells[i],
                                      1, &shells[j], &left_pair, 0);

        for(int k = 0; k < 5; k++)
        for(int l = 0; l <= k; l++)
        {
            int kl = (k*(k+1))/2 + l;
            if(ij < kl) // skip due to permutational symmetry?
                continue;

            simint_create_multi_shellpair(1, &shells[k],
                                          1, &shells[l], &right_pair, 0);


            /////////////////////////////////////
            // Integrals actually computed here
            /////////////////////////////////////
            
            // Should always return one in this case (since we aren't batching)
            simint_compute_eri(&left_pair, &right_pair, 0.0, work, target);

            // keep track of how many shells and integrals we calculated
            ncomputed_shell++;
            ncomputed_integrals += shell_nbf[i] * shell_nbf[j] * shell_nbf[k] * shell_nbf[l];


            // index of where we are in the intermediate "target" buffer
            int target_idx = 0;

            // Now, determine where the integrals go in the final integral storage and put them there
            for(int m = 0; m < shell_nbf[i]; m++)
            for(int n = 0; n < shell_nbf[j]; n++)
            for(int o = 0; o < shell_nbf[k]; o++)
            for(int p = 0; p < shell_nbf[l]; p++)
            {
                int m_idx = shell_start[i] + m;
                int n_idx = shell_start[j] + n;
                int o_idx = shell_start[k] + o;
                int p_idx = shell_start[l] + p;

                int mn_idx = m_idx < n_idx ? (n_idx*(n_idx+1))/2 + m_idx : (m_idx*(m_idx+1))/2 + n_idx;
                int op_idx = o_idx < p_idx ? (p_idx*(p_idx+1))/2 + o_idx : (o_idx*(o_idx+1))/2 + p_idx;
                int mnop_idx = mn_idx < op_idx ? (op_idx*(op_idx+1))/2 + mn_idx : (mn_idx*(mn_idx+1))/2 + op_idx;

                integrals[mnop_idx] = target[target_idx];

                target_idx++;
            }
        }
    }

    //////////////////////////
    // Print out some info
    //////////////////////////
    printf("\n");
    printf("    Number of shell quartets (unpacked): %d\n", 5*5*5*5);
    printf("  Number of shell quartets (calculated): %d\n", ncomputed_shell);
    printf("         Number of integrals (unpacked): %d\n", 7*7*7*7);
    printf("       Number of integrals (calculated): %d\n", ncomputed_integrals);
    printf("\n");


    // cleanup
    simint_free_multi_shellpair(&left_pair);
    simint_free_multi_shellpair(&right_pair);

    for(int i = 0; i < 5; i++)
        simint_free_shell(&shells[i]);

    SIMINT_FREE(target);
    SIMINT_FREE(work);
    SIMINT_FREE(integrals);

    simint_finalize();

    return 0;
}
