#include <math.h> // for fabs()

#include "simint/simint_eri.h"
#include "simint/ostei/ostei_config.h"
#include "simint/constants.h"
#include "simint/shell/shell_screen.h"
#include "simint/vectorization/vectorization.h"

double
simint_shellscreen(struct simint_shell const * A,
                   struct simint_shell const * B,
                   int screen_method)
{
    switch(screen_method)
    {
        case SIMINT_SCREEN_SCHWARZ:
            return simint_shellscreen_schwarz(A, B);
        case SIMINT_SCREEN_FASTSCHWARZ:
            return simint_shellscreen_fastschwarz(A, B);
        default:
            return simint_shellscreen_schwarz(A, B);
    }
}


double
simint_primscreen(struct simint_shell const * A,
                  struct simint_shell const * B,
                  double * out,
                  int screen_method)
{
    switch(screen_method)
    {
        case SIMINT_SCREEN_SCHWARZ:
            return simint_primscreen_schwarz(A, B, out);
        case SIMINT_SCREEN_FASTSCHWARZ:
            return simint_primscreen_fastschwarz(A, B, out);
        default:
            return simint_primscreen_schwarz(A, B, out);
    }
}


///////////////////////////////////////
// Shell screening implementations
///////////////////////////////////////
double
simint_shellscreen_schwarz(struct simint_shell const * A,
                           struct simint_shell const * B)
{
    // workspace
    double * work = (double *)SIMINT_ALLOC(SIMINT_OSTEI_MAX_WORKMEM);

    const int ncart1 = ((A->am+1) * (A->am+2))/2;
    const int ncart2 = ((B->am+1) * (B->am+2))/2;
    const int ncart12 = ncart1*ncart2;
    const int ncart1234 = ncart12*ncart12;

    // holds the calculated integrals
    double integrals[ncart1234] SIMINT_ALIGN_ARRAY_DBL;

    struct simint_multi_shellpair P;
    simint_initialize_multi_shellpair(&P);
    simint_create_multi_shellpair(1, A, 1, B, &P, 0);

    // calculate (ab|ab)
    simint_compute_eri(&P, &P, 0.0, work, integrals);

    // find the max value
    double max = 0;
    for(int n = 0; n < ncart1234; n++)
    {
        double abint = fabs(integrals[n]);
        max = ( abint > max ? abint : max );
    }

    simint_free_multi_shellpair(&P);
    SIMINT_FREE(work);
    return max;
}


double
simint_shellscreen_fastschwarz(struct simint_shell const * A,
                               struct simint_shell const * B)
{
    const int same_shell = compare_shell(A, B);

    int idx = 0;
    double val = 0.0;

    for(int i = 0; i < A->nprim; i++)
    {
        const int Bend = (same_shell ? (i+1) : B->nprim);
        for(int j = 0; j < Bend; j++)
        {
            const double a = A->alpha[i];
            const double b = B->alpha[j];
            const double p = a + b;
            const double oop = 1.0/p;
            const double rho = (a*b)/p;
            const double Rx = (A->x - B->x);
            const double Ry = (A->y - B->y);
            const double Rz = (A->z - B->z);
            const double R2 = (Rx*Rx + Ry*Ry + Rz*Rz);

            // this is actually Gab**2 (where Gab is from from MEST)
            double Gab = SQRT_TWO_TIMES_PI_52 * pow(oop, 2.5) * exp(-2.0 * rho * R2);
            Gab *= A->coef[i] * A->coef[i] * B->coef[j] * B->coef[j];
            val += Gab;
        }
    }

    return val;
}




///////////////////////////////////////////
// Primitive screening implementations
///////////////////////////////////////////
double
simint_primscreen_schwarz(struct simint_shell const * A,
                          struct simint_shell const * B,
                          double * out)
{
    const int ncart1 = ((A->am+1) * (A->am+2))/2;
    const int ncart2 = ((B->am+1) * (B->am+2))/2;
    const int ncart12 = ncart1*ncart2;
    const int ncart1234 = ncart12*ncart12;

    // holds the calculated integrals
    double integrals[ncart1234] SIMINT_ALIGN_ARRAY_DBL;

    // here, we basically uncontract the shells
    struct simint_shell new_A, new_B;
    simint_initialize_shell(&new_A);
    simint_initialize_shell(&new_B);

    const int same_shell = compare_shell(A, B);

    struct simint_multi_shellpair P;
    simint_initialize_multi_shellpair(&P);

    // workspace
    double * work = (double *)SIMINT_ALLOC(SIMINT_OSTEI_MAX_WORKMEM);

    double total_max = 0.0;

    int idx = 0;
    int ij = 0;
    for(int i = 0; i < A->nprim; i++)
    {
        simint_create_shell(1, A->am, A->x, A->y, A->z,
                            A->alpha + i, A->coef + i,
                            &new_A);

        const int Bend = (same_shell ? (i+1) : B->nprim);
        for(int j = 0; j < Bend; j++)
        {
            // I don't think we need a factor of 2 for the off-diagonal,
            // do we?
            simint_create_shell(1, B->am, B->x, B->y, B->z,
                                B->alpha + j, B->coef + j,
                                &new_B);

            // create a shell pair
            // Note - we aren't screening this. That would be a infinite loop.
            simint_create_multi_shellpair(1, &new_A, 1, &new_B, &P, 0);

            // calculate (ab|ab)
            simint_compute_eri(&P, &P, 0.0, work, integrals);

            // find the maximum value and store in the output array
            double max = 0;
            for(int m = 0; m < ncart1; m++)
            {
                const int idx_m = m*ncart2*ncart1*ncart2 +
                                  m*ncart2;

                for(int n = 0; n < ncart2; n++)
                {
                    const int idx_mn = idx_m +
                                    n*ncart1*ncart2 +
                                    n;

                    double abint = fabs(integrals[idx_mn]);
                    max = ( abint > max ? abint : max );
                }
            }

            // we avoid doing the sqrt by squaring the tolerance later
            if(out != NULL)
                out[idx++] = max;

            if(max > total_max)
                total_max = max;

            ij++;

            // pad out, if needed
            if( (ij % SIMINT_NSHELL_SIMD) == 0 || ij >= A->nprim*B->nprim)
            {
                while(idx < SIMINT_SIMD_ROUND(idx))
                    out[idx++] = 0.0;
            }
        }
    }

    simint_free_multi_shellpair(&P);
    simint_free_shell(&new_A);
    simint_free_shell(&new_B);
    SIMINT_FREE(work);

    return total_max;
}


double
simint_primscreen_fastschwarz(struct simint_shell const * A,
                              struct simint_shell const * B,
                              double * restrict out)
{
    // workspace
    double * work = (double *)SIMINT_ALLOC(SIMINT_OSTEI_MAX_WORKMEM);

    const int same_shell = compare_shell(A, B);
    double total_max = 0.0;
    int idx = 0;

    // we manually calculate [00|00] for all primitives
    // (manually, rather than calling to simint_compute_eri)
    // Since this is schwarz screening, we calculate (ij|ij)
    // for each primitive pair
    for(int i = 0; i < A->nprim; i++)
    {
        const int Bend = (same_shell ? (i+1) : B->nprim);
        for(int j = 0; j < Bend; j++)
        {
            const double a = A->alpha[i];
            const double b = B->alpha[j];
            const double p = a + b;
            const double oop = 1.0/p;
            const double rho = (a*b)/p;
            const double Rx = (A->x - B->x);
            const double Ry = (A->y - B->y);
            const double Rz = (A->z - B->z);
            const double R2 = (Rx*Rx + Ry*Ry + Rz*Rz);

            // this is actually Gab**2 (where Gab is from from MEST)
            double Gab = SQRT_TWO_TIMES_PI_52 * pow(oop, 2.5) * exp(-2.0 * rho * R2);
            Gab *= A->coef[i] * A->coef[i] * B->coef[j] * B->coef[j];

            if(out != NULL)
                out[idx] = Gab;

            if(Gab > total_max)
                total_max = Gab;

            idx++; // increment even if we aren't outputting
                   // this helps the compiler vectorize this loop
        }
    }

    SIMINT_FREE(work);
    return total_max;
}
