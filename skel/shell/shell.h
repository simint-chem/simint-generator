#ifndef SIMINT_SHELL_H
#define SIMINT_SHELL_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct gaussian_shell
{
    int am;
    int nprim;

    double x, y, z;

    double * alpha;
    double * coef;
};

struct multishell_pair
{
    int am1, am2;          // angular momentum.
    int nprim;             // Total number of primitives
    int nprim_length;      // Acutal length of alpha, etc, arrays, including padding (for alignment)

    size_t memsize;          // Total memory for doubles (in bytes)

    int nshell1, nshell2;  // number of shells
    int nshell12;          // nshell1 * nshell2
    int * nprim12;    // length nshell12;

    int nbatch;
    int * batchprim;   // length nbatch. number of primitives in a batch, including padding

    // length nshell12
    double * AB_x;
    double * AB_y;
    double * AB_z;

    // these are all of length nprim
    double * x;
    double * y;
    double * z;
    double * PA_x;
    double * PA_y;
    double * PA_z;
    double * bAB_x;
    double * bAB_y;
    double * bAB_z;
    double * alpha;
    double * prefac;
    int * shellidx;
};


void allocate_gaussian_shell(int nprim, struct gaussian_shell * const restrict G);
void free_gaussian_shell(struct gaussian_shell G);
struct gaussian_shell copy_gaussian_shell(const struct gaussian_shell G);
void normalize_gaussian_shells(int n, struct gaussian_shell * const restrict G);


void free_multishell_pair(struct multishell_pair P);

void fill_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                          int nb, struct gaussian_shell const * const restrict B,
                          struct multishell_pair * const restrict P);

struct multishell_pair
create_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                       int nb, struct gaussian_shell const * const restrict B);

#ifdef __cplusplus
}
#endif

#endif

