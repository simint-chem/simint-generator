#ifndef SIMINT_SHELL_H
#define SIMINT_SHELL_H

struct gaussian_shell
{
    int am;
    double x, y, z;
    int nprim;
    double * alpha;
    double * coef;
};

struct shell_pair
{
    int am1, am2;  // angular momentum.
    int nprim;  // Total number of primitives
    int nshell1, nshell2;  // number of shells
    int * nprim1;  // number of primitives in shells (of length nshell1)
    int * nprim2;  // number of primitives in shells (of length nshell2)

    int nshell12;
    int * nprim12; // length nshell12;

    // these are all of length nprim
    double * x;
    double * y;
    double * z;
    double * alpha;
    double * prefac;
};

struct shell_pair 
allocate_shell_pair_from_shells(int na, struct gaussian_shell const * const restrict A,
                                int nb, struct gaussian_shell const * const restrict B);

void free_shell_pair(struct shell_pair P);



void fill_ss_shell_pair(int na, struct gaussian_shell const * const restrict A,
                        int nb, struct gaussian_shell const * const restrict B,
                        struct shell_pair * const restrict P);

struct shell_pair
create_ss_shell_pair(int na, struct gaussian_shell const * const restrict A,
                     int nb, struct gaussian_shell const * const restrict B);

#endif

