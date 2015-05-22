#ifndef SIMINT_SHELL_H
#define SIMINT_SHELL_H

#ifdef __cplusplus
extern "C" {
#endif

struct gaussian_shell
{
    int am;
    double x, y, z;
    int nprim;
    double * alpha;
    double * coef;
};

struct multishell_pair
{
    int am1, am2;          // angular momentum.
    int nprim;             // Total number of primitives
    int nprim_length;      // Actual length of alpha, etc, arrays (!= nprim due to alignment)
    int nshell1, nshell2;  // number of shells

    int nshell12;
    int * nprim12;   // length nshell12;
    int * primstart; // length nshell12
                     // primstart[n] = start of shell pair n
    int * primend;   // length nshell12
                     // primend[n] = end (not inclusive)

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
};


struct multishell_pair_flat
{
    int am1, am2;          // angular momentum.
    int nprim;             // Total number of primitives
    int nshell1, nshell2;  // number of shells
    int nshell12;

    // length nshell12
    double * AB_x;
    double * AB_y;
    double * AB_z;

    // these are all of length nprim
    int * shellidx;
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
};



void allocate_gaussian_shell(int nprim, struct gaussian_shell * const restrict G);
void free_gaussian_shell(struct gaussian_shell G);
void normalize_gaussian_shells(int n, struct gaussian_shell * const restrict G);


void free_multishell_pair(struct multishell_pair P);

void free_multishell_pair_flat(struct multishell_pair_flat P);


void fill_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                          int nb, struct gaussian_shell const * const restrict B,
                          struct multishell_pair * const restrict P);

void fill_multishell_pair_flat(int na, struct gaussian_shell const * const restrict A,
                               int nb, struct gaussian_shell const * const restrict B,
                               struct multishell_pair_flat * const restrict P);

struct multishell_pair
create_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                       int nb, struct gaussian_shell const * const restrict B);

struct multishell_pair_flat
create_multishell_pair_flat(int na, struct gaussian_shell const * const restrict A,
                            int nb, struct gaussian_shell const * const restrict B);

#ifdef __cplusplus
}
#endif

#endif

