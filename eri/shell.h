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

struct multishell_pair
{
    int am1, am2;          // angular momentum.
    int nprim;             // Total number of primitives
    int nprim_length;      // Actual length of alpha, etc, arrays (!= nprim due to alignment)
    int nshell1, nshell2;  // number of shells
    int * nprim1;          // number of primitives in shells (of length nshell1)
    int * nprim2;          // number of primitives in shells (of length nshell2)

    int nshell12;
    int * nprim12;   // length nshell12;
    int * primstart; // length nshell12
                     // primstart[n] = start of shell pair n
    int * primend;   // length nshell12
                     // primend[n] = end (not inclusive)

    // these are all of length nprim
    double * x;
    double * y;
    double * z;
    double * PA_x;
    double * PA_y;
    double * PA_z;
    double * alpha;
    double * prefac;
};


struct shell_pair
{
    int am1, am2;          // angular momentum of each shell
    int nprim;             // Total number of primitives
    int nprim_length;      // Actual length of alpha, etc, arrays (!= nprim due to alignment)
    int nprim1, nprim2;    // number of primitives in each shell

    // these are all of length nprim
    double * x;
    double * y;
    double * z;
    double * PA_x;
    double * PA_y;
    double * PA_z;
    double * alpha;
    double * prefac;
};


void allocate_gaussian_shell(int nprim, struct gaussian_shell * const restrict G);
void free_gaussian_shell(struct gaussian_shell G);
void normalize_gaussian_shells(int n, struct gaussian_shell * const restrict G);


void allocate_shell_pair_from_shells(struct gaussian_shell const A,
                                     struct gaussian_shell const B,
                                     struct shell_pair * const restrict P);

void allocate_multishell_pair_from_shells(int na, struct gaussian_shell const * const restrict A,
                                          int nb, struct gaussian_shell const * const restrict B,
                                          struct multishell_pair * const restrict P);

void free_shell_pair(struct shell_pair P);
void free_multishell_pair(struct multishell_pair P);



void fill_ss_shell_pair(struct gaussian_shell const A,
                        struct gaussian_shell const B,
                        struct shell_pair * const restrict P);

void fill_ss_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                             int nb, struct gaussian_shell const * const restrict B,
                             struct multishell_pair * const restrict P);

struct shell_pair
create_ss_shell_pair(struct gaussian_shell const A,
                     struct gaussian_shell const B);

struct multishell_pair
create_ss_multishell_pair(int na, struct gaussian_shell const * const restrict A,
                          int nb, struct gaussian_shell const * const restrict B);

#endif

