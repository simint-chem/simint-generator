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

    // these are all of length nprim
    double * x;
    double * y;
    double * z;
    double * alpha;
    double * prefac;
};



#endif

