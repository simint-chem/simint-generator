#ifndef SIMINT_ERD_H
#define SIMINT_ERD_H

#include "shell/shell.h"
#include "test/common.hpp"
#include "test/timer.h"

typedef std::vector<gaussian_shell> GaussianVec;

class ERD_ERI
{

    public:
        ERD_ERI(int am, int nprim, int ncgto);

        ERD_ERI(int na, struct gaussian_shell const * const restrict A,
                int nb, struct gaussian_shell const * const restrict B,
                int nc, struct gaussian_shell const * const restrict C,
                int nd, struct gaussian_shell const * const restrict D);

        ERD_ERI(int am1, int nprim1, int ncgto1,
                int am2, int nprim2, int ncgto2,
                int am3, int nprim3, int ncgto3,
                int am4, int nprim4, int ncgto4);

        ~ERD_ERI(void); 



        TimerType Integrals(struct gaussian_shell const * const restrict A, int nshell1,
                            struct gaussian_shell const * const restrict B, int nshell2,
                            struct gaussian_shell const * const restrict C, int nshell3,
                            struct gaussian_shell const * const restrict D, int nshell4,
                            double * const integrals);



    private:
        double * dscratch;
        int * iscratch;

        int i_buffer_size;
        int d_buffer_size;

        double * cc;
        double * alpha;

        void Init_(int am1, int nprim1, int ncgto1,
                   int am2, int nprim2, int ncgto2,
                   int am3, int nprim3, int ncgto3,
                   int am4, int nprim4, int ncgto4);

        TimerType Compute_shell_(struct gaussian_shell const A,
                                 struct gaussian_shell const B,
                                 struct gaussian_shell const C,
                                 struct gaussian_shell const D,
                                 double * integrals);
};

void normalize_gaussian_shells_erd(int n, struct gaussian_shell * const restrict G);



// in the ERD library
extern "C" {

void erd__gener_eri_batch_(const int *imax, const int *zmax, const int *nalpha, const int *ncoeff,
                           const int *ncsum, const int *ncgto1, const int *ncgto2,
                           const int *ncgto3, const int *ncgto4, const int *npgto1,
                           const int *npgto2, const int *npgto3, const int *npgto4,
                           const int *shell1, const int *shell2, const int *shell3,
                           const int *shell4, const double *x1, const double *y1, const double *z1,
                           const double *x2,const double *y2,const double *z2, const double *x3,
                           const double *y3,const double *z3,const double *x4, const double *y4, const double *z4,
                           const double *alpha, const double *cc, const int *ccbeg, const int *ccend,
                           const int *spheric,  const int *screen, int *icore,
                           int *nbatch, int * nfirst, double *zcore );

void erd__memory_eri_batch_(const int *nalpha, const int *ncoeff,
                            const int *ncgto1, const int *ncgto2, const int *ncgto3, const int *ncgto4,
                            const int *npgto1, const int *npgto2, const int *npgto3, const int *npgto4,
                            const int *shell1, const int *shell2, const int *shell3, const int *shell4,
                            const double *x1, const double *y1, const double *z1, const double *x2, const double *y2,
                            const double *z2, const double *x3, const double *y3, const double *z3, const double *x4,
                            const double *y4, const double *z4, const double *alpha, const double *cc, const int *spheric,
                            int *imin, int *iopt, int *zmin, int *zopt);
}


#endif
