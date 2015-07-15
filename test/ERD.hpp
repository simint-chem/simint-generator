#ifndef SIMINT_ERD_H
#define SIMINT_ERD_H

#include "shell/shell.h"
#include "test/common.hpp"
#include "test/timer.hpp"

typedef std::vector<gaussian_shell, AlignedAllocator<gaussian_shell>> AlignedGaussianVec;

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


        TimerInfo Compute_shell(struct gaussian_shell const A,
                                struct gaussian_shell const B,
                                struct gaussian_shell const C,
                                struct gaussian_shell const D,
                                double * integrals);


        TimerInfo Integrals(const AlignedGaussianVec & g1, const AlignedGaussianVec & g2,
                            const AlignedGaussianVec & g3, const AlignedGaussianVec & g4,
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
};


#endif
