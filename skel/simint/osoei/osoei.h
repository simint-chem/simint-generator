#pragma once

#include "simint/shell/shell.h"

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif


/*! \brief Terms required by some one electron integrals */
void simint_osoei_overlap_terms(const double alpha1, const double * xyz1,
                                const double alpha2, const double * xyz2,
                                int nam1, int nam2,
                                double * restrict terms);


void simint_osoei_ke_terms(const double alpha1, const double * xyz1,
                           const double alpha2, const double * xyz2,
                           int nam1, int nam2,
                           double * restrict terms);


int simint_compute_osoei_overlap(struct simint_shell const * sh1,
                                 struct simint_shell const * sh2,
                                 double * restrict integrals);

int simint_compute_osoei_ke(struct simint_shell const * sh1,
                            struct simint_shell const * sh2,
                            double * restrict integrals);

int simint_compute_osoei_potential(int ncenter,
                                   double * Z, double * x, double * y, double * z,
                                   struct simint_shell const * sh1,
                                   struct simint_shell const * sh2,
                                   double * restrict integrals);

#ifdef __cplusplus
}
#endif
