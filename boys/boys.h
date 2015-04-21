#ifndef BOYS_H
#define BOYS_H

// Prototypes
void Boys_Init(double max_x, int max_n);
void Boys_Finalize(void);

double Boys_Max(const int ncenter,
                const double * X, const double * Y, const double * Z,
                const int * n_prim_per_center, const double * alpha);

#endif
