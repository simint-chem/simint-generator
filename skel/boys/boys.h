#ifndef BOYS_H
#define BOYS_H

#ifdef __cplusplus
extern "C" {
#endif

// Prototypes
void Boys_Init(void);
void Boys_Finalize(void);

double Boys_Max(const int ncenter,
                const double * X, const double * Y, const double * Z,
                const int * n_prim_per_center, const double * alpha);

#ifdef __cplusplus
}
#endif

#endif
