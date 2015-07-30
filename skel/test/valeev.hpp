#ifndef VALEEV_H
#define VALEEV_H

#define MAXFAC 200
#define EPS 1e-17

// initialize math stuff
void Valeev_Init(void);
void Valeev_Finalize(void);

// As original as can be
void Valeev_F(long double *F, int n, long double x);

long double Valeev_eri(int l1, int m1, int n1, long double alpha1,
                  const long double* A, int l2, int m2, int n2,
                  long double alpha2, const long double* B, int l3, int m3,
                  int n3, long double alpha3, const long double* C, int l4,
                  int m4, int n4, long double alpha4, const long double* D,
                  int norm_flag);

#endif
