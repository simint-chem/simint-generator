#ifndef VALEEV_H
#define VALEEV_H

#define MAXFAC 200
#define EPS 1e-17

// initialize math stuff
void Valeev_Init(void);
void Valeev_Finalize(void);

// As original as can be
void Valeev_F(double *F, int n, double x);

double Valeev_eri(int l1, int m1, int n1, double alpha1,
                  const double* A, int l2, int m2, int n2,
                  double alpha2, const double* B, int l3, int m3,
                  int n3, double alpha3, const double* C, int l4,
                  int m4, int n4, double alpha4, const double* D,
                  int norm_flag);

#endif
