#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "eri/eri.h"
#include "boys/boys.h"

int main(int argc, char ** argv)
{
    int ntest1 = atoi(argv[1]);
    int ntest2 = atoi(argv[2]);
    int ntest3 = atoi(argv[3]);
    int ntest4 = atoi(argv[4]);

    int ntest12 = ntest1*ntest2;
    int ntest34 = ntest3*ntest4;
    int ntest1234 = ntest12*ntest34;

    // for the shells
    double Ax, Ay, Az;
    double Bx, By, Bz;
    double Cx, Cy, Cz;
    double Dx, Dy, Dz;

    double * Aa = (double *)malloc(ntest1 * sizeof(double));
    double * Ac = (double *)malloc(ntest1 * sizeof(double));

    double * Ba = (double *)malloc(ntest2 * sizeof(double));
    double * Bc = (double *)malloc(ntest2 * sizeof(double));

    double * Ca = (double *)malloc(ntest3 * sizeof(double));
    double * Cc = (double *)malloc(ntest3 * sizeof(double));

    double * Da = (double *)malloc(ntest4 * sizeof(double));
    double * Dc = (double *)malloc(ntest4 * sizeof(double));


    // for shell pairs
    double * Px = (double *)malloc(ntest12 * sizeof(double));
    double * Py = (double *)malloc(ntest12 * sizeof(double));
    double * Pz = (double *)malloc(ntest12 * sizeof(double));
    double * Pa = (double *)malloc(ntest12 * sizeof(double));
    double * Pc = (double *)malloc(ntest12 * sizeof(double));
    int * Pn1 = (int *)malloc(1 * sizeof(int));
    int * Pn2 = (int *)malloc(1 * sizeof(int));

    double * Qx = (double *)malloc(ntest34 * sizeof(double));
    double * Qy = (double *)malloc(ntest34 * sizeof(double));
    double * Qz = (double *)malloc(ntest34 * sizeof(double));
    double * Qa = (double *)malloc(ntest34 * sizeof(double));
    double * Qc = (double *)malloc(ntest34 * sizeof(double));
    int * Qn1 = (int *)malloc(1 * sizeof(int));
    int * Qn2 = (int *)malloc(1 * sizeof(int));


    double * res1_s = (double *)malloc(1 * sizeof(double));
    double * res1_m = (double *)malloc(1 * sizeof(double));
    double * res2_s = (double *)malloc(1 * sizeof(double));
    double * res2_m = (double *)malloc(1 * sizeof(double));
    double * intwork = (double *)malloc(ntest1234 * sizeof(double));

    srand(time(NULL));
    double rmax = RAND_MAX;

    Ax = 1.0 * rand() / rmax - 0.5;
    Ay = 1.0 * rand() / rmax - 0.5;
    Az = 1.0 * rand() / rmax - 0.5;
    Bx = 1.0 * rand() / rmax - 0.5;
    By = 1.0 * rand() / rmax - 0.5;
    Bz = 1.0 * rand() / rmax - 0.5;
    Cx = 1.0 * rand() / rmax - 0.5;
    Cy = 1.0 * rand() / rmax - 0.5;
    Cz = 1.0 * rand() / rmax - 0.5;
    Dx = 1.0 * rand() / rmax - 0.5;
    Dy = 1.0 * rand() / rmax - 0.5;
    Dz = 1.0 * rand() / rmax - 0.5;

    for(int i = 0; i < ntest1; i++)
    {
        Aa[i] = 100.0 * rand() / rmax;
        Ac[i] = 10.0 * rand() / rmax;
    }

    for(int i = 0; i < ntest2; i++)
    {
        Ba[i] = 100.0 * rand() / rmax;
        Bc[i] = 10.0 * rand() / rmax;
    }

    for(int i = 0; i < ntest3; i++)
    {
        Ca[i] = 100.0 * rand() / rmax;
        Cc[i] = 10.0 * rand() / rmax;
    }

    for(int i = 0; i < ntest4; i++)
    {
        Da[i] = 100.0 * rand() / rmax;
        Dc[i] = 10.0 * rand() / rmax;
    }

    Boys_Init();

    // set up gaussian hell structures
    struct gaussian_shell A, B, C, D;
    A.nprim = ntest1;  B.nprim = ntest2;  C.nprim = ntest3;  D.nprim = ntest4;
    A.x = Ax;  A.y = Ay;  A.z = Az;  A.alpha = Aa;  A.coef = Ac;  A.am = 0;
    B.x = Bx;  B.y = By;  B.z = Bz;  B.alpha = Ba;  B.coef = Bc;  B.am = 0;
    C.x = Cx;  C.y = Cy;  C.z = Cz;  C.alpha = Ca;  C.coef = Cc;  C.am = 0;
    D.x = Dx;  D.y = Dy;  D.z = Dz;  D.alpha = Da;  D.coef = Dc;  D.am = 0;

    // set up the shell pairs
    struct shell_pair P, Q;
    P.nprim1 = Pn1;  P.nprim2 = Pn2;
    P.x = Px;  P.y = Py;  P.z = Pz;  P.alpha = Pa;  P.prefac = Pc;
    Q.nprim1 = Qn1;  Q.nprim2 = Qn2;
    Q.x = Qx;  Q.y = Qy;  Q.z = Qz;  Q.alpha = Qa;  Q.prefac = Qc;

    fill_ss_shell_pair(1, &A, 1, &B, &P);
    fill_ss_shell_pair(1, &C, 1, &D, &Q);

    // Actually calculate
    eri_1pair_ssss_multi(1, &A, 1, &B, Q, res1_m, intwork);
    eri_1pair_ssss_single(A, B, Q, res1_s);
    eri_2pair_ssss_multi(P, Q, res2_m, intwork);
    eri_2pair_ssss_single(P, Q, res2_s);


    Boys_Finalize();

    for(int i = 0; i < 1; i++)
    {
        double diff1 = fabs(res1_m[i] - res1_s[i]);
        double diff2 = fabs(res2_m[i] - res2_s[i]);
        if(diff1 > 1e-10 || diff2 > 1e-10)
          printf("%8.4e  %8.4e  %8.4e  %8.4e  --  %8.4e  %8.4e\n", 
                                      res1_s[i], res1_m[i], res2_s[i], res2_m[i],
                                      diff1, diff2);
    }


    free(Aa); free(Ac);
    free(Ba); free(Bc);
    free(Ca); free(Cc);
    free(Da); free(Dc);

    free(Px); free(Py); free(Pz); free(Pa); free(Pc); free(Pn1); free(Pn2);
    free(Qx); free(Qy); free(Qz); free(Qa); free(Qc); free(Qn1); free(Qn2);

    free(res1_s); free(res1_m); free(res2_s); free(res2_m); free(intwork);

    return 0;
}
