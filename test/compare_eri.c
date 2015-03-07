#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "valeev/valeev.h"
#include "eri/eri.h"

int main(int argc, char ** argv)
{

    // These are for valeev comparison
    double A[3];
    double B[3];
    double C[3];
    double D[3];

    int ntest = atoi(argv[1]); 

    double * Ax = (double *)malloc(ntest * sizeof(double));
    double * Ay = (double *)malloc(ntest * sizeof(double));
    double * Az = (double *)malloc(ntest * sizeof(double));
    double * Aa = (double *)malloc(ntest * sizeof(double));

    double * Bx = (double *)malloc(ntest * sizeof(double));
    double * By = (double *)malloc(ntest * sizeof(double));
    double * Bz = (double *)malloc(ntest * sizeof(double));
    double * Ba = (double *)malloc(ntest * sizeof(double));

    double * Cx = (double *)malloc(ntest * sizeof(double));
    double * Cy = (double *)malloc(ntest * sizeof(double));
    double * Cz = (double *)malloc(ntest * sizeof(double));
    double * Ca = (double *)malloc(ntest * sizeof(double));

    double * Dx = (double *)malloc(ntest * sizeof(double));
    double * Dy = (double *)malloc(ntest * sizeof(double));
    double * Dz = (double *)malloc(ntest * sizeof(double));
    double * Da = (double *)malloc(ntest * sizeof(double));

    double * res = (double *)malloc(ntest * sizeof(double));
    double * vres = (double *)malloc(ntest * sizeof(double));

    srand(time(NULL));
    double rmax = RAND_MAX;

    for(int i = 0; i < ntest; i++)
    {
        Ax[i] = 1.0 * rand() / rmax - 0.5;
        Ay[i] = 1.0 * rand() / rmax - 0.5;
        Az[i] = 1.0 * rand() / rmax - 0.5;
        Aa[i] = 100.0 * rand() / rmax;

        Bx[i] = 1.0 * rand() / rmax - 0.5;
        By[i] = 1.0 * rand() / rmax - 0.5;
        Bz[i] = 1.0 * rand() / rmax - 0.5;
        Ba[i] = 100.0 * rand() / rmax;

        Cx[i] = 1.0 * rand() / rmax - 0.5;
        Cy[i] = 1.0 * rand() / rmax - 0.5;
        Cz[i] = 1.0 * rand() / rmax - 0.5;
        Ca[i] = 100.0 * rand() / rmax;

        Dx[i] = 1.0 * rand() / rmax - 0.5;
        Dy[i] = 1.0 * rand() / rmax - 0.5;
        Dz[i] = 1.0 * rand() / rmax - 0.5;
        Da[i] = 10.0 * rand() / rmax;
    }     


    eri_ssss(ntest,
             Ax, Ay, Az,
             Bx, By, Bz,
             Cx, Cy, Cz,
             Dx, Dy, Dz,
             Aa, Ba, Ca, Da,
             res);



    Valeev_Init();

    for(int i = 0; i < ntest; i++)
    {
        A[0] = Ax[i]; A[1] = Ay[i]; A[2] = Az[i];
        B[0] = Bx[i]; B[1] = By[i]; B[2] = Bz[i];
        C[0] = Cx[i]; C[1] = Cy[i]; C[2] = Cz[i];
        D[0] = Dx[i]; D[1] = Dy[i]; D[2] = Dz[i];

        vres[i] = Valeev_eri(0, 0, 0, Aa[i], A,
                             0, 0, 0, Ba[i], B, 
                             0, 0, 0, Ca[i], C, 
                             0, 0, 0, Da[i], D, 1);
    }

    Valeev_Finalize();

    for(int i = 0; i < ntest; i++)
    {
        double diff = fabs(res[i] - vres[i]);
        if(diff > 1e-14)
          printf("%20.15e  %20.15e  %20.15e\n", res[i], vres[i], res[i]-vres[i]);
    }


    free(Ax); free(Ay); free(Az); free(Aa);
    free(Bx); free(By); free(Bz); free(Ba);
    free(Cx); free(Cy); free(Cz); free(Ca);
    free(Dx); free(Dy); free(Dz); free(Da);
    free(res); free(vres);

    return 0;
}
