#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "constants.h"
#include "valeev/valeev.h"
#include "boys/boys.h"

#define MAXN 0

int main(int argc, char ** argv)
{

    double maxx = atof(argv[1]); 

    //const double kfac = F0_KFAC * ONESIX_OVER_SQRT_PI;

    Valeev_Init();
    Boys_Init(maxx, MAXN+7);

    srand(time(NULL));
    //double rmax = RAND_MAX;

    double maxerr_taylor = 0.0;
    double maxerr_cheby = 0.0;
    double maxerr_FO = 0.0;
    double maxerr_FO2 = 0.0;

    double my_taylor = 0.0;
    double my_cheby = 0.0;
    double my_FO = 0.0;
    double my_FO2 = 0.0;
    double valeev = 0.0;
    

    for(double x = 0.0; x < maxx; x += 0.01)
    {
        Valeev_F(&valeev, MAXN, x);

        my_taylor = Boys_F0_taylor(x);
        my_cheby = Boys_F0_cheby(x);
        my_FO = Boys_F0_FO(x);
        my_FO2 = Boys_F0_FO2(x);

        for(int n = 0; n <= MAXN; n++)
        {
            double diff_taylor = fabs(my_taylor - valeev);
            double diff_cheby  = fabs(my_cheby  - valeev);
            double diff_FO     = fabs(my_FO     - valeev);
            double diff_FO2    = fabs(my_FO2    - valeev);

            if(diff_taylor > maxerr_taylor)
                maxerr_taylor = diff_taylor;
            if(diff_cheby > maxerr_cheby)
                maxerr_cheby = diff_cheby;
            if(diff_FO > maxerr_FO)
                maxerr_FO = diff_FO;
            if(diff_FO2 > maxerr_FO2)
                maxerr_FO2 = diff_FO2;
            //printf("%3d %12.8f     %12.8e  %12.8e  %12.8e     %12.8e  %12.8e\n", n, x, my_taylor, my_cheby, valeev, diff_taylor, diff_cheby);
        }
        //printf("----------------------------\n");
    }

    printf("\n");
    printf("***Max error - taylor: %12.8e\n", maxerr_taylor);
    printf("***Max error -  cheby: %12.8e\n", maxerr_cheby);
    printf("***Max error -     FO: %12.8e\n", maxerr_FO);
    printf("***Max error -    FO2: %12.8e\n", maxerr_FO2);
    printf("\n");

    Valeev_Finalize();
    Boys_Finalize();

    return 0;
}
