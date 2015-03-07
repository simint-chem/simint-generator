#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "valeev/valeev.h"
#include "boys/boys.h"

#define NN 5

int main(int argc, char ** argv)
{

    double maxx = atof(argv[1]); 

    Valeev_Init();
    Boys_Init();

    srand(time(NULL));
    //double rmax = RAND_MAX;

    double maxerror = 0.0;

    double * my = (double *)malloc( (NN+1) * sizeof(double));
    double * valeev = (double *)malloc( (NN+1) * sizeof(double));
    

    for(double x = 0.0; x < maxx; x += 0.001)
    {
        //x = 40.0 * (rand() / rmax);

        Valeev_F(valeev, NN, x);
        Boys_F(my, NN, x);

        for(int n = 0; n <= NN; n++)
        {
            double diff = fabs(my[n] - valeev[n]);
            if(diff > maxerror)
                maxerror = diff;
            //printf("%3d %12.8f %12.8e  %12.8e     %12.8e\n", n, x, my[n], valeev[n], diff);
        }
        //printf("----------------------------\n");
    }

    printf("\n***Max error: %12.8e\n", maxerror);

    Valeev_Finalize();
    Boys_Finalize();
    free(my);
    free(valeev);

    return 0;
}
