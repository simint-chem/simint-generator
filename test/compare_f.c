#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "valeev/valeev.h"
#include "boys/boys.h"

#define MAXN 5

int main(int argc, char ** argv)
{

    double maxx = atof(argv[1]); 

    Valeev_Init();
    Boys_Init(maxx, MAXN+7);

    srand(time(NULL));
    //double rmax = RAND_MAX;

    double maxerror = 0.0;

    double * my = (double *)malloc( (MAXN+1) * sizeof(double));
    double * valeev = (double *)malloc( (MAXN+1) * sizeof(double));
    

    for(double x = 0.0; x < maxx; x += 0.01)
    {
        Valeev_F(valeev, MAXN, x);
        Boys_F(my, MAXN, x);

        for(int n = 0; n <= MAXN; n++)
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
