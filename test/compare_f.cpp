#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "constants.h"
#include "valeev/valeev.h"
#include "boys/boys.h"
#include "boys/boys_split.h"
#include "boys/boys_FO.h"

#define MAXN 2

int main(int argc, char ** argv)
{
    int nval = atoi(argv[1]);
    double maxx = atof(argv[2]);

    if(nval > MAXN)
    {
        printf("\nError - N value too high!\n");
        return 1;
    }

    Valeev_Init();
    Boys_Init(maxx, MAXN+7);

    srand(time(NULL));
    //double rmax = RAND_MAX;

    double maxerr_taylor = 0.0;
    double maxerr_FO = 0.0;

    double maxrelerr_taylor = 0.0;
    double maxrelerr_FO = 0.0;

    double valeev[MAXN+1];
    double my_taylor[MAXN+1];
    double my_FO[MAXN+1];

    for(double x = 0.0; x < maxx; x += 0.0001)
    {
        Valeev_F(valeev, nval, x);
        Boys_F_taylor(my_taylor, nval, x);
        Boys_F_FO(my_FO, nval, x);

        double diff_taylor = fabs(my_taylor[nval] - valeev[nval]);
        double diff_FO     = fabs(my_FO[nval]     - valeev[nval]);
        double rel_taylor = diff_taylor / valeev[nval];
        double rel_FO = diff_FO / valeev[nval];

        if(diff_FO > maxerr_FO)
            maxerr_FO = diff_FO;
        if(diff_taylor > maxerr_taylor)
            maxerr_taylor = diff_taylor;

        if(rel_taylor > maxrelerr_taylor)
            maxrelerr_taylor = rel_taylor; 
        if(rel_FO > maxrelerr_FO)
            maxrelerr_FO = rel_FO; 

        //printf("%12.8f     %25.17e  %25.17e  %25.17e\n", x, valeev[nval], my_taylor[nval], my_FO[nval]);
    }

    printf("\n");
    printf("***Max abserr - taylor: %12.8e\n", maxerr_taylor);
    printf("***Max relerr - taylor: %12.7e\n", maxrelerr_taylor);
    printf("\n");
    printf("***Max abserr -     FO: %12.8e\n", maxerr_FO);
    printf("***Max relerr -     FO: %12.8e\n", maxrelerr_FO);
    printf("\n");

    printf("\n");

    Valeev_Finalize();
    Boys_Finalize();

    return 0;
}
