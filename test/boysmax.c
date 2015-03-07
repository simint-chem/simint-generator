/*
 * Tests the Boys_Max function vs. a search via a 4-index loop
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "boys/boys.h"

// Max number of primitives per center
#define MAXPERCEN 20

int main(int argc, char ** argv)
{
    int ncenter = atoi(argv[1]); 
    int ntest = atoi(argv[2]); 
    int i,j,k,l, n;

    int * npercenter = (int *)malloc(ncenter * sizeof(int));
    double * Ax = (double *)malloc(ncenter * sizeof(double));
    double * Ay = (double *)malloc(ncenter * sizeof(double));
    double * Az = (double *)malloc(ncenter * sizeof(double));
    double * Aa = (double *)malloc(MAXPERCEN * ncenter * sizeof(double));
    int * nmap = (int *)malloc(MAXPERCEN * ncenter * sizeof(int));

    srand(time(NULL));
    double rmax = RAND_MAX;

    int failures = 0;

    for(n = 0; n < ntest; n++)
    {
        int pa = 0;
        for(i = 0; i < ncenter; i++)
        {
            Ax[i] = 100.0 * rand() / rmax - 50;
            Ay[i] = 100.0 * rand() / rmax - 50;
            Az[i] = 100.0 * rand() / rmax - 50;

            npercenter[i] = rand() % (MAXPERCEN-1)+1;
            for(j = 0; j < npercenter[i]; j++)
            {
              nmap[pa] = i;
              Aa[pa++] = 1000.0 * rand() / rmax;
            }
        }     
       
        pa = 0; 
        for(i = 0; i < ncenter; i++)
        {
            printf("  %d :  %12.8f %12.8f %12.8f\n", i, Ax[i], Ay[i], Az[i]);
            for(j = 0; j < npercenter[i]; j++)
            {
                printf("  %d     %12.8f\n", pa, Aa[pa]);
                pa++;
            }
        }

        int nprim = pa;

        // from my function
        double fmax = Boys_Max(ncenter, Ax, Ay, Az, npercenter, Aa);  
 
        // brute force parameter to boys
        double maxparam = 0;
        int maxi, maxj, maxk, maxl;

        for(i = 0; i < nprim; i++)
        for(j = 0; j < nprim; j++)
        for(k = 0; k < nprim; k++)
        for(l = 0; l < nprim; l++)
        {
            double a_x = Ax[nmap[i]];
            double a_y = Ay[nmap[i]];
            double a_z = Az[nmap[i]];
            double b_x = Ax[nmap[j]];
            double b_y = Ay[nmap[j]];
            double b_z = Az[nmap[j]];
            double c_x = Ax[nmap[k]];
            double c_y = Ay[nmap[k]];
            double c_z = Az[nmap[k]];
            double d_x = Ax[nmap[l]];
            double d_y = Ay[nmap[l]];
            double d_z = Az[nmap[l]];

            double Pab_x = (a_x*Aa[i] + b_x*Aa[j])/(Aa[i]+Aa[j]);
            double Pab_y = (a_y*Aa[i] + b_y*Aa[j])/(Aa[i]+Aa[j]);
            double Pab_z = (a_z*Aa[i] + b_z*Aa[j])/(Aa[i]+Aa[j]);
            double Pcd_x = (c_x*Aa[k] + d_x*Aa[l])/(Aa[k]+Aa[l]);
            double Pcd_y = (c_y*Aa[k] + d_y*Aa[l])/(Aa[k]+Aa[l]);
            double Pcd_z = (c_z*Aa[k] + d_z*Aa[l])/(Aa[k]+Aa[l]);
    
            // Calculate R
            double Rx = Pab_x - Pcd_x;
            double Ry = Pab_y - Pcd_y;
            double Rz = Pab_z - Pcd_z;
            double R2 = Rx*Rx + Ry*Ry + Rz*Rz;
    
            // Parameter to boys function
            double param = R2 * ((Aa[i]+Aa[j])*(Aa[k]+Aa[l]))/(Aa[i]+Aa[j]+Aa[k]+Aa[l]);
            if(param > maxparam)
            {
                maxi = i;
                maxj = j;
                maxk = k;
                maxl = l;
                maxparam = param;
            }
        }

        printf("Maximum parameter from 4-index search: ( %d %d | %d %d )\n", maxi, maxj, maxk, maxl);
        printf("  From long search: %12.8e\n", maxparam);
        printf("     From Boys_Max: %12.8e\n", fmax);
        printf("              Diff: %12.8e\n", fabs(maxparam-fmax));
    }

    free(Ax); free(Ay); free(Az); free(Aa);
    free(nmap);
    free(npercenter);

    return failures;
}
