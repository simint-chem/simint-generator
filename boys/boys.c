// in Boys_taylorgrid.c
void Boys_taylorgrid_Init(double max_x, int max_n);
void Boys_taylorgrid_Finalize(void);



void Boys_Init(double max_x, int max_n)
{
    Boys_taylorgrid_Init(max_x, max_n);
}

void Boys_Finalize(void)
{
    Boys_taylorgrid_Finalize();
}


// Maximum value needed for Boys function
double Boys_Max(const int ncenter,
                const double * X, const double * Y, const double * Z,
                const int * n_prim_per_center, const double * alpha)
{
    int i, j, a, b;

    /*
    int maxpa = 0;
    int maxpb = 0;
    */

    double maxd = 0;
    double d, d2;

    int pa = 0; // counts over alpha
    int pb = 0; // counts over alpha

    // these hold the start of this center in the alpha array
    int paorig = 0;
    int pborig = 0;

    for(i = 0; i < ncenter; i++)
    {
        pborig = 0;
        for(j = 0; j < ncenter; j++)
        {
            d = (X[i]-X[j])*(X[i]-X[j])
              + (Y[i]-Y[j])*(Y[i]-Y[j])
              + (Z[i]-Z[j])*(Z[i]-Z[j]);

            pa = paorig;
            for(a = 0; a < n_prim_per_center[i]; a++)
            {
                pb = pborig;
                for(b = 0; b < n_prim_per_center[j]; b++)
                {

                    d2 = d * (alpha[pa]*alpha[pb]) / (alpha[pa] + alpha[pb]);
                    if(d2 > maxd)
                    {
                        /*
                        maxpa = pa;
                        maxpb = pb;
                        */
                        maxd = d2;
                    }
                    pb++;
                }
                pa++;
            }
            pborig += n_prim_per_center[j];
        }
        paorig += n_prim_per_center[i];
    }

    //printf("Max i, j: %d %d\n", maxpa, maxpb);

    // missing a factor of 2
    return 2.0 * maxd;;
}
