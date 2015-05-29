#include <cstdio>
#include <fstream>

#include "eri/eri.h"
#include "boys/boys.h"
#include "test/common.hpp"
#include "test/valeev.hpp"

int main(int argc, char ** argv)
{
    // parse command line
    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
        return 1;
    }

    // basis file to read
    std::string basfile(argv[1]);

    // read in the shell info
    ShellMap shellmap = ReadBasis(basfile);

    // normalize
    for(auto & it : shellmap)
        normalize_gaussian_shells(it.second.size(), it.second.data());


    // find the max dimensions
    std::array<int, 3> maxparams = FindMapMaxParams(shellmap);
    const int maxam = maxparams[0];
    //const int maxnprim = maxparams[1];
    const int maxsize = maxparams[2];

    /* Storage of integrals results */
    double * res_valeev = (double *)ALLOC(maxsize * sizeof(double));

    // initialize stuff
    Valeev_Init();

    // open the output file
    std::string reffile = basfile + ".ref";
    std::ofstream out(reffile.c_str(), std::ofstream::binary | std::ofstream::trunc);

    // loop over all quartets, choosing only valid ones
    for(int i = 0; i <= maxam; i++)
    for(int j = 0; j <= maxam; j++)
    for(int k = 0; k <= maxam; k++)
    for(int l = 0; l <= maxam; l++)
    {
        if(!ValidQuartet(i, j, k, l))
            continue;
 
        printf("( %d %d | %d %d ) ... ", i, j, k, l);

        const int ncart1234 = NCART(i) * NCART(j) * NCART(k) * NCART(l);   
        const int nshell1234 = shellmap[i].size() * shellmap[j].size() * shellmap[k].size() * shellmap[l].size();
        const int nsize = ncart1234 * nshell1234;

        // calculate
        ValeevIntegrals(shellmap[i], shellmap[j],
                        shellmap[k], shellmap[l],
                        res_valeev, false);


        // write to output file
        out.write(reinterpret_cast<const char *>(res_valeev), nsize * sizeof(double));

        //for(int m = 0; m < ncart1234; m++)
        //    printf("%20.8e\n", res_valeev[m]);

        printf(" done\n");
    }

    // done with the output file
    out.close();

    FreeShellMap(shellmap);
    Valeev_Finalize();
    FREE(res_valeev);

    return 0;
}
