#include <cstdio>
#include <cmath>

#include "boys/boys.h"
#include "test/common.hpp"
#include "test/valeev.hpp"
#include "test/erd_interface.hpp"


int main(int argc, char ** argv)
{
    // set up the function pointers
    Init_Test();

    // parse command line
    if(argc != 2)
    {
        printf("Give me 1 argument! I got %d\n", argc-1);
        return 1;
    }

    // basis functions file to read
    std::string basfile(argv[1]);

    // read in the shell info
    std::map<int, AlignedGaussianVec> shellmap = ReadBasis(basfile);

// todo
// Make a deep copy of the map, and normalize here rather than in the loop
//    // normalize
//    for(const auto & it : shellmap)
//        normalize_gaussian_shells(it.second.size(), it.second.data());

    // find the max dimensions
    int maxnprim = 0;
    int maxam = 0;
    int maxel = 0;
    for(auto & it : shellmap)
    {
        const int nca = NCART(it.first);
        const int nsh = it.second.size();
        const int n = nca * nsh;
        if(n > maxel)
            maxel = n;

        if(it.first > maxam)
            maxam = it.first;

        for(auto & it2 : it.second)
        {
            if(it2.nprim > maxnprim)
                maxnprim = it2.nprim;
        }
    }

    maxel = pow(maxel, 4);

    /* Storage of integrals */
    double * res_FO              = (double *)ALLOC(maxel * sizeof(double));
    double * res_vref            = (double *)ALLOC(maxel * sizeof(double));
    double * res_vref_flat       = (double *)ALLOC(maxel * sizeof(double));
    double * res_liberd          = (double *)ALLOC(maxel * sizeof(double));
    double * res_valeev          = (double *)ALLOC(maxel * sizeof(double));

    // initialize stuff
    Valeev_Init();
    Boys_Init();
    ERD_Init(maxam, maxnprim, 1);
             
             
    printf("\n");
    printf("%17s    %10s  %10s  %10s  %10s    %10s  %10s  %10s  %10s\n", "Quartet", "MaxErr", "", "", "", 
                                                                "MaxRelErr", "", "", "");
    printf("%17s    %10s  %10s  %10s  %10s    %10s  %10s  %10s  %10s\n", "", "FO", "vref", "vrefF", "erd", 
                                                                "FO", "vref", "vrefF", "erd");

    // loop over all quartets, choosing only valid ones
    for(const auto & it_i : shellmap)
    for(const auto & it_j : shellmap)
    for(const auto & it_k : shellmap)
    for(const auto & it_l : shellmap)
    {
        std::array<int, 4> am{it_i.first, it_j.first, it_k.first, it_l.first};

        if(!ValidQuartet(am))
            continue;
 
        const int ncart = NCART(am[0]) * NCART(am[1]) * NCART(am[2]) * NCART(am[3]);   

        // create a vec quartet
        const VecQuartet gshells_orig{it_i.second, it_j.second,
                                      it_k.second, it_l.second};

        // number of shells
        const std::array<int, 4> nshell{gshells_orig[0].size(), gshells_orig[1].size(),
                                        gshells_orig[2].size(), gshells_orig[3].size()};

        const int nshell1234 = nshell[0] * nshell[1] * nshell[2] * nshell[3];


        // make a copy the shells
        VecQuartet gshells(CopyQuartets(gshells_orig));
        VecQuartet gshells_erd(CopyQuartets(gshells_orig));

        
        /////////////////////////////////
        // Normalize the shells
        /////////////////////////////////
        for(auto & it : gshells)
            normalize_gaussian_shells(it.size(), it.data());
        for(auto & it : gshells_erd)
            normalize_gaussian_shells_erd(it.size(), it.data());
       
 
        /////////////////////////////////
        // Calculate or read valeev
        // reference integrals
        /////////////////////////////////
        ValeevIntegrals(gshells, res_valeev, false);
        //ReadValeevIntegrals(basedir, am, res_valeev);


        /////////////////////////////////
        // Calculate with liberd
        /////////////////////////////////
        ERDIntegrals(gshells_erd, res_liberd);


        /////////////////////////////////
        // calculate with my code
        /////////////////////////////////

        // set up shell pairs
        struct multishell_pair P = create_multishell_pair(nshell[0], gshells[0].data(),
                                                          nshell[1], gshells[1].data());
        struct multishell_pair Q = create_multishell_pair(nshell[2], gshells[2].data(),
                                                          nshell[3], gshells[3].data());

        struct multishell_pair_flat Pf = create_multishell_pair_flat(nshell[0], gshells[0].data(),
                                                                     nshell[1], gshells[1].data());
        struct multishell_pair_flat Qf = create_multishell_pair_flat(nshell[2], gshells[2].data(),
                                                                     nshell[3], gshells[3].data());

        // acutally calculate
        Integral_FO(P, Q, res_FO);
        Integral_vref(P, Q, res_vref);
        Integral_vref_flat(Pf, Qf, res_vref_flat);


        // chop
        const int arrlen = nshell1234 * ncart;
        Chop(res_valeev, arrlen);
        Chop(res_FO, arrlen);
        Chop(res_vref, arrlen);
        Chop(res_vref_flat, arrlen);


        // Analyze
        std::pair<double, double> err_FO        = CalcError(res_FO,        res_valeev,  arrlen);
        std::pair<double, double> err_vref      = CalcError(res_vref,      res_valeev,  arrlen);
        std::pair<double, double> err_vref_flat = CalcError(res_vref_flat, res_valeev,  arrlen);
        std::pair<double, double> err_erd       = CalcError(res_liberd,    res_valeev,  arrlen);

        printf("( %2d %2d | %2d %2d )    %10.3e  %10.3e  %10.3e  %10.3e    %10.3e  %10.3e  %10.3e  %10.3e\n", am[0], am[1], am[2], am[3],
                                                      err_FO.first, err_vref.first, err_vref_flat.first, err_erd.first,
                                                      err_FO.second, err_vref.second, err_vref_flat.second, err_erd.second);


        // For debugging
        //for(int i = 0; i < ncart * nshell1234; i++)
        //    printf("%25.15e  %25.15e  %25.15e\n", res_valeev[i], res_liberd[i], res_vref[i]);

        free_multishell_pair(P);
        free_multishell_pair(Q);
        free_multishell_pair_flat(Pf);
        free_multishell_pair_flat(Qf);
        FreeQuartets(gshells);
        FreeQuartets(gshells_erd);

    }

    printf("\n");

    FreeShellMap(shellmap);

    Valeev_Finalize();
    Boys_Finalize();
    ERD_Finalize();

    FREE(res_valeev);
    FREE(res_liberd);
    FREE(res_vref);
    FREE(res_vref_flat);
    FREE(res_FO);

    return 0;
}
