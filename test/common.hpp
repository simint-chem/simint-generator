#include <array>
#include <vector>
#include <string>
    
#include "test/cppvectorization.hpp"
#include "eri/shell.h"

#define MAX_COORD 0.5
#define MAX_EXP   50.0
#define MAX_COEF  2.0

#define MAXAM 2

#define NTEST 10

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

typedef std::vector<gaussian_shell, AlignedAllocator<gaussian_shell>> AlignedGaussianVec;
typedef std::array<AlignedGaussianVec, 4> VecQuartet;

typedef int (*erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);
typedef int (*eriflatfunc)(struct multishell_pair_flat const, struct multishell_pair_flat const, double * const restrict);

/*
struct gaussian_shell random_shell(int nprim, int am);

VecQuartet CreateRandomQuartets(std::array<int, 4> nshell,
                                std::array<int, 4> nprim,
                                std::array<int, 4> am,
                                bool normalize = true);
*/


VecQuartet ReadQuartets(std::array<int, 4> am,
                        std::string basedir,
                        bool normalize);

VecQuartet CopyQuartets(const VecQuartet & orig);


void FreeQuartets(VecQuartet & arr);


bool IsValidGaussian(const std::array<int, 3> & g);

bool IterateGaussian(std::array<int, 3> & g);

int ReadValeevIntegrals(std::string basedir,
                        const std::array<int, 4> & am,
                        double * res);

void ValeevIntegrals(const VecQuartet & gshells, double * integrals, bool normalize);


void ERDIntegrals(const VecQuartet & gshells, double * integrals);


bool ValidQuartet(std::array<int, 4> am);


void Init_Test(void);


std::pair<double, double> CalcError(double const * const calc, double const * const ref, int ncalc);

int eri_notyetimplemented(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict dummy);

int eriflat_notyetimplemented(struct multishell_pair_flat const P,
                              struct multishell_pair_flat const Q,
                              double * const restrict dummy);


int Integral_FO(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals);
int Integral_vref(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals);
int Integral_vref_flat(struct multishell_pair_flat const P, struct multishell_pair_flat const Q, double * const restrict integrals);
