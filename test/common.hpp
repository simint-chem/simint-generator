#include <array>
#include <vector>
#include <string>
#include <map>
    
#include "test/cppvectorization.hpp"
#include "eri/shell.h"

#define MAXAM 2
#define NTEST 10

#define RND_ZERO 1.0e-15

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

typedef std::vector<gaussian_shell, AlignedAllocator<gaussian_shell>> AlignedGaussianVec;
typedef std::array<AlignedGaussianVec, 4> VecQuartet;
typedef std::map<int, AlignedGaussianVec> ShellMap;

typedef int (*erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);
typedef int (*eriflatfunc)(struct multishell_pair_flat const, struct multishell_pair_flat const, double * const restrict);

ShellMap ReadBasis(const std::string & file);


AlignedGaussianVec CopyAlignedGaussianVec(const AlignedGaussianVec & v);
VecQuartet CopyQuartets(const VecQuartet & orig);
ShellMap CopyShellMap(const ShellMap & m);
void FreeAlignedGaussianVec(AlignedGaussianVec & agv);
void FreeShellMap(ShellMap & m);
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


void Chop(double * const restrict calc, int ncalc);
std::pair<double, double> CalcError(double const * const restrict calc, double const * const restrict ref, int ncalc);

std::array<int, 3> FindMapMaxParams(const ShellMap & m);

int eri_notyetimplemented(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict dummy);

int eriflat_notyetimplemented(struct multishell_pair_flat const P,
                              struct multishell_pair_flat const Q,
                              double * const restrict dummy);


int Integral_FO(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals);
int Integral_vref(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals);
int Integral_vref_flat(struct multishell_pair_flat const P, struct multishell_pair_flat const Q, double * const restrict integrals);
