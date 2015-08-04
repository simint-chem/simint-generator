#ifndef SIMINT_TEST_COMMON_HPP
#define SIMINT_TEST_COMMON_HPP

#include <array>
#include <vector>
#include <string>
#include <map>
#include <fstream>
    
#include "shell/shell.h"
#include "test/timer.hpp"

#define MAXAM 2

#define RND_ZERO 1.0e-15

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

typedef std::vector<gaussian_shell> GaussianVec;
typedef std::map<int, GaussianVec> ShellMap;

// Function pointer typedefs
typedef int (*erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);
typedef int (*eriflatfunc)(struct multishell_pair_flat const, struct multishell_pair_flat const, double * const restrict);


// Some gaussian & quartet handling stuff
bool IsValidGaussian(const std::array<int, 3> & g);
bool IterateGaussian(std::array<int, 3> & g);
bool ValidQuartet(std::array<int, 4> am);
bool ValidQuartet(int i, int j, int k, int l);


// Deep copying and memory management
GaussianVec CopyGaussianVec(const GaussianVec & v);
void FreeGaussianVec(GaussianVec & agv);

ShellMap CopyShellMap(const ShellMap & m);
void FreeShellMap(ShellMap & m);


// Reads a basis set from a file
ShellMap ReadBasis(const std::string & file);

// Max parameters (maxam, maxnprim, maxel**4) in a basis set
std::array<int, 3> FindMapMaxParams(const ShellMap & m);


// Calculating reference integrals
void ValeevIntegrals(const GaussianVec & g1, const GaussianVec & g2,
                     const GaussianVec & g3, const GaussianVec & g4,
                     double * const integrals, bool normalize);


void ERDIntegrals(const GaussianVec & g1, const GaussianVec & g2,
                  const GaussianVec & g3, const GaussianVec & g4,
                  double * const integrals);


// Error analysis
void Chop(double * const restrict calc, int ncalc);

std::pair<double, double> CalcError(double const * const restrict calc, double const * const restrict ref, int ncalc);



// Setting up function pointers and calculating integrals
// using my code
void Init_Test(void);

int eri_notyetimplemented(struct multishell_pair const P,
                          struct multishell_pair const Q,
                          double * const restrict dummy);

int eriflat_notyetimplemented(struct multishell_pair_flat const P,
                              struct multishell_pair_flat const Q,
                              double * const restrict dummy);


TimerInfo Integral(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals);


#endif
