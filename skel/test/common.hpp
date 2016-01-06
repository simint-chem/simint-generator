#ifndef SIMINT_TEST_COMMON_HPP
#define SIMINT_TEST_COMMON_HPP

#include <array>
#include <vector>
#include <string>
#include <map>
#include <fstream>
    
#include "shell/shell.h"
#include "test/timer.h"

#define MAXAM 3

#define RND_ZERO 1.0e-15

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

typedef std::vector<gaussian_shell> GaussianVec;
typedef std::map<int, GaussianVec> ShellMap;



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
void ValeevIntegrals(gaussian_shell const * const A, int nshell1,
                     gaussian_shell const * const B, int nshell2,
                     gaussian_shell const * const C, int nshell3,
                     gaussian_shell const * const D, int nshell4,
                     double * const integrals, bool normalize);


// Error analysis
void Chop(double * const restrict calc, int ncalc);

std::pair<double, double> CalcError(double const * const restrict calc, double const * const restrict ref, int ncalc);




#endif
