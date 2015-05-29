#include <array>
#include <vector>
#include <string>
#include <map>
#include <fstream>
    
#include "test/cppvectorization.hpp"
#include "eri/shell.h"

#define MAXAM 2
#define NTEST 10

#define RND_ZERO 1.0e-15

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

typedef std::vector<gaussian_shell, AlignedAllocator<gaussian_shell>> AlignedGaussianVec;
typedef std::map<int, AlignedGaussianVec> ShellMap;

// Function pointer typedefs
typedef int (*erifunc)(struct multishell_pair const, struct multishell_pair const, double * const restrict);
typedef int (*eriflatfunc)(struct multishell_pair_flat const, struct multishell_pair_flat const, double * const restrict);


// Class for reading reference integrals
class RefIntegralReader
{
    public:
        RefIntegralReader(const std::string & basfile);
        RefIntegralReader(const RefIntegralReader & rhs) = delete;
        RefIntegralReader(const RefIntegralReader && rhs) = delete;

        ~RefIntegralReader() = default;

        void ReadNext(double * out, int nsize);
        void Reset(void);

    private:
        std::ifstream file_;
};



// Some gaussian & quartet handling stuff
bool IsValidGaussian(const std::array<int, 3> & g);
bool IterateGaussian(std::array<int, 3> & g);
bool ValidQuartet(std::array<int, 4> am);
bool ValidQuartet(int i, int j, int k, int l);


// Deep copying and memory management
AlignedGaussianVec CopyAlignedGaussianVec(const AlignedGaussianVec & v);
void FreeAlignedGaussianVec(AlignedGaussianVec & agv);

ShellMap CopyShellMap(const ShellMap & m);
void FreeShellMap(ShellMap & m);


// Reads a basis set from a file
ShellMap ReadBasis(const std::string & file);

// Max parameters (maxam, maxnprim, maxel**4) in a basis set
std::array<int, 3> FindMapMaxParams(const ShellMap & m);


// Reads reference integrals from a file
int ReadValeevIntegrals(std::string basfile,
                        double * res);


// Calculating reference integrals
void ValeevIntegrals(const AlignedGaussianVec & g1, const AlignedGaussianVec & g2,
                     const AlignedGaussianVec & g3, const AlignedGaussianVec & g4,
                     double * const integrals, bool normalize);


void ERDIntegrals(const AlignedGaussianVec & g1, const AlignedGaussianVec & g2,
                  const AlignedGaussianVec & g3, const AlignedGaussianVec & g4,
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


int Integral_FO(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals);
int Integral_vref(struct multishell_pair const P, struct multishell_pair const Q, double * const restrict integrals);
int Integral_vref_flat(struct multishell_pair_flat const P, struct multishell_pair_flat const Q, double * const restrict integrals);
