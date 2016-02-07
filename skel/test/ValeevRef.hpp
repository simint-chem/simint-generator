#ifndef VALEEV_REF_HPP
#define VALEEV_REF_HPP

#define MAXFAC 200
#define EPS 1e-17

// forward declaration
struct gaussian_shell;

// initialize math stuff
void ValeevRef_Init(void);
void ValeevRef_Finalize(void);


// Calculating reference integrals
void ValeevRef_Integrals(gaussian_shell const * const A, int nshell1,
                         gaussian_shell const * const B, int nshell2,
                         gaussian_shell const * const C, int nshell3,
                         gaussian_shell const * const D, int nshell4,
                         double * const integrals, bool normalize);
#endif
