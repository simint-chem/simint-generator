#pragma once

#define MAXFAC 200
#define EPS 1e-17

// forward declaration
struct simint_shell;

// initialize math stuff
void ValeevRef_Init(void);
void ValeevRef_Finalize(void);


// Calculating reference integrals
void ValeevRef_Integrals(simint_shell const * const A, int nshellA,
                         simint_shell const * const B, int nshellB,
                         simint_shell const * const C, int nshellC,
                         simint_shell const * const D, int nshellD,
                         double * const integrals, int deriv, bool normalize);

