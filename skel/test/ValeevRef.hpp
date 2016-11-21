#pragma once

#define MAXFAC 200
#define EPS 1e-17

// forward declaration
struct simint_shell;

// initialize math stuff
void ValeevRef_Init(void);
void ValeevRef_Finalize(void);


// Calculating reference integrals
void ValeevRef_Integrals(simint_shell const * const A, int nshell1,
                         simint_shell const * const B, int nshell2,
                         simint_shell const * const C, int nshell3,
                         simint_shell const * const D, int nshell4,
                         double * const integrals, int deriv, bool normalize);

