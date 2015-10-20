#ifndef TEST_LIBINT2_HPP
#define TEST_LIBINT2_HPP

#include <vector>

#include "shell/shell.h"
#include "test/timer.h"

// Disable intel diagnostics for libint
// These happen in the libint header, so there's not much I can do about them
// 193 : zero used for undefined preprocessing identifier "LIBINT2_DEFINED__aB_s___0__s___1___TwoPRep_s___0__s___1___Ab__up_18
// 1418: external function definition with no prior declaration
//   82: storage class is not first
// 2259: non-pointer conversion from "double" to "float" may lose significant bits  
#ifdef __INTEL_COMPILER
    #pragma warning(push)
    #pragma warning(disable:193)
    #pragma warning(disable:1418)
    #pragma warning(disable:82)
    #pragma warning(disable:2259)
#endif
#include <libint2.h>
#ifdef __INTEL_COMPILER
    #pragma warning(pop)
#endif

class Libint2_ERI
{

    public:
        Libint2_ERI(int maxam, size_t maxnprim);
        ~Libint2_ERI();

        TimerType Integrals(struct multishell_pair P, struct multishell_pair Q, double * integrals);

    private:
        std::vector<Libint_eri_t> erival_;

};

#endif
