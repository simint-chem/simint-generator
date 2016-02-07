#ifndef TEST_LIBINT2_HPP
#define TEST_LIBINT2_HPP

#include <vector>

#include "simint/shell/shell.h"
#include "test/Timer.h"

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

        Libint2_ERI(const Libint2_ERI & rhs) = delete;
        Libint2_ERI(Libint2_ERI && rhs)      = default;
        Libint2_ERI & operator=(const Libint2_ERI & rhs) = delete;
        Libint2_ERI & operator=(Libint2_ERI && rhs)      = delete;

        ~Libint2_ERI();

        std::pair<TimerType, TimerType>
        Integrals(struct multishell_pair P,
                  struct multishell_pair Q,
                  double * integrals);

    private:
        std::vector<Libint_eri_t> erival_;

        TimerType IntegralsScalar_(struct multishell_pair P, struct multishell_pair Q, double * integrals);

        #ifdef TESTS_LIBINT2_SIMD
        TimerType IntegralsSIMD_(struct multishell_pair P, struct multishell_pair Q, double * integrals);

        size_t worksize_;
        double *work_;
        double *tmp_AB_x_, *tmp_AB_y_, *tmp_AB_z_;
        double *tmp_CD_x_, *tmp_CD_y_, *tmp_CD_z_;
        double *tmp_PA_x_, *tmp_PA_y_, *tmp_PA_z_;
        double *tmp_QC_x_, *tmp_QC_y_, *tmp_QC_z_;
        double *tmp_WP_x_, *tmp_WP_y_, *tmp_WP_z_;
        double *tmp_WQ_x_, *tmp_WQ_y_, *tmp_WQ_z_;
        double *tmp_oo2z_, *tmp_oo2e_, *tmp_oo2ze_, *tmp_roz_, *tmp_roe_;
        double *tmp_vecF_;
        #endif

};

#endif
