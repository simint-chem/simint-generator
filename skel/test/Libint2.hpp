#ifndef TEST_LIBINT2_HPP
#define TEST_LIBINT2_HPP

#include "shell/shell.h"
#include "test/timer.hpp"
#include "test/libint2/libint2.h"

class Libint2_ERI
{

    public:
        Libint2_ERI(int maxam, size_t maxnprim);
        ~Libint2_ERI();

        TimerInfo Integrals(struct multishell_pair P, struct multishell_pair Q, double * integrals);

    private:
        std::vector<Libint_eri_t> erival_;

};

#endif
