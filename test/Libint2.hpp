#ifndef TEST_LIBINT2_HPP
#define TEST_LIBINT2_HPP

#include "shell/shell.h"
#include "test/timer.hpp"

class Libint2_ERI
{
    public:
        TimerInfo Integrals(struct multishell_pair P, struct multishell_pair Q, double * integrals);

};

#endif
