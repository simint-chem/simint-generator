#pragma once

#include "simint/shell/shell.h"
#include "simint/ostei/ostei_config.h"

#ifdef __cplusplus
#include "simint/cpp_restrict.hpp"
extern "C" {
#endif


int simint_compute_overlap(struct simint_shell * a,
                           struct simint_shell * b, 
                           double * restrict integrals);

int simint_compute_ke(struct simint_shell * a,
                      struct simint_shell * b, 
                      double * restrict integrals);

#ifdef __cplusplus
}
#endif
