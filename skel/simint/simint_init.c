#include "simint/simint_init.h"
#include "simint/ostei/ostei_init.h"

void simint_init(void)
{
    simint_ostei_init();
    simint_ostei_deriv1_init();
}


void simint_finalize(void)
{
    simint_ostei_finalize();
    simint_ostei_deriv1_finalize();
}

