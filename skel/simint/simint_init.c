#include "simint/simint_init.h"
#include "simint/ostei/ostei_init.h"

void simint_init(void)
{
    simint_ostei_init();
}


void simint_finalize(void)
{
    simint_ostei_finalize();
}

