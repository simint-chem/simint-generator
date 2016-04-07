#include "simint/simint_init.h"
#include "simint/eri/eri_init.h"
#include "simint/boys/boys_init.h"

void simint_init(void)
{
    simint_eri_init();
    simint_boys_init();
}


void simint_finalize(void)
{
    simint_boys_finalize();
    simint_eri_finalize();
}

