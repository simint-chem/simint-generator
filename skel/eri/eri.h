#ifndef ERI_H
#define ERI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "eri/gen/eri_generated.h"


// Disable this warning:
// remark #2620: attribute appears more than once
// No idea what it means, and from what I can tell the
// syntax is correct
#ifdef __INTEL_COMPILER
    #pragma warning(push)
    #pragma warning(disable:2620)
#endif

#include "eri/gen/hrr.h"
#include "eri/gen/vrr.h"

#ifdef __INTEL_COMPILER
    #pragma warning(pop)
#endif



#ifdef __cplusplus
}
#endif

#endif
