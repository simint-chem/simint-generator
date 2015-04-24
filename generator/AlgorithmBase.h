#ifndef ALGORITHM_BASE_H
#define ALGORITHM_BASE_H

#include "Classes.h"

class HRR_Algorithm_Base
{
    public:
        virtual ~HRR_Algorithm_Base() { }

        virtual HRRStep operator() (const Quartet & q, DoubletType steptype) = 0;
};


#endif
