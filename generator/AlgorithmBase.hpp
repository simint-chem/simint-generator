#ifndef ALGORITHMBASE_H
#define ALGORITHMBASE_H

#include "generator/Classes.hpp"

class HRR_Algorithm_Base
{
    public:
        virtual ~HRR_Algorithm_Base() { }

        virtual HRRStep step(const Quartet & q, DoubletType steptype) = 0;
};


#endif
