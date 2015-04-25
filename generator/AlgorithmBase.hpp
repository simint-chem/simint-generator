#ifndef ALGORITHMBASE_HPP
#define ALGORITHMBASE_HPP

#include "generator/Classes.hpp"

class HRR_Algorithm_Base
{
    public:
        virtual ~HRR_Algorithm_Base() { }

        virtual HRRStep step(const Quartet & q, DoubletType steptype) = 0;
};


#endif
