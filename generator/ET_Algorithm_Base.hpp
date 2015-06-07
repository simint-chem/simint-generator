#ifndef ET_ALGORITHMBASE_HPP
#define ET_ALGORITHMBASE_HPP

#include "generator/Classes.hpp"

class ET_Algorithm_Base
{
    public:
        virtual ETStepList Create_ETStepList(const QuartetSet & inittargets);

        QAMSet TopAM(void) const;

        virtual ~ET_Algorithm_Base() = default;

    private:
        QuartetSet ettop_;
        QAMSet ettopam_;

        virtual ETStep etstep(const Quartet & target) = 0;

        virtual void ETStepLoop(ETStepList & etsl,
                                const QuartetSet & inittargets,
                                QuartetSet & solvedquartets);

        void ETAddWithDependencies(std::vector<QAM> & amorder, QAM am);

};


#endif
