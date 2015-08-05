#ifndef ET_ALGORITHMBASE_HPP
#define ET_ALGORITHMBASE_HPP

#include "generator/Classes.hpp"

class ET_Algorithm_Base
{
    public:
        virtual void Create_ETStepList(const QuartetSet & inittargets);

        ETStepList ETSteps(void) const;
        QAMSet TopQAM(void) const;
        QuartetSet TopQuartets(void) const;

        virtual ~ET_Algorithm_Base() = default;

    private:
        ETStepList etsteps_;
        QuartetSet ettop_;
        QAMSet ettopam_;

        virtual ETStep ETStep_(const Quartet & target) = 0;

        virtual void ETStepLoop_(ETStepList & etsl,
                                const QuartetSet & inittargets,
                                QuartetSet & solvedquartets, QuartetSet & pruned);

        void ETAddWithDependencies_(std::vector<QAM> & amorder, QAM am);
        static void PruneQuartets_(QuartetSet & qs, QuartetSet & pruned);

};


#endif
