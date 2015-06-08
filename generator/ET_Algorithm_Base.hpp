#ifndef ET_ALGORITHMBASE_HPP
#define ET_ALGORITHMBASE_HPP

#include "generator/Classes.hpp"

class ET_Algorithm_Base
{
    public:
        virtual ETStepList Create_ETStepList(const QuartetSet & inittargets);

        QAMSet TopQAM(void) const;
        QuartetSet TopQuartets(void) const;
        GaussianMap TopGaussians(void) const;

        virtual ~ET_Algorithm_Base() = default;

    private:
        QuartetSet ettop_;
        QAMSet ettopam_;
        GaussianMap ettopgauss_;

        virtual ETStep ETStep_(const Quartet & target) = 0;

        virtual void ETStepLoop_(ETStepList & etsl,
                                const QuartetSet & inittargets,
                                QuartetSet & solvedquartets, QuartetSet & pruned);

        void ETAddWithDependencies_(std::vector<QAM> & amorder, QAM am);
        static void PruneQuartets_(QuartetSet & qs, QuartetSet & pruned);

};


#endif
