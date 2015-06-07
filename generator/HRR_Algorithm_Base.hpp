#ifndef HRR_ALGORITHM_BASE_HPP
#define HRR_ALGORITHM_BASE_HPP

#include "generator/Classes.hpp"

class HRR_Algorithm_Base
{
    public:
        HRRBraKetStepList Create_DoubletStepLists(QAM amlist);

        std::pair<DAMSet, DAMSet> TopBraKetAM(void) const; 

        QAMSet TopQAM(void) const;
        QuartetSet TopQuartets(void) const;

        virtual ~HRR_Algorithm_Base() = default;

    private:
        DoubletSet bratop_, kettop_;
        DAMSet bratopam_, kettopam_;

        QAMSet topqam_;
        QuartetSet topquartets_;

        virtual HRRDoubletStep doubletstep(const Doublet & target) = 0;

        void HRRDoubletLoop(HRRDoubletStepList & hrrlist,
                            const DoubletSet & inittargets,
                            DoubletSet & solveddoublets);
};


#endif
