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

        virtual HRRDoubletStep DoubletStep_(const Doublet & target) = 0;

        void HRRDoubletLoop_(HRRDoubletStepList & hrrlist,
                             const DoubletSet & inittargets,
                             DoubletSet & solveddoublets,
                             DoubletSet & pruned);

        static void PruneDoublets_(DoubletSet & d, DoubletSet & pruned);
};


#endif
