#ifndef ET_ALGORITHMBASE_HPP
#define ET_ALGORITHMBASE_HPP

#include "generator/Classes.hpp"

typedef std::vector<ETStep> ETStepList;

typedef std::map<QAM, ETStepList> ET_StepMap;

// Which AM are required for a given target AM
typedef std::map<QAM, QAMSet> ET_AMReqMap;

// Maps integer constant requirements for a target AM
typedef std::map<QAM, IntSet> ET_IntReqMap;

class ET_Algorithm_Base
{
    public:
        void Create(const QuartetSet & inittargets);
        void Create(QAM am);

        ETStepList GetSteps(QAM am) const;

        QAMList GetAMOrder(void) const;
        QAMSet GetAllAM(void) const;
        QAMSet GetAMReq(QAM am) const;
        QuartetSet TopQuartets(void) const;
        QAMSet TopAM(void) const;

        IntSet GetIntReq(QAM am) const;
        IntSet GetAllInt(void) const;

        virtual ~ET_Algorithm_Base() = default;

    private:
        ET_StepMap etsteps_;
        ET_AMReqMap etreq_;
        QAMList amorder_;

        QAMSet allam_;

        ET_IntReqMap intreq_;
        IntSet allintreq_;

        QuartetSet ettop_;
        QAMSet ettopam_;

        virtual ETStep ETStep_(const Quartet & target) = 0;

        virtual void ETStepLoop_(ETStepList & etsl,
                                 const QuartetSet & inittargets,
                                 QuartetSet & solvedquartets, QuartetSet & pruned);

        static void PruneQuartets_(QuartetSet & qs, QuartetSet & pruned);

        void AMOrder_AddWithDependencies_(std::vector<QAM> & amorder, QAM am) const;

};


#endif
