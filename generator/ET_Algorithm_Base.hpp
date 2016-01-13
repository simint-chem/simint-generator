#ifndef ET_ALGORITHMBASE_HPP
#define ET_ALGORITHMBASE_HPP

#include "generator/Classes.hpp"
#include "generator/Options.hpp"

typedef std::vector<ETStep> ETStepList;

typedef std::map<QAM, ETStepList> ET_StepMap;

// Which AM are required for a given target AM
typedef std::map<QAM, QAMSet> ET_AMReqMap;

// Maps integer constant requirements for a target AM
typedef std::map<QAM, IntSet> ET_IntReqMap;


class ET_Algorithm_Base
{
    public:
        ET_Algorithm_Base(const OptionMap & options);

        void Create(const QuartetSet & inittargets, DoubletType direction);
        void Create(QAM am, DoubletType direction);

        ETStepList GetSteps(QAM am) const;

        QAMList GetAMOrder(void) const;
        QAMSet GetAllAM(void) const;
        QAMSet GetAMReq(QAM am) const;
        QuartetSet TopQuartets(void) const;
        QAMSet TopAM(void) const;

        IntSet GetIntReq(QAM am) const;
        IntSet GetAllInt(void) const;
        IntSet GetAllInt_p(void) const;
        IntSet GetAllInt_q(void) const;

        DoubletType GetDirection(void) const;

        bool HasET(void) const;
        bool HasBraET(void) const;
        bool HasKetET(void) const;

        virtual ~ET_Algorithm_Base() = default;

    protected:
       int GetOption(Option opt) const; 

    private:
        // Options
        OptionMap options_;

        DoubletType direction_;

        ET_StepMap etsteps_;
        ET_AMReqMap etreq_;
        QAMList amorder_;

        QAMSet allam_;

        ET_IntReqMap intreq_;
        IntSet allintreq_, allintreq_p_, allintreq_q_;

        QuartetSet ettop_;
        QAMSet ettopam_;

        virtual ETStep ETStep_(const Quartet & target) = 0;

        virtual void ETStepLoop_(ETStepList & etsl,
                                 const QuartetSet & inittargets,
                                 QuartetSet & solvedquartets, QuartetSet & pruned);

        void PruneQuartets_(QuartetSet & qs, QuartetSet & pruned) const;

        void AMOrder_AddWithDependencies_(std::vector<QAM> & amorder, QAM am) const;

};


#endif
