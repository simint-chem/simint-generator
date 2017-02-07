#pragma once


#include "generator/ostei/OSTEI_Types.hpp"
#include "generator/Options.hpp"


typedef std::vector<HRRDoubletStep> HRRDoubletStepList;
typedef std::map<DAM, HRRDoubletStepList> HRR_StepMap;
typedef std::map<DAM, DAMSet> HRR_AMReqMap;

class OSTEI_HRR_Algorithm_Base
{
    public:
        OSTEI_HRR_Algorithm_Base(const OptionMap & options);

        void Create(QAM am);
        void Create(std::set<QAM> am);
        void Create(QAM am, RRStepType brasteptype, RRStepType ketsteptype);
        void Create(std::set<QAM> am, RRStepType brasteptype, RRStepType ketsteptype);
        
        QAMSet TopAM(void) const;
        QuartetSet TopQuartets(void) const;
        QAMList GetAMOrder(void) const;

        RRStepType GetBraRRStep(DAM am) const;
        RRStepType GetKetRRStep(DAM am) const;
        DoubletType GetDoubletStep(QAM am) const;

        DAMSet GetBraAMReq(DAM am) const;
        DAMSet GetKetAMReq(DAM am) const;

        QAMList GenerateAMReq(QAM am, RRStepType rrstep) const;

        HRRDoubletStepList GetBraSteps(DAM am) const;
        HRRDoubletStepList GetKetSteps(DAM am) const;

        bool HasHRR(void) const;
        bool HasBraHRR(void) const;
        bool HasKetHRR(void) const;


        virtual ~OSTEI_HRR_Algorithm_Base() = default;

    protected:
       int GetOption(Option opt) const; 

    private:
        // Options
        OptionMap options_;

        DAMSet allbraam_, allketam_;

        HRR_StepMap brasteps_, ketsteps_;
        HRR_AMReqMap brareq_, ketreq_;

        QuartetSet topquartets_;
        QAMSet topqam_;

        QAMList amorder_;
        std::map<QAM, DoubletType> stepmap_;

        virtual HRRDoubletStep DoubletStep_(const Doublet & target, RRStepType steptype) = 0;

        void HRRDoubletLoop_(HRRDoubletStepList & hrrlist,
                             const DoubletSet & inittargets,
                             DoubletSet & solveddoublets,
                             DoubletSet & pruned,
                             RRStepType steptype);

        static void PruneDoublets_(DoubletSet & d, DoubletSet & pruned, RRStepType steptype);

        void AMOrder_AddWithDependencies_(QAMList & order,
                                          QAMSet & topqam,
                                          std::map<QAM, DoubletType> & smap,
                                          QAMSet am) const;
 
};



