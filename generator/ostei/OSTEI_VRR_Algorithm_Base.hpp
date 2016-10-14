#pragma once


#include "generator/ostei/OSTEI_Types.hpp"
#include "generator/Options.hpp"

typedef std::set<VRRStep> VRR_StepSet;

// Maps an AM quartet to the steps needed to create it
typedef std::map<QAM, VRR_StepSet> VRR_StepMap;

// Gives the maximum m value for a given quartet
typedef std::map<QAM, int> VRR_MReqMap;

// Which AM are required for a given target AM
typedef std::map<QAM, QAMSet> VRR_AMReqMap;

// Maps integer constant requirements for a target AM
typedef std::map<QAM, IntSet> VRR_IntReqMap;

// Maps factors/variables (as strings) for a target AM
typedef std::map<QAM, StringSet> VRR_VarReqMap;

// Order we are doing the centers
typedef std::array<int, 4> IdxOrder;

class OSTEI_VRR_Algorithm_Base
{
    public:
        OSTEI_VRR_Algorithm_Base(const OptionMap & options);

        void Create(const QuartetSet & q);
        void Create(QAM q);
        void Create_WithOrder(const QuartetSet & q, IdxOrder idx_order);
        void Create_WithOrder(const QAM q, IdxOrder idx_order);

        QAMList GetAMOrder(void) const;
        QAMSet GetAllAM(void) const;
        RRStepType GetRRStep(QAM am) const;

        int GetMaxFm(void) const;

        VRR_StepSet GetSteps(QAM am) const;
        int GetMReq(QAM am) const;
        QAMList GetAMReq(QAM am) const;
        IntSet GetIntReq_2p(QAM am) const; 
        IntSet GetIntReq_2q(QAM am) const; 
        IntSet GetIntReq_2pq(QAM am) const; 

        IntSet GetAllInt_2p(void) const;
        IntSet GetAllInt_2q(void) const;
        IntSet GetAllInt_2pq(void) const;

        StringSet GetAllVarReq(void) const;
        StringSet GetVarReq(QAM am) const;

        QAMList GenerateAMReq(QAM am, RRStepType rrstep, bool all = false) const;
        StringSet GenerateVarReq(RRStepType rrstep) const;
        StringSet GenerateVarReq(QAM am, RRStepType rrstep) const;

        bool HasVRROfType(RRStepType steptype) const;

        bool HasVRR(void) const;
        bool HasBraVRR(void) const;
        bool HasKetVRR(void) const;
        bool HasVRR_I(void) const;
        bool HasVRR_J(void) const;
        bool HasVRR_K(void) const;
        bool HasVRR_L(void) const;

        int GetMaxInt(void) const;

        virtual ~OSTEI_VRR_Algorithm_Base() = default; 

    protected:
       int GetOption(Option opt) const;

    private:
        // Options
        OptionMap options_;

        // VRR_StepMap maps a AM quartet to its steps
        VRR_StepMap vrrmap_;

        // Maximum/minimum m value needed for a quartet
        VRR_MReqMap vrrmreq_max_;
        //VRR_MReqMap vrrmreq_min_;

        QAMList amorder_;

        // Factors, etc, required for a given QAM
        //  qamint_2p_ = constants multiplied by 1/2p
        //  qamint_2q_ = constants multiplied by 1/2q
        //  qamint_2pq_ = constants multiplied by 1/2(p+q)
        VRR_IntReqMap qamint_2p_, qamint_2q_, qamint_2pq_;

        // max int constant needed
        int maxint_;

        QAMSet allam_;

        void PruneQuartets_(QuartetSet & q) const;
        QuartetSet PruneQuartets_(QuartetSet & q, RRStepType rrstep) const;
        void AMOrder_AddWithDependencies_(QAMList & order, QAM am) const;

        virtual VRRStep VRRStep_(const Quartet & q, RRStepType rrstep) = 0;
        VRRStep Create_Single_(Quartet q, RRStepType rrstep);
};



