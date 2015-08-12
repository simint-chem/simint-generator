#ifndef VRR_ALGORITHM_BASE_HPP
#define VRR_ALGORITHM_BASE_HPP

#include "generator/Classes.hpp"
#include "generator/Options.hpp"

typedef std::vector<VRRStep> VRR_StepList;

// Maps an AM quartet to the steps needed to create it
typedef std::map<QAM, VRR_StepList> VRR_StepMap;

// Gives the maximum m value for a given quartet
typedef std::map<QAM, int> VRR_MReqMap;

// Which AM are required for a given target AM
typedef std::map<QAM, QAMSet> VRR_AMReqMap;

// Maps integer constant requirements for a target AM
typedef std::map<QAM, IntSet> VRR_IntReqMap;

// Maps factors/variables (as strings) for a target AM
typedef std::map<QAM, StringSet> VRR_VarReqMap;

class VRR_Algorithm_Base
{
    public:
        VRR_Algorithm_Base(const OptionsMap & options);

        void Create(const QuartetSet & q);
        void Create(QAM q);

        QAM TargetAM(void) const;
        QAMList GetAMOrder(void) const;
        QAMSet GetAllAM(void) const;

        int GetMaxFm(void) const;

        VRR_StepList GetSteps(QAM am) const;
        int GetMReq(QAM am) const;
        QAMSet GetAMReq(QAM am) const;
        IntSet GetIntReq_2p(QAM am) const; 
        IntSet GetIntReq_2q(QAM am) const; 
        IntSet GetIntReq_2pq(QAM am) const; 

        IntSet GetAllInt_2p(void) const;
        IntSet GetAllInt_2q(void) const;
        IntSet GetAllInt_2pq(void) const;

        StringSet GetVarReq(QAM am) const;
        StringSet GetAllVarReq(void) const;

        bool HasBraVRR(void) const;
        bool HasKetVRR(void) const;

        int GetMaxInt(void) const;

        virtual ~VRR_Algorithm_Base() = default; 

    protected:
       int GetOption(int opt) const; 

    private:
        // Options
        OptionsMap options_;

        // VRR_StepMap maps a AM quartet to its steps
        VRR_StepMap vrrmap_;

        // Maximum m value needed for a quartet
        VRR_MReqMap vrrmreq_;

        // QAM required for a given QAM
        VRR_AMReqMap qamreq_;

        QAMList amorder_;

        // Factors, etc, required for a given QAM
        //  qamint_2p_ = constants multiplied by 1/2p
        //  qamint_2q_ = constants multiplied by 1/2q
        //  qamint_2pq_ = constants multiplied by 1/2(p+q)
        VRR_VarReqMap varreq_;
        VRR_IntReqMap qamint_2p_, qamint_2q_, qamint_2pq_;

        // max int constant needed
        int maxint_;

        QAMSet allam_;
        QAM targetam_;

        void PruneQuartets_(QuartetSet & q) const;
        void AMOrder_AddWithDependencies_(QAMList & order, QAM am) const;

        virtual VRRStep VRRStep_(const Quartet & q) = 0;
};


#endif
