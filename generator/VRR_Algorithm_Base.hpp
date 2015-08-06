#ifndef VRR_ALGORITHM_BASE_HPP
#define VRR_ALGORITHM_BASE_HPP

#include "generator/Classes.hpp"

class VRR_Algorithm_Base
{
    public:
        void Create(const QuartetSet & q);
        void Create(const QAM & q);

        VRRMap GetVRRMap(void) const;
        int GetMaxFm(void) const;
        VRRMReq GetVRRMReq(void) const;
        std::map<QAM, QAMSet> GetQAMReq(void) const;
        std::map<QAM, std::set<int>> GetIntReq_2p(void) const; 
        std::map<QAM, std::set<int>> GetIntReq_2q(void) const; 
        std::map<QAM, std::set<int>> GetIntReq_2pq(void) const; 

        std::map<QAM, std::set<std::string>> GetVarReq(void) const;
        std::set<std::string> GetAllVarReq(void) const;

        int GetMaxInt(void) const;

        virtual ~VRR_Algorithm_Base() = default; 

    private:
        // VRRMap maps a AM quartet to its steps
        VRRMap vrrmap_;

        // Maximum m value needed for a quartet
        VRRMReq vrrmreq_;

        // QAM required for a given QAM
        std::map<QAM, QAMSet> qamreq_;

        // Factors, etc, required for a given QAM
        std::map<QAM, std::set<std::string>> varreq_;
        std::map<QAM, std::set<int>> qamint_2p_, qamint_2q_, qamint_2pq_;

        // max int constant needed
        int maxint_;


        void PruneQuartets_(QuartetSet & q) const;
        virtual VRRStep VRRStep_(const Quartet & q) = 0;
};


#endif
