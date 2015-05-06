#ifndef ALGORITHMBASE_HPP
#define ALGORITHMBASE_HPP

#include "generator/Classes.hpp"

class HRR_Algorithm_Base
{
    public:
        virtual ~HRR_Algorithm_Base() { }

        HRRBraKetStepList Create_DoubletStepLists(QAMList amlist);
        HRRQuartetStepList Create_QuartetStepList(QAMList amlist);

    private:
        virtual HRRDoubletStep doubletstep(const Doublet & target) = 0;
        virtual HRRQuartetStep quartetstep(const Quartet & target, DoubletType steptype);

        void HRRDoubletLoop(HRRDoubletStepList & hrrlist,
                            const DoubletSet & inittargets,
                            DoubletSet & solveddoublets);

        void HRRQuartetLoop(HRRQuartetStepList & hrrlist,
                            const QuartetSet & inittargets,
                            QuartetSet & solvedquartets,
                            DoubletType type);
};


class ET_Algorithm_Base
{
    public:
        virtual ETStepList Create_ETStepList(const QuartetSet & inittargets);

    private:    
        virtual ETStep etstep(const Quartet & target) = 0;
};


class VRR_Algorithm_Base
{
    public:
        virtual ~VRR_Algorithm_Base() { }
        virtual VRRMap CreateVRRMap(int am) = 0;

        // this will create a map for all possible
        // components, but hopefully only some
        // will be needed
        std::pair<VRRMap, VRRReqMap> CreateAllMaps(const GaussianSet & greq);
};


#endif
