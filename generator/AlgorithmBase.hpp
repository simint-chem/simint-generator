#ifndef ALGORITHMBASE_HPP
#define ALGORITHMBASE_HPP

#include "generator/Classes.hpp"

class HRR_Algorithm_Base
{
    public:
        virtual ~HRR_Algorithm_Base() { }

        virtual HRRDoubletStep doubletstep(const Doublet & target) = 0;
        virtual HRRQuartetStep quartetstep(const Quartet & target, DoubletType steptype)
        {
            // Call doublet step, then create the Quartet 

            Doublet d = target.get(steptype);

            // call doubletstep
            HRRDoubletStep hds = this->doubletstep(d);

            // Create the HRR step
            // (this preserves the flags of target
            if(steptype == DoubletType::BRA)
            {
                HRRQuartetStep hrr{target, 
                                  {hds.src1, target.ket, target.m, 0},   // src1 quartet
                                  {hds.src2, target.ket, target.m, 0},   // src2 quartet
                                  steptype,                                 // bra or ket being stepped
                                  hds.xyz};                                 // cartesian direction being stepped
                return hrr;
            }
            else
            {
                HRRQuartetStep hrr{target, 
                                  {target.bra, hds.src1, target.m, 0},   // src1 quartet
                                  {target.bra, hds.src2, target.m, 0},   // src2 quartet
                                  steptype,                                 // bra or ket being stepped
                                  hds.xyz};                                 // cartesian direction being stepped
                return hrr;
            }
        }
};


class VRR_Algorithm_Base
{
    public:
        virtual ~VRR_Algorithm_Base() { }
        virtual VRRMap CreateVRRMap(int am) = 0;

        // this will create a map for all possible
        // components, but hopefully only some
        // will be needed
        virtual VRRMap CreateAllMaps(int maxam)
        {
            VRRMap vm;
            for(int i = 0; i <= maxam; i++)
            {
                VRRMap vm2 = CreateVRRMap(i);
                vm.insert(vm2.begin(), vm2.end());
            }
           
            return vm; 
        }
};


#endif
