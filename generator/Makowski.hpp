#ifndef MAKOWSKI_HPP
#define MAKOWSKI_HPP

#include <stdexcept>
#include <algorithm>

#include "generator/AlgorithmBase.hpp" 

class Makowski_HRR : public HRR_Algorithm_Base
{
    public:
        virtual HRRDoubletStep doubletstep(const Doublet & target)
        {
            if(target.am() == 0)
                throw std::runtime_error("Cannot HRR step to an s doublet!");

            // idx is the xyz index
            ExpList ijk = target.right.ijk;
            std::sort(ijk.begin(), ijk.end());
            auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
            auto it = std::find(target.right.ijk.rbegin(), target.right.ijk.rend(), *v); 
            int idx = 2 - std::distance(target.right.ijk.rbegin(), it);  // remember we are working with reverse iterators

            // gaussian common to both src1 and src2
            Gaussian common(target.right.StepDown(idx, 1));

            // for src1
            Gaussian src1g(target.left.StepUp(idx, 1));

            // create new doublets
            Doublet src1d{target.type, src1g, common};
            Doublet src2d{target.type, target.left, common};
           
            XYZStep xyzstep = IdxToXYZStep(idx);

            // Create the HRR doublet step
            HRRDoubletStep hrr{target, 
                               src1d, src2d,
                               xyzstep};
            return hrr;
        }
};

class Makowski_VRR : public VRR_Algorithm_Base
{
    public:
        virtual VRRMap CreateVRRMap(int am)
        {
            VRRMap vm;

            if(am == 0)
                return vm;  // empty!

            Gaussian g{am, 0, 0};

            do {      
                // lowest exponent, favoring the far right if equal
                // idx is the xyz index
                ExpList ijk = g.ijk;
                std::sort(ijk.begin(), ijk.end());
                auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
                auto it = std::find(g.ijk.rbegin(), g.ijk.rend(), *v); 
                int idx = 2 - std::distance(g.ijk.rbegin(), it);  // remember we are working with reverse iterators

                vm[g] = IdxToXYZStep(idx);
            } while(g.Iterate());

            return vm;
            
        }
};


#endif
