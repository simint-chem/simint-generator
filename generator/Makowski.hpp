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
            Gaussian common(target.right);
            common.ijk[idx] -= 1;

            // for src1
            Gaussian src1g(target.left);
            src1g.ijk[idx] += 1;

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

#endif
