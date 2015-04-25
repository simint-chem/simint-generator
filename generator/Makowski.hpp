#ifndef MAKOWSKI_HPP
#define MAKOWSKI_HPP

#include <stdexcept>
#include <algorithm>

#include "generator/AlgorithmBase.hpp" 

class Makowski_HRR : public HRR_Algorithm_Base
{
    public:
        virtual HRRStep step(const Quartet & target, DoubletType steptype)
        {
            //Makowski: Recurse in the lowest non-zero angular
            //component of the second function
            Doublet d = target.get(steptype);

            if(d.am() == 0)
                throw std::runtime_error("Cannot HRR step to an s doublet!");

            // idx is the xyz index
            std::array<int, 3> ijk = d.right.ijk;
            std::sort(ijk.begin(), ijk.end());
            auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
            auto it = std::find(d.right.ijk.rbegin(), d.right.ijk.rend(), *v); 
            int idx = 2 - std::distance(d.right.ijk.rbegin(), it);  // remember we are working with reverse iterators

            // gaussian common to both src1 and src2
            Gaussian common(d.right);
            common.ijk[idx] -= 1;

            // the left for src1
            Gaussian src1g(d.left);
            src1g.ijk[idx] += 1;

            // create new doublets
            Doublet src1d{steptype, src1g, common};
            Doublet src2d{steptype, d.left, common};
           
            XYZStep xyzstep = IdxToXYZStep(idx);

            // Create the HRR step
            if(steptype == DoubletType::BRA)
            {
                HRRStep hrr{target, 
                           {src1d, target.ket, target.m},   // src1 quartet
                           {src2d, target.ket, target.m},   // src2 quartet
                           steptype,                        // bra or ket being stepped
                           xyzstep};                        // cartesian direction being stepped
                return hrr;
            }
            else
            {
                HRRStep hrr{target, 
                           {target.bra, src1d, target.m},   // src1 quartet
                           {target.bra, src2d, target.m},   // src2 quartet
                           steptype,                        // bra or ket being stepped
                           xyzstep};                        // cartesian direction being stepped
                return hrr;
            }
        }
};

#endif
