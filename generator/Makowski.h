#ifndef MAKOWSKI_H
#define MAKOWSKI_H

#include <stdexcept>
#include <algorithm>

#include "AlgorithmBase.h" 

class Makowski_HRR : public HRR_Algorithm_Base
{
    public:
        virtual HRRStep operator() (const Quartet & target, DoubletType steptype)
        {
            //Makowski: Recurse in the lowest non-zero angular
            //component of the second function
            Doublet const * dptr;

            if(steptype == DoubletType::BRA)
                dptr = &target.bra;
            else
                dptr = &target.ket;

            if(dptr->am() == 0)
                throw std::runtime_error("Cannot HRR step to an s doublet!");

            // idx is the xyz index
            std::array<int, 3> ijk = dptr->right.ijk;
            std::sort(ijk.begin(), ijk.end());
            auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
            auto it = std::find(dptr->right.ijk.rbegin(), dptr->right.ijk.rend(), *v); 
            int idx = 2 - std::distance(dptr->right.ijk.rbegin(), it);  // remember we are working with reverse iterators

            // gaussian common to both src1 and src2
            Gaussian common(dptr->right);
            common.ijk[idx] -= 1;

            // the left for src1
            Gaussian src1g(dptr->left);
            src1g.ijk[idx] += 1;

            // create new doublets
            Doublet src1d{steptype, src1g, common};
            Doublet src2d{steptype, dptr->left, common};
           
            XYZStep xyzstep = IdxToXYZStep(idx);
            std::cout << "IDX: " << idx << "\n";

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
