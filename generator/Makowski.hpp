#ifndef MAKOWSKI_HPP
#define MAKOWSKI_HPP

#include <stdexcept>
#include <algorithm>

#include "generator/VRR_Algorithm_Base.hpp" 
#include "generator/ET_Algorithm_Base.hpp" 
#include "generator/HRR_Algorithm_Base.hpp" 

class Makowski_HRR : public HRR_Algorithm_Base
{
    public:
        Makowski_HRR(const OptionMap & options)
            : HRR_Algorithm_Base(options)
        { }


    private:
        virtual HRRDoubletStep DoubletStep_(const Doublet & target, RRStepType steptype)
        {
            if(target.am() == 0)
                throw std::runtime_error("Cannot HRR step to an s doublet!");

            const Gaussian * t;
            const Gaussian * t2;

            if(steptype == RRStepType::I || steptype == RRStepType::K)  // Going J->I or L->K
            {
                t = &target.left;
                t2 = &target.right;
            }
            else // Going I->J or K->L
            {
                t = &target.right;
                t2 = &target.left;
            }
                

            // idx is the xyz index
            ExpList ijk = t->ijk;
            std::sort(ijk.begin(), ijk.end());
            auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
            auto it = std::find(t->ijk.rbegin(), t->ijk.rend(), *v); 
            int idx = 2 - std::distance(t->ijk.rbegin(), it);  // remember we are working with reverse iterators

            // gaussian common to both src1 and src2
            Gaussian common(t->StepDown(idx, 1));

            // for src1
            Gaussian src1g(t2->StepUp(idx, 1));

            // create new doublets
            Doublet src1d, src2d;
            if(steptype == RRStepType::I || steptype == RRStepType::K)
            {
                src1d = Doublet{target.type, common, src1g};
                src2d = Doublet{target.type, common, target.right};
            }
            else
            {
                src1d = Doublet{target.type, src1g, common};
                src2d = Doublet{target.type, target.left, common};
            }
           
            XYZStep xyzstep = IdxToXYZStep(idx);

            // Create the HRR doublet step
            HRRDoubletStep hrr{steptype,
                               target, 
                               src1d, src2d,
                               xyzstep};
            return hrr;
        }
};



class Makowski_VRR : public VRR_Algorithm_Base
{
    public:
        Makowski_VRR(const OptionMap & options)
            : VRR_Algorithm_Base(options)
        { }


    private:
        virtual VRRStep VRRStep_(const Quartet & q)
        {
            // if there is a ket part, we do that
            if(q.ket.am() > 0)
            {
                const Gaussian * g;
                VRRStep vs;

                // which center?
                if(q.ket.right.am() > 0)
                {
                    vs.type = RRStepType::L;
                    g = &q.ket.right;
                }
                else
                {
                    vs.type = RRStepType::K;
                    g = &q.ket.left;
                }

                // lowest exponent, favoring the far right if equal
                // idx is the xyz index
                ExpList ijk = g->ijk;
                std::sort(ijk.begin(), ijk.end());
                auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
                auto vit = std::find(g->ijk.rbegin(), g->ijk.rend(), *v); 
                int idx = 2 - std::distance(g->ijk.rbegin(), vit);  // remember we are working with reverse iterators


                vs.xyz = IdxToXYZStep(idx);
                vs.target = q;
                vs.src = {q, q, q, q, q, q, q, q};

                if(vs.type == RRStepType::K)
                {
                    vs.src[0].ket.left.ijk[idx]--;
                    vs.src[1].ket.left.ijk[idx]--;
                    vs.src[1].m++;

                    vs.src[2].ket.left.ijk[idx] -= 2;
                    vs.src[3].ket.left.ijk[idx] -= 2;
                    vs.src[3].m++;

                    vs.src[4].ket.right.ijk[idx]--;
                    vs.src[4].ket.left.ijk[idx]--;
                    vs.src[5].ket.right.ijk[idx]--;
                    vs.src[5].ket.left.ijk[idx]--;
                    vs.src[5].m++;

                    vs.src[6].bra.left.ijk[idx]--;
                    vs.src[6].ket.left.ijk[idx]--;
                    vs.src[6].m++;
                    vs.src[7].bra.right.ijk[idx]--;
                    vs.src[7].ket.left.ijk[idx]--;
                    vs.src[7].m++;
                }
                else
                {
                    vs.src[0].ket.right.ijk[idx]--;
                    vs.src[1].ket.right.ijk[idx]--;
                    vs.src[1].m++;

                    vs.src[2].ket.left.ijk[idx]--;
                    vs.src[2].ket.right.ijk[idx]--;
                    vs.src[3].ket.left.ijk[idx]--;
                    vs.src[3].ket.right.ijk[idx]--;
                    vs.src[3].m++;

                    vs.src[4].ket.right.ijk[idx] -= 2;
                    vs.src[5].ket.right.ijk[idx] -= 2;
                    vs.src[5].m++;

                    vs.src[6].bra.left.ijk[idx]--;
                    vs.src[6].ket.right.ijk[idx]--;
                    vs.src[6].m++;
                    vs.src[7].bra.right.ijk[idx]--;
                    vs.src[7].ket.right.ijk[idx]--;
                    vs.src[7].m++;
                }

                return vs;
            }
            else
            {
                const Gaussian * g;
                VRRStep vs;

                // which center?
                if(q.bra.right.am() > 0)
                {
                    vs.type = RRStepType::J;
                    g = &q.bra.right;
                }
                else
                {
                    vs.type = RRStepType::I;
                    g = &q.bra.left;
                }

                // lowest exponent, favoring the far right if equal
                // idx is the xyz index
                ExpList ijk = g->ijk;
                std::sort(ijk.begin(), ijk.end());
                auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
                auto vit = std::find(g->ijk.rbegin(), g->ijk.rend(), *v); 
                int idx = 2 - std::distance(g->ijk.rbegin(), vit);  // remember we are working with reverse iterators

                vs.xyz = IdxToXYZStep(idx);
                vs.target = q;
                vs.src = {q, q, q, q, q, q, q, q};

                if(vs.type == RRStepType::I)
                {
                    vs.src[0].bra.left.ijk[idx]--;
                    vs.src[1].bra.left.ijk[idx]--;
                    vs.src[1].m++;

                    vs.src[2].bra.left.ijk[idx] -= 2;
                    vs.src[3].bra.left.ijk[idx] -= 2;
                    vs.src[3].m++;

                    vs.src[4].bra.right.ijk[idx]--;
                    vs.src[4].bra.left.ijk[idx]--;
                    vs.src[5].bra.right.ijk[idx]--;
                    vs.src[5].bra.left.ijk[idx]--;
                    vs.src[5].m++;

                    vs.src[6].ket.left.ijk[idx]--;
                    vs.src[6].bra.left.ijk[idx]--;
                    vs.src[6].m++;
                    vs.src[7].ket.right.ijk[idx]--;
                    vs.src[7].bra.left.ijk[idx]--;
                    vs.src[7].m++;
                }
                else
                {
                    vs.src[0].bra.right.ijk[idx]--;
                    vs.src[1].bra.right.ijk[idx]--;
                    vs.src[1].m++;

                    vs.src[2].bra.left.ijk[idx]--;
                    vs.src[2].bra.right.ijk[idx]--;
                    vs.src[3].bra.left.ijk[idx]--;
                    vs.src[3].bra.right.ijk[idx]--;
                    vs.src[3].m++;

                    vs.src[4].bra.right.ijk[idx] -= 2;
                    vs.src[5].bra.right.ijk[idx] -= 2;
                    vs.src[5].m++;

                    vs.src[6].ket.left.ijk[idx]--;
                    vs.src[6].bra.right.ijk[idx]--;
                    vs.src[6].m++;
                    vs.src[7].ket.right.ijk[idx]--;
                    vs.src[7].bra.right.ijk[idx]--;
                    vs.src[7].m++;
                }

                return vs;
            }
        }
};



class Makowski_ET : public ET_Algorithm_Base
{
    public:
        Makowski_ET(const OptionMap & options)
            : ET_Algorithm_Base(options)
        { }


    private:
        virtual ETStep ETStep_(const Quartet & target)
        {
            if(target.am() == 0)
                throw std::runtime_error("Cannot ET step to an s doublet!");

            if(GetDirection() == DoubletType::KET)
            {
                // idx is the xyz index
                ExpList ijk = target.ket.left.ijk;
                std::sort(ijk.begin(), ijk.end());
                auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
                auto it = std::find(target.ket.left.ijk.rbegin(), target.ket.left.ijk.rend(), *v); 
                int idx = 2 - std::distance(target.ket.left.ijk.rbegin(), it);  // remember we are working with reverse iterators

                Gaussian src1g_1(target.bra.left);
                Gaussian src1g_2(target.ket.left.StepDown(idx, 1));

                Gaussian src2g_1(target.bra.left.StepDown(idx, 1));
                Gaussian src2g_2(target.ket.left.StepDown(idx, 1));

                Gaussian src3g_1(target.bra.left);
                Gaussian src3g_2(target.ket.left.StepDown(idx, 2));

                Gaussian src4g_1(target.bra.left.StepUp(idx, 1));
                Gaussian src4g_2(target.ket.left.StepDown(idx, 1));
                

                // create new doublets
                Gaussian s{0,0,0};

                Doublet src1d_bra{DoubletType::BRA, src1g_1, s};
                Doublet src1d_ket{DoubletType::KET, src1g_2, s};

                Doublet src2d_bra{DoubletType::BRA, src2g_1, s};
                Doublet src2d_ket{DoubletType::KET, src2g_2, s};

                Doublet src3d_bra{DoubletType::BRA, src3g_1, s};
                Doublet src3d_ket{DoubletType::KET, src3g_2, s};

                Doublet src4d_bra{DoubletType::BRA, src4g_1, s};
                Doublet src4d_ket{DoubletType::KET, src4g_2, s};


                XYZStep xyzstep = IdxToXYZStep(idx);

                // Create the electron transfer step
                ETStep et{target, 
                          {{
                            {src1d_bra, src1d_ket, 0},
                            {src2d_bra, src2d_ket, 0},
                            {src3d_bra, src3d_ket, 0},
                            {src4d_bra, src4d_ket, 0}
                          }},
                          {{
                            target.bra.left.ijk[idx],
                            target.ket.left.ijk[idx]-1,
                          }},
                          xyzstep
                         };

                return et;
            }
            else
            {
                // idx is the xyz index
                ExpList ijk = target.bra.left.ijk;
                std::sort(ijk.begin(), ijk.end());
                auto v = std::find_if(ijk.begin(), ijk.end(), [](int i) { return i != 0; });
                auto it = std::find(target.bra.left.ijk.rbegin(), target.bra.left.ijk.rend(), *v); 
                int idx = 2 - std::distance(target.bra.left.ijk.rbegin(), it);  // remember we are working with reverse iterators

                Gaussian src1g_1(target.bra.left.StepDown(idx, 1));
                Gaussian src1g_2(target.ket.left);

                Gaussian src2g_1(target.bra.left.StepDown(idx, 2));
                Gaussian src2g_2(target.ket.left.StepDown(idx, 0));

                Gaussian src3g_1(target.bra.left.StepDown(idx, 1));
                Gaussian src3g_2(target.ket.left.StepDown(idx, 1));

                Gaussian src4g_1(target.bra.left.StepDown(idx, 1));
                Gaussian src4g_2(target.ket.left.StepUp(idx, 1));
                

                // create new doublets
                Gaussian s{0,0,0};

                Doublet src1d_bra{DoubletType::BRA, src1g_1, s};
                Doublet src1d_ket{DoubletType::KET, src1g_2, s};

                Doublet src2d_bra{DoubletType::BRA, src2g_1, s};
                Doublet src2d_ket{DoubletType::KET, src2g_2, s};

                Doublet src3d_bra{DoubletType::BRA, src3g_1, s};
                Doublet src3d_ket{DoubletType::KET, src3g_2, s};

                Doublet src4d_bra{DoubletType::BRA, src4g_1, s};
                Doublet src4d_ket{DoubletType::KET, src4g_2, s};


                XYZStep xyzstep = IdxToXYZStep(idx);

                // Create the electron transfer step
                ETStep et{target, 
                          {{
                            {src1d_bra, src1d_ket, 0},
                            {src2d_bra, src2d_ket, 0},
                            {src3d_bra, src3d_ket, 0},
                            {src4d_bra, src4d_ket, 0}
                          }},
                          {{
                            target.bra.left.ijk[idx]-1,
                            target.ket.left.ijk[idx],
                          }},
                          xyzstep
                         };

                return et;
            }
        }
};

#endif
