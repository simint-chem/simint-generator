#include "generator/VRRWriter.hpp"
#include "generator/WriterBase.hpp"


VRRWriter::VRRWriter(const VRRMap & vrrmap, const VRRReqMap & vrrreqmap)
          : vrrmap_(vrrmap), vrrreqmap_(vrrreqmap)
{ }



bool VRRWriter::HasVRR(void) const
{
    return ( (vrrreqmap_.size() > 1) || (vrrreqmap_.size() == 1 && vrrreqmap_.begin()->first != 0) );
}



void VRRWriter::DeclarePointers(std::ostream & os, const WriterBase & base) const
{
    if(vrrreqmap_.size() > 0)
    {
        os << "            // set up pointers to the contracted integrals - VRR\n";

        for(const auto & it : vrrreqmap_)
        {
            int vam = it.first;
            if(base.IsContArray({vam, 0, 0, 0}))
                os << "            double * const restrict PRIM_" << base.ArrVarName({vam, 0, 0, 0}) << " = "
                   << base.ArrVarName({vam, 0, 0, 0}) << " + (abcd * " << NCART(vam) << ");\n";
        }
    }
}



void VRRWriter::DeclareAuxArrays(std::ostream & os, const WriterBase & base) const
{
    if(vrrreqmap_.size())
    {
        os << "                    // Holds the auxiliary integrals ( i 0 | 0 0 )^m in the primitive basis\n";
        os << "                    // with m as the slowest index\n";

        for(const auto & greq : vrrreqmap_)
        {
            os << "                    // AM = " << greq.first << ": Needed from this AM: " << greq.second.size() << "\n";
            os << "                    double " << base.AuxName(greq.first) << "[" << (base.L()-greq.first+1) << " * " << NCART(greq.first) << "];\n";
            os << "\n";
        }

        os << "\n\n";

    }
}



void VRRWriter::WriteVRRInline(std::ostream & os, const WriterBase & base) const
{
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Vertical recurrance\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "\n";


    std::string indent1(20, ' ');
    std::string indent2(24, ' ');

    // iterate over increasing am
    for(const auto & it3 : vrrreqmap_)
    {
        int am = it3.first;

        // don't do zero - that is handled by the boys function stuff
        if(am == 0)
            continue;

        // greq is what is actually required from this am
        const GaussianSet & greq = it3.second;

        // requirements for this am
        os << indent1 << "// Forming " << base.AuxName(am) << "[" << (base.L()-am+1) << " * " << NCART(am) << "];\n";
        os << indent1 << "// Needed from this AM:\n";
        for(const auto & it : greq)
            os << indent1 << "//    " << it << "\n";
        os << indent1 << "for(m = 0; m < " << (base.L()-am+1) << "; m++)  // loop over orders of boys function\n";
        os << indent1 << "{\n";

        // iterate over the requirements
        // should be in order since it's a set
        for(const auto & it : greq)
        {
            // Get the stepping
            XYZStep step = vrrmap_.at(it);

            // and then step to the the required gaussians
            Gaussian g1 = it.StepDown(step, 1);
            Gaussian g2 = it.StepDown(step, 2);

            // the value of i in the VRR eqn
            // = value of exponent of g1 in the position of the step
            int vrr_i = g1.ijk[XYZStepToIdx(step)];

            os << indent2 << "//" << it <<  " : STEP: " << step << "\n";
            os << indent2 << base.AuxName(am) << "[m * " << NCART(am) << " + " << it.idx() << "] = P_PA_" << step << " * ";

            if(g1)
                os << base.AuxName(am-1) << "[m * " << NCART(am-1) << " + " << g1.idx() 
                   << "] - aop_PQ_" << step << " * " << base.AuxName(am-1) << "[(m+1) * " << NCART(am-1) << " + " << g1.idx() << "]";
            if(g2)
            {
                os << "\n"
                   << indent2 << "              + " << vrr_i 
                   << " * one_over_2p * ( " << base.AuxName(am-2) << "[m * " << NCART(am-2) << " +  " << g2.idx() << "]"
                   << " - a_over_p * " << base.AuxName(am-2) << "[(m+1) * " << NCART(am-2) << " + " << g2.idx() << "] )";
            }
            os << ";\n\n"; 
        }

        os << indent1 << "}\n";


        // if this target is also a contracted array, accumulate there
        if(base.IsContArray({am, 0, 0, 0}))
        {
            os << "\n";
            os << indent1 << "// Accumulating in contracted workspace\n";

            if(greq.size() == NCART(am))  // only do if wr calculated all of them?
            {
                os << indent1 << "for(n = 0; n < " << NCART(am) << "; n++)\n";
                os << indent2 << "PRIM_" << base.ArrVarName({am, 0, 0, 0}) << "[n] += " << base.AuxName(am) << "[n];\n";
            }
            else
            { 
                for(const auto & it : greq)
                    os << indent1 << "PRIM_" << base.AuxName(am) << "[" << it.idx() << "] += " << base.AuxName(am) << "[" << it.idx() << "];\n";
            }
        }

        os << "\n\n";
    }

    // accumulate INT__0_0_0_0 if needed
    if(base.IsContArray({0, 0, 0, 0}))
    {
        os << "\n";
        os << indent1 << "// Accumulating INT__0_0_0_0 in contracted workspace\n";
        os << indent1 << "*PRIM_" << base.ArrVarName({0, 0, 0, 0}) << " += *" << base.AuxName(0) << ";\n";
        os << "\n";
    }

    os << "\n";
}



void VRRWriter::WriteVRRFile(std::ostream & os, const WriterBase & base) const
{
    std::string indent1(24, ' ');
    std::string indent2(4, ' ');

    // iterate over increasing am
    for(const auto & it3 : vrrreqmap_)
    {
        int am = it3.first;

        // don't do zero - no VRR!
        if(am == 0)
            continue;

        // greq is what is actually required from this am
        const GaussianSet & greq = it3.second;

        os << "\n\n\n";
        os << "// VRR to obtain " << base.AuxName(am) << "\n";
        os << "void VRR_" << base.AuxName(am) << "(const int m,\n";
        os << indent1 << "const double P_PA_x, const double P_PA_y, const double P_PA_z,\n";
        os << indent1 << "const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,\n";
        os << indent1 << "const double a_over_p, const double one_over_2p,\n"; 
        os << indent1 << "double * const restrict " << base.AuxName(am) << ",\n";
        os << indent1 << "double const * const restrict " << base.AuxName(am-1);
        if(am > 1)
        {
            os << ",\n";
            os << indent1 << "double const * const restrict " << base.AuxName(am-2);
        }
        
        os << ")\n";
        os << "{\n";
        os << indent2 << "for(m = 0; m < " << (base.L()-am+1) << "; m++)  // loop over orders of boys function\n";
        os << indent2 << "{\n";

        // iterate over the requirements
        // should be in order since it's a set
        for(const auto & it : greq)
        {
            // Get the stepping
            XYZStep step = vrrmap_.at(it);

            // and then step to the the required gaussians
            Gaussian g1 = it.StepDown(step, 1);
            Gaussian g2 = it.StepDown(step, 2);

            // the value of i in the VRR eqn
            // = value of exponent of g1 in the position of the step
            int vrr_i = g1.ijk[XYZStepToIdx(step)];

            os << indent2 << "//" << it <<  " : STEP: " << step << "\n";
            os << indent2 << base.AuxName(am) << "[m * " << NCART(am) << " + " << it.idx() << "] = P_PA_" << step << " * ";

            if(g1)
                os << base.AuxName(am-1) << "[m * " << NCART(am-1) << " + " << g1.idx() 
                   << "] - aop_PQ_" << step << " * " << base.AuxName(am-1) << "[(m+1) * " << NCART(am-1) << " + " << g1.idx() << "]";
            if(g2)
            {
                os << "\n"
                   << indent2 << "              + " << vrr_i 
                   << " * one_over_2p * ( " << base.AuxName(am-2) << "[m * " << NCART(am-2) << " +  " << g2.idx() << "]"
                   << " - a_over_p * " << base.AuxName(am-2) << "[(m+1) * " << NCART(am-2) << " + " << g2.idx() << "] )";
            }
            os << ";\n\n"; 
        }

        os << "    }\n";
        os << "}\n";
    }
}

