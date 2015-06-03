#include "generator/VRRWriter.hpp"
#include "generator/WriterBase.hpp"


VRRWriter::VRRWriter(const VRRMap & vrrmap, const VRRReqMap & vrrreqmap)
          : vrrmap_(vrrmap), vrrreqmap_(vrrreqmap)
{ }



bool VRRWriter::HasVRR(void) const
{
    return ( (vrrreqmap_.size() > 1) || (vrrreqmap_.size() == 1 && vrrreqmap_.begin()->first != 0) );
}



void VRRWriter::WriteIncludes(std::ostream & os, const WriterBase & base) const
{
    if(base.GetOption(OPTION_INLINEVRR) == 0)
        os << "#include \"eri/vrr/vrr.h\"\n";
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



void VRRWriter::WriteVRRSteps_(std::ostream & os, const WriterBase & base, const GaussianSet & greq, const std::string & num_m) const
{
    std::string indent1(20, ' ');
    std::string indent2(24, ' ');

    int am = greq.begin()->am();
    os << indent1 << "// Forming " << base.AuxName(am) << "[" << num_m << " * " << NCART(am) << "];\n";

    //os << indent1 << "// Needed from this AM:\n";
    //for(const auto & it : greq)
    //    os << indent1 << "//    " << it << "\n";

    os << indent1 << "for(m = 0; m < " << num_m << "; m++)  // loop over orders of auxiliary function\n";
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
}



void VRRWriter::WriteVRRInline_(std::ostream & os, const WriterBase & base) const
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

        // Write out the steps
        std::stringstream ss;
        ss << (greq, base.L()-am+1);
        WriteVRRSteps_(os, base, greq, ss.str());

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
    std::string indent1(26, ' ');

    os << "//////////////////////////////////////////////\n";
    os << "// VRR functions\n";
    os << "//////////////////////////////////////////////\n";
    os << "\n";
    os << "#include \"vectorization.h\"\n";
    os << "\n";

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
        os << "#pragma omp declare simd simdlen(SIMD_LEN)\n";
        os << "void VRR_" << base.AuxName(am) << "(const int num_m,\n";
        os << indent1 << "const double P_PA_x, const double P_PA_y, const double P_PA_z,\n";
        os << indent1 << "const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,\n";
        os << indent1 << "const double a_over_p,";
        if(am > 1)
            os << " const double one_over_2p,";
        os << "\n"; 
        os << indent1 << "double * const restrict " << base.AuxName(am) << ",\n";
        os << indent1 << "double const * const restrict " << base.AuxName(am-1);
        if(am > 1)
        {
            os << ",\n";
            os << indent1 << "double const * const restrict " << base.AuxName(am-2);
        }
        
        os << ")\n";
        os << "{\n";

        os << "    int m = 0;\n";

        // Write out the steps
        WriteVRRSteps_(os, base, greq, "num_m");

        os << "}\n";
    }
}




void VRRWriter::WriteVRRHeaderFile(std::ostream & os, const WriterBase & base) const
{
    std::string indent1(24, ' ');
    std::string indent2(4, ' ');

    os << "#ifndef VRR__H\n";
    os << "#define VRR__H\n";
    os << "\n";
    os << "//////////////////////////////////////////////\n";
    os << "// VRR functions\n";
    os << "//////////////////////////////////////////////\n";
    os << "\n";
    os << "#include \"vectorization.h\"\n";
    os << "\n";

    // iterate over increasing am
    for(const auto & it3 : vrrreqmap_)
    {
        int am = it3.first;

        // don't do zero - no VRR!
        if(am == 0)
            continue;

        os << "\n\n\n";
        os << "// VRR to obtain " << base.AuxName(am) << "\n";
        os << "#pragma omp declare simd simdlen(SIMD_LEN)\n";
        os << "void VRR_" << base.AuxName(am) << "(const int num_m,\n";
        os << indent1 << "const double P_PA_x, const double P_PA_y, const double P_PA_z,\n";
        os << indent1 << "const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,\n";
        os << indent1 << "const double a_over_p,";
        if(am > 1)
            os << " const double one_over_2p,";
        os << "\n"; 
        os << indent1 << "double * const restrict " << base.AuxName(am) << ",\n";
        os << indent1 << "double const * const restrict " << base.AuxName(am-1);
        if(am > 1)
        {
            os << ",\n";
            os << indent1 << "double const * const restrict " << base.AuxName(am-2);
        }
        
        os << ");\n";
    }

    os << "#endif\n";
    
}


void VRRWriter::WriteVRRExternal_(std::ostream & os, const WriterBase & base) const
{
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Vertical recurrance\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "\n";


    std::string indent1(20, ' ');
    std::string indent2(40, ' ');

    // iterate over increasing am
    for(const auto & it3 : vrrreqmap_)
    {
        int am = it3.first;

        // don't do zero - that is handled by the boys function stuff
        if(am == 0)
            continue;

        // greq is what is actually required from this am
        const GaussianSet & greq = it3.second;

        // call the function
        os << indent1 << "VRR_" << base.AuxName(am) << "(" << (base.L()-am+1) << ",\n";
        os << indent2 << "P_PA_x, P_PA_y, P_PA_z, aop_PQ_x, aop_PQ_y, aop_PQ_z,\n";
        os << indent2 << "a_over_p,";
        if(am > 1)
            os << "one_over_2p, ";
        os << "\n"; 

        os << base.AuxName(am) << ",\n";
        os << indent2 << base.AuxName(am-1);
        if(am > 1)
        {
            os << ",\n";
            os << indent2 << base.AuxName(am-2);
        }

        os << ");\n";


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



void VRRWriter::WriteVRR(std::ostream & os, const WriterBase & base) const
{
    if(base.GetOption(OPTION_INLINEVRR) > 0)
        WriteVRRInline_(os, base);
    else
        WriteVRRExternal_(os, base);
}

