#include "generator/VRR_Writer.hpp"
#include "generator/WriterBase.hpp"
#include "generator/VRR_Algorithm_Base.hpp"
#include "generator/Helpers.hpp"


VRR_Writer::VRR_Writer(const VRR_Algorithm_Base & vrr_algo)
{ 
    vrrmap_ = vrr_algo.GetVRRMap();
    vrramreq_ = vrr_algo.GetAMReq();
}


void VRR_Writer::WriteIncludes(std::ostream & os, const WriterBase & base) const
{
    if(base.GetOption(OPTION_INLINEVRR) == 0)
        os << "#include \"eri/vrr.gen/vrr.h\"\n";
}



void VRR_Writer::DeclarePrimArrays(std::ostream & os, const WriterBase & base) const
{
    if(vrramreq_.size())
    {
        os << indent6 << "// Holds the auxiliary integrals ( i 0 | 0 0 )^m in the primitive basis\n";
        os << indent6 << "// with m as the slowest index\n";

        for(const auto & greq : vrramreq_)
        {
            //os << indent4 << "// AM = " << greq.first << ": Needed from this AM: " << greq.second.size() << "\n";
            os << indent6 << base.DoubleType() << " " << base.PrimVarName({greq.first, 0, 0, 0})
               << "[" << (base.L()-greq.first+1) << " * " << NCART(greq.first) << "];\n";
            //os << "\n";
        }

        os << "\n\n";

    }
}


void VRR_Writer::DeclarePrimPointers(std::ostream & os, const WriterBase & base) const
{
    if(vrramreq_.size())
    {
        for(const auto & greq : vrramreq_)
        {
            QAM qam({greq.first, 0, 0, 0});
            if(base.IsContArray(qam))
                os << indent4 << "double * const restrict " << base.PrimPtrName(qam)
                   << " = " << base.ArrVarName(qam) << " + abcd * " << NCART(qam[0]) << ";\n";
        }

        os << "\n\n";

    }
}


void VRR_Writer::WriteVRRSteps_(std::ostream & os, const WriterBase & base, const GaussianSet & greq, const std::string & num_n) const
{
    int am = greq.begin()->am();
    QAM qam{am, 0, 0, 0};
    QAM qam1{am-1, 0, 0, 0};
    QAM qam2{am-2, 0, 0, 0};

    os << indent6 << "// Forming " << base.PrimVarName(qam) << "[" << num_n << " * " << NCART(am) << "];\n";

    //os << indent6 << "// Needed from this AM:\n";
    //for(const auto & it : greq)
    //    os << indent6 << "//    " << it << "\n";

    os << indent6 << "for(n = 0; n < " << num_n << "; ++n)  // loop over orders of auxiliary function\n";
    os << indent6 << "{\n";

    os << indent7 << "const int idx = n * " << NCART(am) << ";    // n * NCART(am)\n";
    os << indent7 << "const int idx1 = n * " << NCART(am-1) << ";    // n * NCART(am-1)\n"; 
    os << indent7 << "const int idx11 = idx1 + " << NCART(am-1) << ";    // (n+1) * NCART(am-1)\n"; 
    if(am > 1)
    {
        os << indent7 << "const int idx2 = n * " << NCART(am-2) << ";    // n * NCART(am-2)\n";
        os << indent7 << "const int idx21 = idx2 + " << NCART(am-2) << ";   // (n+1) * NCART(am-2)\n";
    }
        
    os << "\n";

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
        std::stringstream vrr_i;
        vrr_i << g1.ijk[XYZStepToIdx(step)];

        os << indent7 << "//" << it <<  " : STEP: " << step << "\n";
        os << indent7 << base.PrimVarName(qam) << "[idx + " << it.idx() << "] = P_PA_" << step << " * ";

        if(g1)
            os << base.PrimVarName(qam1) << "[idx1 + " << g1.idx() 
               << "] - aop_PQ_" << step << " * " << base.PrimVarName(qam1) << "[idx11 + " << g1.idx() << "]";
        if(g2)
        {
            os << "\n"
               << indent8 << "+ " << base.DoubleSet(vrr_i.str())
               << " * one_over_2p * ( " << base.PrimVarName(qam2) << "[idx2 + " << g2.idx() << "]"
               << " - a_over_p * " << base.PrimVarName(qam2) << "[idx21 + " << g2.idx() << "] )";
        }
        os << ";\n\n"; 
    }

    os << indent6 << "}\n";
}



void VRR_Writer::WriteVRRInline_(std::ostream & os, const WriterBase & base) const
{
    os << "\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << indent6 << "// Primitive integrals: Vertical recurrance\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << "\n";

    // iterate over increasing am
    for(const auto & it3 : vrramreq_)
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
        // (the check for contracted array happens in function)
        WriteAccumulate_(os, am, base);
        os << "\n\n";
    }

    // accumulate ssss if needed
    WriteAccumulate_(os, 0, base);
    os << "\n\n";
}



void VRR_Writer::WriteVRRFile(std::ostream & os, const WriterBase & base) const
{
    os << "//////////////////////////////////////////////\n";
    os << "// VRR functions\n";
    os << "//////////////////////////////////////////////\n";
    os << "\n";
    os << "#include \"vectorization/vectorization.h\"\n";
    os << "\n";

    // iterate over increasing am
    for(const auto & it3 : vrramreq_)
    {
        int am = it3.first;
        QAM qam{am, 0, 0, 0};
        QAM qam1{am-1, 0, 0, 0};
        QAM qam2{am-2, 0, 0, 0};

        // don't do zero - no VRR!
        if(am == 0)
            continue;

        // greq is what is actually required from this am
        const GaussianSet & greq = it3.second;


        os << "\n\n\n";
        os << "// VRR to obtain " << base.PrimVarName(qam) << "\n";
        os << "#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_n)\n";
        os << "void VRR_" << amchar[am] << "(const int num_n,\n";
        os << indent3 << "const double P_PA_x, const double P_PA_y, const double P_PA_z,\n";
        os << indent3 << "const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,\n";
        os << indent3 << "const double a_over_p,";
        if(am > 1)
            os << " const double one_over_2p,";
        os << "\n"; 
        os << indent3 << "double * const restrict " << base.PrimVarName(qam) << ",\n";
        os << indent3 << "double const * const restrict " << base.PrimVarName(qam1);
        if(am > 1)
        {
            os << ",\n";
            os << indent3 << "double const * const restrict " << base.PrimVarName(qam2);
        }
        
        os << ")\n";
        os << "{\n";

        os << "    int n = 0;\n";

        // Write out the steps
        WriteVRRSteps_(os, base, greq, "num_n");

        os << "}\n";
    }
}




void VRR_Writer::WriteVRRHeaderFile(std::ostream & os, const WriterBase & base) const
{
    os << "#ifndef VRR__H\n";
    os << "#define VRR__H\n";
    os << "\n";
    os << "//////////////////////////////////////////////\n";
    os << "// VRR functions\n";
    os << "//////////////////////////////////////////////\n";
    os << "\n";
    os << "#include \"vectorization/vectorization.h\"\n";
    os << "\n";

    // iterate over increasing am
    for(const auto & it3 : vrramreq_)
    {
        int am = it3.first;
        QAM qam{am, 0, 0, 0};
        QAM qam1{am-1, 0, 0, 0};
        QAM qam2{am-2, 0, 0, 0};

        // don't do zero - no VRR!
        if(am == 0)
            continue;

        os << "\n\n\n";
        os << "// VRR to obtain " << base.PrimVarName(qam) << "\n";
        os << "#pragma omp declare simd simdlen(SIMD_LEN) uniform(num_n)\n";
        os << "void VRR_" << amchar[am] << "(const int num_n,\n";
        os << indent3 << "const double P_PA_x, const double P_PA_y, const double P_PA_z,\n";
        os << indent3 << "const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,\n";
        os << indent3 << "const double a_over_p,";
        if(am > 1)
            os << " const double one_over_2p,";
        os << "\n"; 
        os << indent3 << "double * const restrict " << base.PrimVarName(qam) << ",\n";
        os << indent3 << "double const * const restrict " << base.PrimVarName(qam1);
        if(am > 1)
        {
            os << ",\n";
            os << indent3 << "double const * const restrict " << base.PrimVarName(qam2);
        }
        
        os << ");\n";
    }

    os << "#endif\n";
    
}



void VRR_Writer::WriteAccumulate_(std::ostream & os, int am, const WriterBase & base) const
{
    QAM qam{am, 0, 0, 0};

    // if this target is also a contracted array, accumulate there
    if(base.IsContArray(qam))
    {
        os << "\n";
        os << indent6 << "// Accumulating in contracted workspace\n";

        os << indent6 << "for(n = 0; n < " << NCART(am) << "; n++)\n";

        if(base.Intrinsics())
        {
            os << indent6 << "{\n";
            
            /*
            os << indent7 << "double vec[" << base.SimdLen() << "] __attribute__((aligned(" << base.ByteAlign() << ")));\n";
            os << indent7 << base.DoubleStore(base.PrimVarName(qam) + "[n]", "vec", "") << ";\n";
            os << indent7 << base.PrimPtrName(qam) << "[n] += vec[0]";
            */

            os << indent7 << "union double4 vec = (union double4)" << base.PrimVarName(qam) << "[n];\n";    
            os << indent7 << base.PrimPtrName(qam) << "[n] += vec.v[0]";

            for(int i = 1; i < base.SimdLen(); i++)
                os << " + vec.v[" << i << "]";
                //os << " + vec[" << i << "]";
            os << ";\n";
            os << indent6 << "}\n";
        }
        else
            os << indent7 << base.PrimPtrName(qam) << "[n] += " << base.PrimVarName(qam) << "[n];\n";

        // TODO
        // if I don't need all, then don't accumulate all
        //for(const auto & it : greq)
        //    os << indent1 << base.ArrVarName(qam) << "[" << it.idx() << " * nshell1234 + abcd" << "] += " << base.PrimVarName(qam) << "[" << it.idx() << "];\n";
    }
}


void VRR_Writer::WriteVRRExternal_(std::ostream & os, const WriterBase & base) const
{
    os << "\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << indent6 << "// Primitive integrals: Vertical recurrance\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << "\n";

    // iterate over increasing am
    for(const auto & it3 : vrramreq_)
    {
        int am = it3.first;
        QAM qam{am, 0, 0, 0};
        QAM qam1{am-1, 0, 0, 0};
        QAM qam2{am-2, 0, 0, 0};

        // don't do zero - that is handled by the boys function stuff
        if(am == 0)
            continue;

        // greq is what is actually required from this am
        //const GaussianSet & greq = it3.second;

        // call the function
        os << indent6 << "VRR_" << amchar[am] << "(" << (base.L()-am+1) << ",\n";
        os << indent7 << "P_PA_x, P_PA_y, P_PA_z, aop_PQ_x, aop_PQ_y, aop_PQ_z,\n";
        os << indent7 << "a_over_p,";
        if(am > 1)
            os << "one_over_2p, ";
        os << "\n"; 

        os << indent7 << base.PrimVarName(qam) << ",\n";
        os << indent7 << base.PrimVarName(qam1);
        if(am > 1)
        {
            os << ",\n";
            os << indent7 << base.PrimVarName(qam2);
        }

        os << ");\n";


        // if this target is also a contracted array, accumulate there
        WriteAccumulate_(os, am, base);
        os << "\n\n";
    }

    // accumulate ssss if needed
    WriteAccumulate_(os, 0, base);
    os << "\n\n";
}



void VRR_Writer::WriteVRR(std::ostream & os, const WriterBase & base) const
{
    if(base.GetOption(OPTION_INLINEVRR) > 0)
        WriteVRRInline_(os, base);
    else
        WriteVRRExternal_(os, base);
}

