#include "generator/VRR_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/VRR_Algorithm_Base.hpp"
#include "generator/Helpers.hpp"
#include "generator/Ncart.hpp"


VRR_Writer::VRR_Writer(const VRR_Algorithm_Base & vrr_algo)
    : vrr_algo_(vrr_algo)
{ 
}


void VRR_Writer::AddConstants(void) const
{
    if(WriterInfo::GetOption(OPTION_INLINEVRR) > 0)
    {
        for(int i = 1; i <= vrr_algo_.GetMaxInt(); i++)
            WriterInfo::AddIntConstant(i);
    }
}


void VRR_Writer::DeclarePrimArrays(std::ostream & os) const
{
    os << indent5 << "// Holds the auxiliary integrals ( i 0 | k 0 )^m in the primitive basis\n";
    os << indent5 << "// with m as the slowest index\n";

    QAMSet allam = vrr_algo_.GetAllAM();

    for(const auto & qam : allam)
    {
        // add +1 fromm required m values to account for 0
        os << indent5 << "// AM = (" << qam[0] << " " << qam[1] << " | " << qam[2] << " " << qam[3] << " )\n";
        os << indent5 << WriterInfo::DoubleType() << " " << WriterInfo::PrimVarName(qam)
           << "[" << (vrr_algo_.GetVRR_MReq(qam)+1) << " * " 
           << NCART(qam) << "] SIMINT_ALIGN_ARRAY_DBL;\n";
    }

    os << "\n\n";
}


void VRR_Writer::DeclarePrimPointers(std::ostream & os) const
{
    QAMSet allam = vrr_algo_.GetAllAM();

    for(const auto & qam : allam)
    {
        if(WriterInfo::IsContArray(qam))
            os << indent4 << "double * restrict " << WriterInfo::PrimPtrName(qam)
               << " = " << WriterInfo::ArrVarName(qam) << " + abcd * " << NCART(qam) << ";\n";
    }

    os << "\n\n";
}


void VRR_Writer::WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepList & vs, const std::string & num_n) const
{
    os << indent5 << "// Forming " << WriterInfo::PrimVarName(qam) << "[" << num_n << " * " << NCART(qam) << "];\n";

    os << indent5 << "for(n = 0; n < " << num_n << "; ++n)  // loop over orders of auxiliary function\n";
    os << indent5 << "{\n";

    os << "\n";

    // TODO - FMA
    // iterate over the requirements
    // should be in order since it's a set
    for(const auto & it : vs)
    {
        // Get the stepping
        XYZStep step = it.xyz;

        std::stringstream primname, srcname[8];
        primname << WriterInfo::PrimVarName(qam) << "[n * " << NCART(qam) << " + " << it.target.idx() << "]";
        if(it.src[0])
            srcname[0] << WriterInfo::PrimVarName(it.src[0].amlist()) << "[n * " << it.src[0].ncart() << " + " << it.src[0].idx() << "]";
        if(it.src[1])
            srcname[1] << WriterInfo::PrimVarName(it.src[1].amlist()) << "[(n+1) * " << it.src[1].ncart() << " + " << it.src[1].idx() << "]";
        if(it.src[2])
            srcname[2] << WriterInfo::PrimVarName(it.src[2].amlist()) << "[n * " << it.src[2].ncart() << " + " << it.src[2].idx() << "]";
        if(it.src[3])
            srcname[3] << WriterInfo::PrimVarName(it.src[3].amlist()) << "[(n+1) * " << it.src[3].ncart() << " + " << it.src[3].idx() << "]";
        if(it.src[4])
            srcname[4] << WriterInfo::PrimVarName(it.src[4].amlist()) << "[n * " << it.src[4].ncart() << " + " << it.src[4].idx() << "]";
        if(it.src[5])
            srcname[5] << WriterInfo::PrimVarName(it.src[5].amlist()) << "[(n+1) * " << it.src[5].ncart() << " + " << it.src[5].idx() << "]";
        if(it.src[6])
            srcname[6] << WriterInfo::PrimVarName(it.src[6].amlist()) << "[(n+1) * " << it.src[6].ncart() << " + " << it.src[6].idx() << "]";
        if(it.src[7])
            srcname[7] << WriterInfo::PrimVarName(it.src[7].amlist()) << "[(n+1) * " << it.src[7].ncart() << " + " << it.src[7].idx() << "]";

        os << indent6 << "//" << it.target <<  " : STEP: " << step << "\n";

        if(it.type == DoubletType::BRA)
        {
            os << indent6 << primname.str() << " = P_PA_" << step << " * " << srcname[0].str();

            if(it.src[1])
                os << " - aop_PQ_" << step << " * " << srcname[1].str(); 
            if(it.src[2] && it.src[3])
                os << " + vrr_const_" << it.ijkl[0] << "_over_2p * (" << srcname[2].str() << " - a_over_p * " << srcname[3].str() << ")"; 
            if(it.src[4] && it.src[5])
                os << " + vrr_const_" << it.ijkl[1] << "_over_2p * (" << srcname[4].str() << " - a_over_p * " << srcname[5].str() << ")"; 
            if(it.src[6])
                os << " + vrr_const_" << it.ijkl[2] << "_over_2pq * " << srcname[6].str();
            if(it.src[7])
                os << " + vrr_const_" << it.ijkl[3] << "_over_2pq * " << srcname[7].str();
            os << ";\n";

        }
        else
        {
            os << indent6 << primname.str() << " = Q_PA_" << step << " * " << srcname[0].str();

            // is this really supposed to be + aoq_PQ and not - ? That's how it's writtein in MEST
            if(it.src[1])
                os << " + aoq_PQ_" << step << " * " << srcname[1].str(); 
            if(it.src[2] && it.src[3])
                os << " + vrr_const_" << it.ijkl[2] << "_over_2q * (" << srcname[2].str() << " - a_over_q * " << srcname[3].str() << ")"; 
            if(it.src[4] && it.src[5])
                os << " + vrr_const_" << it.ijkl[3] << "_over_2q * (" << srcname[4].str() << " - a_over_q * " << srcname[5].str() << ")"; 
            if(it.src[6])
                os << " + vrr_const_" << it.ijkl[0] << "_over_2pq * " << srcname[6].str();
            if(it.src[7])
                os << " + vrr_const_" << it.ijkl[1] << "_one_over_2pq * " << srcname[7].str();
            os << ";\n";
        }
        
        os << "\n";

    }

    os << indent5 << "}\n";
}



void VRR_Writer::WriteVRRInline_(std::ostream & os) const
{
    if(WriterInfo::HasVRR())
    {
        os << "\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << indent5 << "// Primitive integrals: Vertical recurrance\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << "\n";

        // constants
        for(const auto & it : vrr_algo_.GetAllInt_2p())
        {
            if(it == 1)
                os << indent5 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2p = one_over_2p;\n"; 
            else
                os << indent5 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2p = " << WriterInfo::IntConstant(it) << " * one_over_2p;\n"; 
        }

        for(const auto & it : vrr_algo_.GetAllInt_2q())
        {
            if(it == 1)
                os << indent5 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2q = one_over_2q;\n"; 
            else
                os << indent5 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2q = " << WriterInfo::IntConstant(it) << " * one_over_2q;\n"; 
        }

        for(const auto & it : vrr_algo_.GetAllInt_2pq())
        {
            if(it == 1)
                os << indent5 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2pq = one_over_2pq;\n"; 
            else
                os << indent5 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2pq = " << WriterInfo::IntConstant(it) << " * one_over_2pq;\n"; 
        }


        // iterate over increasing am
        for(const auto & qam : vrr_algo_.GetAllAM())
        {
            // Write out the steps
            if(qam != QAM{0,0,0,0})
                WriteVRRSteps_(os, qam, vrr_algo_.GetVRR_Steps(qam), 
                               std::to_string(vrr_algo_.GetVRR_MReq(qam)+1));

            os << "\n\n";
        }

        os << "\n\n";
    }
}



void VRR_Writer::WriteVRRFile(std::ostream & os) const
{
    QAM am = WriterInfo::FinalAM();

    os << "//////////////////////////////////////////////\n";
    os << "// VRR: ( " << amchar[am[0]] << " " << amchar[am[1]] << " | " << amchar[am[2]] << " " << amchar[am[3]] << " )\n";
    os << "//////////////////////////////////////////////\n";

    // we only do the final
    //if(!WriterInfo::Scalar())
    //    os << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(num_n)\n";
    os << "void VRR_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

    // final target
    os << indent3 << WriterInfo::DoubleType() << " * const restrict " << WriterInfo::PrimVarName(am) << ",\n";

    for(const auto & it : vrr_algo_.Get_AMReq(am))
        os << indent3 << WriterInfo::ConstDoubleType() << " * const restrict " << WriterInfo::PrimVarName(it) << ",\n";
    
    for(const auto & it : vrr_algo_.GetVarReq(am))
        os << indent3 << WriterInfo::ConstDoubleType() << " " << it << ",\n";

    os << indent3 << "const int num_n)\n";

    os << "{\n";

    os << "    int n = 0;\n";

    for(const auto & it : vrr_algo_.GetIntReq_2p(am))
        os << indent1 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2p = " << WriterInfo::DoubleSet1(std::to_string(it)) << " * one_over_2p;\n"; 

    for(const auto & it : vrr_algo_.GetIntReq_2q(am))
        os << indent1 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2q = " << WriterInfo::DoubleSet1(std::to_string(it)) << " * one_over_2p;\n"; 

    for(const auto & it : vrr_algo_.GetIntReq_2pq(am))
        os << indent1 << WriterInfo::ConstDoubleType() << " vrr_const_" << it << "_over_2pq = " << WriterInfo::DoubleSet1(std::to_string(it)) << " * one_over_2pq;\n"; 

    // Write out the steps
    WriteVRRSteps_(os, am, vrr_algo_.GetVRR_Steps(am), "num_n");

    os << "\n";
    os << "}\n";
    os << "\n\n";
}




void VRR_Writer::WriteVRRHeaderFile(std::ostream & os) const
{
}



void VRR_Writer::WriteVRRExternal_(std::ostream & os) const
{
/*
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
        os << indent6 << "VRR_" << amchar[am] << "(" << (WriterInfo::L()-am+1) << ",\n";
        os << indent7 << "P_PA_x, P_PA_y, P_PA_z, aop_PQ_x, aop_PQ_y, aop_PQ_z,\n";
        os << indent7 << "a_over_p,";
        if(am > 1)
            os << "one_over_2p, ";
        os << "\n"; 

        os << indent7 << WriterInfo::PrimVarName(qam) << ",\n";
        os << indent7 << WriterInfo::PrimVarName(qam1);
        if(am > 1)
        {
            os << ",\n";
            os << indent7 << WriterInfo::PrimVarName(qam2);
        }

        os << ");\n";

    }
*/
}



void VRR_Writer::WriteVRR(std::ostream & os) const
{
    if(WriterInfo::GetOption(OPTION_INLINEVRR) > 0)
        WriteVRRInline_(os);
    else
        WriteVRRExternal_(os);
}

