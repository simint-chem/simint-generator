#include "generator/VRR_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/VRR_Algorithm_Base.hpp"
#include "generator/Helpers.hpp"
#include "generator/Ncart.hpp"


VRR_Writer::VRR_Writer(const VRR_Algorithm_Base & vrr_algo)
{ 
    vrrmap_ = vrr_algo.GetVRRMap();
    maxFm_ = vrr_algo.GetMaxFm();
    vrrmreq_ = vrr_algo.GetVRRMReq();

    // determine the constant integers
    // needed

    for(const auto & it : vrrmap_)
    {
        const VRRStepList & it2 = it.second;


        for(const auto & its : it2)
        {
            int istep = XYZStepToIdx(its.xyz);
            if(its.type == DoubletType::BRA)
            {
                if(its.src[2] || its.src[3])
                    vrr_bra_i_.insert(its.target.bra.left.ijk[istep]-1);
                if(its.src[4] || its.src[5])
                    vrr_bra_j_.insert(its.target.bra.right.ijk[istep]);
                if(its.src[6])
                    vrr_bra_k_.insert(its.target.ket.left.ijk[istep]);
                if(its.src[7])
                    vrr_bra_l_.insert(its.target.ket.right.ijk[istep]);
            }
            else
            {
                if(its.src[2] || its.src[3])
                    vrr_ket_k_.insert(its.target.bra.left.ijk[istep]-1);
                if(its.src[4] || its.src[5])
                    vrr_ket_l_.insert(its.target.bra.right.ijk[istep]);
                if(its.src[6])
                    vrr_ket_i_.insert(its.target.ket.left.ijk[istep]);
                if(its.src[7])
                    vrr_ket_j_.insert(its.target.ket.right.ijk[istep]);
            }
        }
    }

    // these factors get multiplied by 1 / 2(p+q)
    vrr_2pq_.insert(vrr_bra_k_.begin(), vrr_bra_k_.end());
    vrr_2pq_.insert(vrr_bra_l_.begin(), vrr_bra_l_.end());
    vrr_2pq_.insert(vrr_ket_i_.begin(), vrr_ket_i_.end());
    vrr_2pq_.insert(vrr_ket_j_.begin(), vrr_ket_j_.end());
}


void VRR_Writer::WriteIncludes(std::ostream & os) const
{
    if(WriterInfo::GetOption(OPTION_INLINEVRR) == 0)
        os << "#include \"eri/vrr.gen/vrr.h\"\n";
}



void VRR_Writer::AddConstants(void) const
{
    for(const auto & it : vrr_bra_i_)
        WriterInfo::AddIntConstant(it);
    for(const auto & it : vrr_bra_j_)
        WriterInfo::AddIntConstant(it);
    for(const auto & it : vrr_bra_k_)
        WriterInfo::AddIntConstant(it);
    for(const auto & it : vrr_bra_l_)
        WriterInfo::AddIntConstant(it);

    for(const auto & it : vrr_ket_i_)
        WriterInfo::AddIntConstant(it);
    for(const auto & it : vrr_ket_j_)
        WriterInfo::AddIntConstant(it);
    for(const auto & it : vrr_ket_k_)
        WriterInfo::AddIntConstant(it);
    for(const auto & it : vrr_ket_l_)
        WriterInfo::AddIntConstant(it);
}


void VRR_Writer::DeclarePrimArrays(std::ostream & os) const
{
    os << indent5 << "// Holds the auxiliary integrals ( i 0 | k 0 )^m in the primitive basis\n";
    os << indent5 << "// with m as the slowest index\n";
    for(const auto & it : vrrmap_)
    {
        QAM qam = it.first;

        // add +1 fromm required m values to account for 0
        os << indent5 << "// AM = (" << qam[0] << " " << qam[1] << " | " << qam[2] << " " << qam[3] << " )\n";
        os << indent5 << WriterInfo::DoubleType() << " " << WriterInfo::PrimVarName(qam)
           << "[" << (vrrmreq_.at(qam)+1) << " * " << NCART(qam) << "] SIMINT_ALIGN_ARRAY_DBL;\n";
        //os << "\n";
    }

    os << "\n\n";
}


void VRR_Writer::DeclarePrimPointers(std::ostream & os) const
{
    for(const auto & it : vrrmap_)
    {
        QAM qam = it.first;
        if(WriterInfo::IsContArray(qam))
            os << indent4 << "double * restrict " << WriterInfo::PrimPtrName(qam)
               << " = " << WriterInfo::ArrVarName(qam) << " + abcd * " << NCART(qam) << ";\n";
    }

    os << "\n\n";
}


void VRR_Writer::WriteVRRSteps_(std::ostream & os, QAM qam, const VRRStepList & vs, const std::string & num_n) const
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
        int istep = XYZStepToIdx(step);

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
            int vrr_i = it.target.bra.left.ijk[istep]-1;
            int vrr_j = it.target.bra.right.ijk[istep];
            int vrr_k = it.target.ket.left.ijk[istep];
            int vrr_l = it.target.ket.right.ijk[istep];

            os << indent6 << primname.str() << " = P_PA_" << step << " * " << srcname[0].str();

            if(it.src[1])
                os << " - aop_PQ_" << step << " * " << srcname[1].str(); 
            if(it.src[2] && it.src[3])
                os << " + const_" << vrr_i << " * one_over_2p * (" << srcname[2].str() << " - a_over_p * " << srcname[3].str() << ")"; 
            if(it.src[4] && it.src[5])
                os << " + const_" << vrr_j << " * one_over_2p * (" << srcname[4].str() << " - a_over_p * " << srcname[5].str() << ")"; 
            if(it.src[6])
                os << " + const_" << vrr_k << " * one_over_2pq * " << srcname[6].str();
            if(it.src[7])
                os << " + const_" << vrr_l << " * one_over_2pq * " << srcname[7].str();
            os << ";\n";

        }
        else
        {
            int vrr_i = it.target.bra.left.ijk[istep];
            int vrr_j = it.target.bra.right.ijk[istep];
            int vrr_k = it.target.ket.left.ijk[istep]-1;
            int vrr_l = it.target.ket.right.ijk[istep];

            os << indent6 << primname.str() << " = Q_PA_" << step << " * " << srcname[0].str();

            // is this really supposed to be + aoq_PQ and not - ? That's how it's writtein in MEST
            if(it.src[1])
                os << " + aoq_PQ_" << step << " * " << srcname[1].str(); 
            if(it.src[2] && it.src[3])
                os << " + const_" << vrr_k << " * one_over_2q * (" << srcname[2].str() << " - a_over_q * " << srcname[3].str() << ")"; 
            if(it.src[4] && it.src[5])
                os << " + const_" << vrr_l << " * one_over_2q * (" << srcname[4].str() << " - a_over_q * " << srcname[5].str() << ")"; 
            if(it.src[6])
                os << " + const_" << vrr_i << " * one_over_2pq * " << srcname[6].str();
            if(it.src[7])
                os << " + const_" << vrr_j << " * one_over_2pq * " << srcname[7].str();
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

        // iterate over increasing am
        for(const auto & it : vrrmap_)
        {
            // Write out the steps
            QAM qam = it.first;
            if(qam != QAM{0,0,0,0})
                WriteVRRSteps_(os, qam, it.second, std::to_string(vrrmreq_.at(qam)+1));

            os << "\n\n";
        }

        os << "\n\n";
    }
}



void VRR_Writer::WriteVRRFile(std::ostream & os) const
{
/*
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
        os << "// VRR to obtain " << WriterInfo::PrimVarName(qam) << "\n";

        if(!WriterInfo::Scalar())
            os << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(num_n)\n";

        os << "void VRR_" << amchar[am] << "(const int num_n,\n";
        os << indent3 << "const double P_PA_x, const double P_PA_y, const double P_PA_z,\n";
        os << indent3 << "const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,\n";
        os << indent3 << "const double a_over_p,";
        if(am > 1)
            os << " const double one_over_2p,";
        os << "\n"; 
        os << indent3 << "double * const restrict " << WriterInfo::PrimVarName(qam) << ",\n";
        os << indent3 << "double const * const restrict " << WriterInfo::PrimVarName(qam1);
        if(am > 1)
        {
            os << ",\n";
            os << indent3 << "double const * const restrict " << WriterInfo::PrimVarName(qam2);
        }
        
        os << ")\n";
        os << "{\n";

        os << "    int n = 0;\n";

        // Write out the steps
        WriteVRRSteps_(os, greq, "num_n");

        os << "}\n";
    }
*/
}




void VRR_Writer::WriteVRRHeaderFile(std::ostream & os) const
{
/*
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
        if(!WriterInfo::Scalar())
            os << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(num_n)\n";
        os << "void VRR_" << amchar[am] << "(const int num_n,\n";
        os << indent3 << "const double P_PA_x, const double P_PA_y, const double P_PA_z,\n";
        os << indent3 << "const double aop_PQ_x, const double aop_PQ_y, const double aop_PQ_z,\n";
        os << indent3 << "const double a_over_p,";
        if(am > 1)
            os << " const double one_over_2p,";
        os << "\n"; 
        os << indent3 << "double * const restrict " << WriterInfo::PrimVarName(qam) << ",\n";
        os << indent3 << "double const * const restrict " << WriterInfo::PrimVarName(qam1);
        if(am > 1)
        {
            os << ",\n";
            os << indent3 << "double const * const restrict " << WriterInfo::PrimVarName(qam2);
        }
        
        os << ");\n";
    }

    os << "#endif\n";
*/  
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

