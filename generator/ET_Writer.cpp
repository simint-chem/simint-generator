#include "generator/ET_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/ET_Algorithm_Base.hpp"
#include "generator/Helpers.hpp"


ET_Writer::ET_Writer(const ET_Algorithm_Base & et_algo) 
{ 
    etsl_ = et_algo.ETSteps();

    // see what we need for arrays
    for(const auto & it : etsl_)
    {
        if(it.src1)
            etint_.insert(it.src1.amlist());
        if(it.src2)
            etint_.insert(it.src2.amlist());
        if(it.src3)
            etint_.insert(it.src3.amlist());
        if(it.src3)
            etint_.insert(it.src4.amlist());
        if(it.target)
            etint_.insert(it.target.amlist());
    }

    // determine the integers that get multiplied by 1/2q
    for(const auto & et : etsl_)
    {
        int stepidx = XYZStepToIdx(et.xyz);

        int ival = et.target.bra.left.ijk[stepidx];
        int kval = (et.target.ket.left.ijk[stepidx]-1);

        // ok to add duplicates. They are stored as a map
        if(ival > 0)
            et_i_.insert(ival);
        if(kval > 0)
            et_i_.insert(kval);
    }
}



void ET_Writer::WriteIncludes(std::ostream & os) const
{
}



void ET_Writer::AddConstants(void) const
{
    for(const auto & it : et_i_)
        WriterInfo::AddIntConstant(it);
}



void ET_Writer::DeclarePrimArrays(std::ostream & os) const
{
    if(etint_.size())
    {
        os << indent6 << "// Holds temporary integrals for electron transfer\n";

        for(const auto & it : etint_)
        {
            // only if these aren't from vrr
            if(it[1] > 0 || it[2] > 0 || it[3] > 0)
                os << indent6 << WriterInfo::DoubleType() << " " << WriterInfo::PrimVarName(it) << "[" << NCART(it[0]) * NCART(it[2]) << "] SIMINT_ALIGN_ARRAY;\n";
        } 

        os << "\n\n";

    }
}


void ET_Writer::DeclarePrimPointers(std::ostream & os) const
{
    if(etint_.size())
    {
        for(const auto & it : etint_)
        {
            if(WriterInfo::IsContArray(it))
                os << indent4  << "double * const restrict " << WriterInfo::PrimPtrName(it)
                   << " = " << WriterInfo::ArrVarName(it) << " + abcd * " << NCART(it[0]) * NCART(it[2]) << ";\n";
        }

        os << "\n\n";

    }
}


void ET_Writer::WriteETInline(std::ostream & os) const
{
    os << "\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << indent6 << "// Primitive integrals: Electron transfer\n";
    os << indent6 << "//////////////////////////////////////////////\n";
    os << "\n";
    os << indent6 << "// Precompute (integer) * 1/2q\n";
    for(const auto & it : et_i_)
        os << indent6 << WriterInfo::ConstDoubleType() << " et_const_" << it << " = " << WriterInfo::IntConstant(it) << " * one_over_2q;\n";

    os << "\n";

    if(etsl_.size() == 0)
        os << indent6 << "//...nothing to do...\n";
    else
    {
        for(const auto & it : etsl_)
        {
            os << indent6 << "// " << it << "\n";
            os << ETStepString(it);
            os << "\n";
        }
    }

    // add to needed contracted integrals
    for(const auto & it : etint_)
    {
        if(WriterInfo::IsContArray(it))
        {
            int ncart = NCART(it[0])*NCART(it[2]);

            os << "\n";
            os << indent6 << "// Accumulating in contracted workspace\n";

            os << indent6 << "for(n = 0; n < " << ncart << "; n++)\n";

            if(WriterInfo::Intrinsics())
            {
                os << indent6 << "{\n";

                os << indent7 << WriterInfo::UnionType() << " vec = (" << WriterInfo::UnionType() << ")" << WriterInfo::PrimVarName(it) << "[n];\n";    
                os << indent7 << WriterInfo::PrimPtrName(it) << "[n] += vec.d[0]";

                for(int i = 1; i < WriterInfo::SimdLen(); i++)
                    os << " + vec.d[" << i << "]";
                    //os << " + vec[" << i << "]";
                os << ";\n";
                os << indent6 << "}\n";
            }
            else
                os << indent7 << WriterInfo::PrimPtrName(it) << "[n] += " << WriterInfo::PrimVarName(it) << "[n];\n";
        }
    }
}



std::string ET_Writer::ETStepVar(const Quartet & q)
{
    std::stringstream ss; 
    ss << WriterInfo::PrimVarName(q.amlist()) << "[" << q.idx() << "]";
    return ss.str();
}



std::string ET_Writer::ETStepString(const ETStep & et)
{
    int stepidx = XYZStepToIdx(et.xyz);
    int ival = et.target.bra.left.ijk[stepidx];
    int kval = (et.target.ket.left.ijk[stepidx]-1);

    // for output
    std::stringstream ss;


    std::stringstream etconst_i, etconst_k;
    etconst_i << "et_const_" << ival;
    etconst_k << "et_const_" << kval;

    std::stringstream etfac;
    etfac << "etfac[" << stepidx << "]";

    if(WriterInfo::HasFMA())
    {
        ss << indent6 << ETStepVar(et.target) << " = " << etfac.str() << " * " << ETStepVar(et.src1) << ";\n";
        if(et.src2.bra.left && et.src2.ket.left)
            ss << indent6 << ETStepVar(et.target) << " = " << WriterInfo::FMAdd(etconst_i.str(), ETStepVar(et.src2), ETStepVar(et.target)) << ";\n";
        if(et.src3.bra.left && et.src3.ket.left)
            ss << indent6 << ETStepVar(et.target) << " = " << WriterInfo::FMAdd(etconst_k.str(), ETStepVar(et.src3), ETStepVar(et.target)) << ";\n";
        if(et.src4.bra.left && et.src4.ket.left)
            ss << indent6 << ETStepVar(et.target) << " = " << WriterInfo::FMAdd("-p_over_q", ETStepVar(et.src4), ETStepVar(et.target)) << ";\n";

    }
    else
    {
        ss << indent6 << ETStepVar(et.target);


        ss << " = ";

        ss << etfac.str() << " * " << ETStepVar(et.src1);

        if(et.src2.bra.left && et.src2.ket.left)
            ss << " + " << etconst_i.str() << " * " << ETStepVar(et.src2);
        if(et.src3.bra.left && et.src3.ket.left)
            ss << " + " << etconst_k.str() << " * " << ETStepVar(et.src3);
        if(et.src4.bra.left && et.src4.ket.left)
            ss << " - p_over_q * " << ETStepVar(et.src4);
        ss << ";\n";
    }

    return ss.str();
}

