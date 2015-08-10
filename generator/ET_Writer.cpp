#include "generator/ET_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/ET_Algorithm_Base.hpp"
#include "generator/Helpers.hpp"
#include "generator/Ncart.hpp"


ET_Writer::ET_Writer(const ET_Algorithm_Base & et_algo) 
{ 
    etsl_ = et_algo.ETSteps();

    // see what we need for arrays
    for(const auto & it : etsl_)
    {
        if(it.src[0])
            etint_.insert(it.src[0].amlist());
        if(it.src[1])
            etint_.insert(it.src[1].amlist());
        if(it.src[2])
            etint_.insert(it.src[2].amlist());
        if(it.src[2])
            etint_.insert(it.src[3].amlist());
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



void ET_Writer::AddConstants(void) const
{
    for(const auto & it : et_i_)
        WriterInfo::AddIntConstant(it);
}



void ET_Writer::DeclarePrimArrays(std::ostream & os) const
{
    if(etint_.size())
    {
        os << indent5 << "// Holds temporary integrals for electron transfer\n";

        for(const auto & it : etint_)
        {
            // only if these aren't from vrr
            if(it[1] > 0 || it[2] > 0 || it[3] > 0)
                os << indent5 << WriterInfo::DoubleType() << " " << WriterInfo::PrimVarName(it) << "[" << NCART(it[0], it[2]) << "] SIMINT_ALIGN_ARRAY_DBL;\n";
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
                os << indent4  << "double * restrict " << WriterInfo::PrimPtrName(it)
                   << " = " << WriterInfo::ArrVarName(it) << " + abcd * " << NCART(it[0], it[2]) << ";\n";
        }

        os << "\n\n";

    }
}


void ET_Writer::WriteET(std::ostream & os) const
{
    if(WriterInfo::GetOption(OPTION_INLINEHRR) > 0)
        WriteETInline_(os);
    else
        WriteETExternal_(os);
}



void ET_Writer::WriteETInline_(std::ostream & os) const
{
    if(et_i_.size())
    {
        os << "\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << indent5 << "// Primitive integrals: Electron transfer\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << "\n";

        os << indent5 << "// Precompute (integer) * 1/2q\n";
        for(const auto & it : et_i_)
        {
            if(it != 1)
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_" << it << " = " << WriterInfo::IntConstant(it) << " * one_over_2q;\n";
            else
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_1 = one_over_2q;\n";
        }

        for(const auto & it : etsl_)
        {
            os << indent5 << "// " << it << "\n";
            os << ETStepString_(it);
            os << "\n";
        }
    }
}


void ET_Writer::WriteETExternal_(std::ostream & os) const
{
    if(et_i_.size())
    {
        os << "\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << indent5 << "// Primitive integrals: Electron transfer\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << "\n";

        os << indent5 << "// Precompute (integer) * 1/2q\n";
        for(const auto & it : et_i_)
        {
            if(it != 1)
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_" << it << " = " << WriterInfo::IntConstant(it) << " * one_over_2q;\n";
            else
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_1 = one_over_2q;\n";
        }

        for(const auto & it : etsl_)
        {
            os << indent5 << "// " << it << "\n";
            os << ETStepString_(it);
            os << "\n";
        }
    }
}


void ET_Writer::WriteETFile(std::ostream & os, std::ostream & osh) const
{
}


std::string ET_Writer::ETStepVar_(const Quartet & q)
{
    std::stringstream ss; 
    ss << WriterInfo::PrimVarName(q.amlist()) << "[" << q.idx() << "]";
    return ss.str();
}



std::string ET_Writer::ETStepString_(const ETStep & et)
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
        ss << indent5 << ETStepVar_(et.target) << " = " << etfac.str() << " * " << ETStepVar_(et.src[0]) << ";\n";
        if(et.src[1].bra.left && et.src[1].ket.left)
            ss << indent5 << ETStepVar_(et.target) << " = " << WriterInfo::FMAdd(etconst_i.str(), ETStepVar_(et.src[1]), ETStepVar_(et.target)) << ";\n";
        if(et.src[2].bra.left && et.src[2].ket.left)
            ss << indent5 << ETStepVar_(et.target) << " = " << WriterInfo::FMAdd(etconst_k.str(), ETStepVar_(et.src[2]), ETStepVar_(et.target)) << ";\n";
        if(et.src[3].bra.left && et.src[3].ket.left)
            ss << indent5 << ETStepVar_(et.target) << " = " << WriterInfo::FMAdd("-p_over_q", ETStepVar_(et.src[3]), ETStepVar_(et.target)) << ";\n";

    }
    else
    {
        ss << indent5 << ETStepVar_(et.target);


        ss << " = ";

        ss << etfac.str() << " * " << ETStepVar_(et.src[0]);

        if(et.src[1].bra.left && et.src[1].ket.left)
            ss << " + " << etconst_i.str() << " * " << ETStepVar_(et.src[1]);
        if(et.src[2].bra.left && et.src[2].ket.left)
            ss << " + " << etconst_k.str() << " * " << ETStepVar_(et.src[2]);
        if(et.src[3].bra.left && et.src[3].ket.left)
            ss << " - p_over_q * " << ETStepVar_(et.src[3]);
        ss << ";\n";
    }

    return ss.str();
}

