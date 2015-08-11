#include "generator/ET_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/ET_Algorithm_Base.hpp"
#include "generator/Helpers.hpp"
#include "generator/Ncart.hpp"


ET_Writer::ET_Writer(const ET_Algorithm_Base & et_algo) 
    : et_algo_(et_algo)
{ 
}


bool ET_Writer::HasBraET(void) const
{
    return et_algo_.HasBraET();
}


bool ET_Writer::HasKetET(void) const
{
    return et_algo_.HasKetET();
}


void ET_Writer::AddConstants(void) const
{
    if(WriterInfo::GetOption(OPTION_INLINEET) > 0)
    {
        for(const auto & it : et_algo_.GetAllInt())
            WriterInfo::AddIntConstant(it);
    }
}



void ET_Writer::DeclarePrimArrays(std::ostream & os) const
{
    QAMSet allam = et_algo_.GetAllAM();
 
    if(allam.size())
    {
        os << indent5 << "// Holds temporary integrals for electron transfer\n";

        for(const auto & it : allam)
            os << indent5 << WriterInfo::DoubleType() << " " << WriterInfo::PrimVarName(it) << "[" << NCART(it[0], it[2]) << "] SIMINT_ALIGN_ARRAY_DBL;\n";

        os << "\n\n";

    }
}


void ET_Writer::DeclarePrimPointers(std::ostream & os) const
{
    QAMSet allam = et_algo_.GetAllAM();
 
    if(allam.size())
    {
        for(const auto & it : allam)
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
    if(WriterInfo::GetOption(OPTION_INLINEET) > 0)
        WriteETInline_(os);
    else
        WriteETExternal_(os);
}



void ET_Writer::WriteETInline_(std::ostream & os) const
{
    if(WriterInfo::HasET())
    {
        os << "\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << indent5 << "// Primitive integrals: Electron transfer\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << "\n";

        os << indent5 << "// Precompute (integer) * 1/2{pq}\n";
        for(const auto & it : et_algo_.GetAllInt_p())
        {
            if(it != 1)
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_p_" << it << " = " << WriterInfo::DoubleSet1(std::to_string(it)) << " * one_over_2p;\n";
            else
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_p_" << "1 = one_over_2p;\n";
        }

        for(const auto & it : et_algo_.GetAllInt_q())
        {
            if(it != 1)
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_q_" << it << " = " << WriterInfo::DoubleSet1(std::to_string(it)) << " * one_over_2q;\n";
            else
                os << indent5 << WriterInfo::ConstDoubleType() << " et_const_q_" << "1 = one_over_2q;\n";
        }

        os << "\n\n";

        // loop over the stpes in order
        for(const auto & am : et_algo_.GetAMOrder())
        {
            ETStepList etsl = et_algo_.GetSteps(am);

            for(const auto & it : etsl)
            {
                os << indent5 << "// " << it << "\n";
                os << ETStepString_(it);
                os << "\n";
            }
        }
    }
}


void ET_Writer::WriteETExternal_(std::ostream & os) const
{
    if(WriterInfo::HasET())
    {
        os << "\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << indent5 << "// Primitive integrals: Electron transfer\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << "\n";

        // loop over the steps in order
        for(const auto & am : et_algo_.GetAMOrder())
        {
            os << indent5 << "ET_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n"; 
            os << indent6 << WriterInfo::PrimVarName(am) << ",\n";

            for(const auto & it : et_algo_.GetAMReq(am))
                os << indent6 << WriterInfo::PrimVarName(it) << ",\n";

            if(et_algo_.GetDirection(am) == DoubletType::KET)
               os << indent5 << "etfac_k, one_over_2q, p_over_q);\n\n";
            else
               os << indent5 << "etfac_b, one_over_2p, q_over_p);\n\n";
        }
    }
}


void ET_Writer::WriteETFile(std::ostream & os, std::ostream & osh) const
{
    QAM am = WriterInfo::FinalAM();

    os << "//////////////////////////////////////////////\n";
    os << "// ET: ( " << amchar[am[0]] << " " << amchar[am[1]] << " | " << amchar[am[2]] << " " << amchar[am[3]] << " )\n";
    os << "//////////////////////////////////////////////\n";


    os << "void ET_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

    os << indent3 << WriterInfo::DoubleType() << " * const restrict " << WriterInfo::PrimVarName(am) << ",\n";

    for(const auto & it : et_algo_.GetAMReq(am))
        os << indent3 << WriterInfo::ConstDoubleType() << " * const restrict " << WriterInfo::PrimVarName(it) << ",\n";

    if(et_algo_.GetDirection(am) == DoubletType::KET)
        os << indent3 << WriterInfo::DoubleType() << " const * const restrict etfac_k, "
           << WriterInfo::ConstDoubleType() << " one_over_2q, "
           << WriterInfo::ConstDoubleType() << " p_over_q)\n";
    else
        os << indent3 << WriterInfo::DoubleType() << " const * const restrict etfac_b, "
           << WriterInfo::ConstDoubleType() << " one_over_2p, "
           << WriterInfo::ConstDoubleType() << " q_over_p)\n";

    os << "{\n";
    os << indent2 << "// Precompute (integer) * 1/2{pq}\n";

    for(const auto & it : et_algo_.GetIntReq(am))
    {
        DoubletType dir = et_algo_.GetDirection(am);
        char pq = (dir == DoubletType::KET) ? 'q' : 'p';

        if(it != 1)
            os << indent2 << WriterInfo::ConstDoubleType() << " et_const_" << pq << "_" << it << " = " << WriterInfo::DoubleSet1(std::to_string(it)) << " * one_over_2" << pq << ";\n";
        else
            os << indent2 << WriterInfo::ConstDoubleType() << " et_const_" << pq << "_" << "1 = one_over_2" << pq << ";\n";
    }

    os << "\n\n";

    ETStepList etsl = et_algo_.GetSteps(am);

    for(const auto & it : etsl)
    {
        os << indent2 << "// " << it << "\n";
        os << ETStepString_(it);
        os << "\n";
    }

    os << "}\n\n\n";


    // add to header
    osh << "void ET_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

    osh << indent3 << WriterInfo::DoubleType() << " * const restrict " << WriterInfo::PrimVarName(am) << ",\n";

    for(const auto & it : et_algo_.GetAMReq(am))
        osh << indent3 << WriterInfo::ConstDoubleType() << " * const restrict " << WriterInfo::PrimVarName(it) << ",\n";

    if(et_algo_.GetDirection(am) == DoubletType::KET)
        osh << indent3 << WriterInfo::DoubleType() << " const * const restrict etfac_k, "
            << WriterInfo::ConstDoubleType() << " one_over_2q, "
            << WriterInfo::ConstDoubleType() << " p_over_q);\n";
    else
        osh << indent3 << WriterInfo::DoubleType() << " const * const restrict etfac_b, "
            << WriterInfo::ConstDoubleType() << " one_over_2p, "
            << WriterInfo::ConstDoubleType() << " q_over_p);\n";

    osh << "\n";
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

    // for output
    char kb = (et.direction == DoubletType::KET) ? 'k' : 'b';

    std::stringstream etconst_i, etconst_k;
    std::stringstream ss;

    etconst_i << "et_const_" << et.pq << "_" << et.ik[0];
    etconst_k << "et_const_" << et.pq << "_" << et.ik[1];

    std::stringstream etfac, aoverb;
    etfac << "etfac_" << kb << "[" << stepidx << "]";

    std::string poverq = (et.direction == DoubletType::KET) ? "-p_over_q" : "-q_over_p";
    

    if(WriterInfo::HasFMA())
    {
        ss << indent5 << ETStepVar_(et.target) << " = " << etfac.str() << " * " << ETStepVar_(et.src[0]) << ";\n";
        if(et.src[1])
            ss << indent5 << ETStepVar_(et.target) << " = " << WriterInfo::FMAdd(etconst_i.str(), ETStepVar_(et.src[1]), ETStepVar_(et.target)) << ";\n";
        if(et.src[2])
            ss << indent5 << ETStepVar_(et.target) << " = " << WriterInfo::FMAdd(etconst_k.str(), ETStepVar_(et.src[2]), ETStepVar_(et.target)) << ";\n";
        if(et.src[3])
            ss << indent5 << ETStepVar_(et.target) << " = " << WriterInfo::FMAdd(poverq, ETStepVar_(et.src[3]), ETStepVar_(et.target)) << ";\n";

    }
    else
    {
        ss << indent5 << ETStepVar_(et.target);

        ss << " = ";

        ss << etfac.str() << " * " << ETStepVar_(et.src[0]);

        if(et.src[1])
            ss << " + " << etconst_i.str() << " * " << ETStepVar_(et.src[1]);
        if(et.src[2])
            ss << " + " << etconst_k.str() << " * " << ETStepVar_(et.src[2]);
        if(et.src[3])
            ss << " " << poverq << " * " << ETStepVar_(et.src[3]);
        ss << ";\n";
    }

    return ss.str();
}

