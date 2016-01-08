#include "generator/ET_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/ET_Algorithm_Base.hpp"
#include "generator/Helpers.hpp"
#include "generator/Ncart.hpp"
#include "generator/Naming.hpp"


ET_Writer::ET_Writer(const ET_Algorithm_Base & et_algo) 
    : et_algo_(et_algo)
{ 
}

bool ET_Writer::HasET(void) const
{
    return et_algo_.HasET();
}

bool ET_Writer::HasBraET(void) const
{
    return et_algo_.HasBraET();
}


bool ET_Writer::HasKetET(void) const
{
    return et_algo_.HasKetET();
}


void ET_Writer::AddConstants(ERIGeneratorInfo & info) const
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
            os << indent5 << WriterInfo::DoubleType() << " " << PrimVarName(it) << "[" << NCART(it[0], it[2]) << "] SIMINT_ALIGN_ARRAY_DBL;\n";

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
                os << indent4  << "double * restrict " << PrimPtrName(it)
                   << " = " << ArrVarName(it) << " + abcd * " << NCART(it[0], it[2]) << ";\n";
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
        std::string braket = (et_algo_.GetDirection() == DoubletType::KET) ? "KET" : "BRA";

        os << "\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << indent5 << "// Primitive integrals: Electron transfer\n";
        os << indent5 << "//////////////////////////////////////////////\n";
        os << "\n";

        // loop over the steps in order
        for(const auto & am : et_algo_.GetAMOrder())
        {
            os << indent5 << "ET_" << braket << "_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n"; 
            os << indent6 << PrimVarName(am) << ",\n";

            for(const auto & it : et_algo_.GetAMReq(am))
                os << indent6 << PrimVarName(it) << ",\n";

            if(et_algo_.GetDirection() == DoubletType::KET)
            {
                os << indent5 << "etfac_k, p_over_q";
                if( (am[0] + am[1] + am[2] + am[3]) > 1 )
                    os << ", one_over_2q";
                os << ");\n\n";
            }
            else
            {
                os << indent5 << "etfac_b, q_over_p";
                if( (am[0] + am[1] + am[2] + am[3]) > 1 )
                    os << ", one_over_2p";
                os << ");\n\n";
            }
        }
    }
}


void ET_Writer::WriteETFile(std::ostream & os, std::ostream & osh) const
{
    QAM am = WriterInfo::FinalAM();
    std::string braket = (et_algo_.GetDirection() == DoubletType::KET) ? "KET" : "BRA";

    os << "//////////////////////////////////////////////\n";
    os << "// ET: ( " << amchar[am[0]] << " " << amchar[am[1]] << " | " << amchar[am[2]] << " " << amchar[am[3]] << " )\n";
    os << "//////////////////////////////////////////////\n";

    std::stringstream prototype;
    prototype << "void ET_" << braket << "_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";
    prototype << indent3 << WriterInfo::DoubleType() << " * const restrict " << PrimVarName(am) << ",\n";

    for(const auto & it : et_algo_.GetAMReq(am))
        prototype << indent3 << WriterInfo::ConstDoubleType() << " * const restrict " << PrimVarName(it) << ",\n";

    if(et_algo_.GetDirection() == DoubletType::KET)
    {
        prototype << indent3 << WriterInfo::DoubleType() << " const * const restrict etfac_k, "
           << WriterInfo::ConstDoubleType() << " p_over_q";
        if( (am[0] + am[1] + am[2] + am[3]) > 1 )
            prototype << ", " << WriterInfo::ConstDoubleType() << " one_over_2q";
        prototype << ")\n";
    }
    else
    {
        prototype << indent3 << WriterInfo::DoubleType() << " const * const restrict etfac_b, "
           << WriterInfo::ConstDoubleType() << " q_over_p";
        if( (am[0] + am[1] + am[2] + am[3]) > 1 )
            prototype << ", " << WriterInfo::ConstDoubleType() << " one_over_2p";
        prototype << ")";
    }


    os << prototype.str() << "\n";
    os << "{\n";
    os << indent2 << "// Precompute (integer) * 1/2{pq}\n";

    for(const auto & it : et_algo_.GetIntReq(am))
    {
        DoubletType dir = et_algo_.GetDirection();
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

    osh << prototype.str() << ";\n\n";
}


std::string ET_Writer::ETStepVar_(const Quartet & q) const
{
    std::stringstream ss; 
    ss << PrimVarName(q.amlist()) << "[" << q.idx() << "]";
    return ss.str();
}



std::string ET_Writer::ETStepString_(const ETStep & et) const
{
    int stepidx = XYZStepToIdx(et.xyz);

    DoubletType direction = et_algo_.GetDirection();

    // for output
    char kb = (direction == DoubletType::KET) ? 'k' : 'b';
    char pq = (direction == DoubletType::KET) ? 'q' : 'p';

    std::stringstream etconst_i, etconst_k;
    std::stringstream ss;

    etconst_i << "et_const_" << pq << "_" << et.ik[0];
    etconst_k << "et_const_" << pq << "_" << et.ik[1];

    std::stringstream etfac, aoverb;
    etfac << "etfac_" << kb << "[" << stepidx << "]";

    std::string poverq = (direction == DoubletType::KET) ? "-p_over_q" : "-q_over_p";
    

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

