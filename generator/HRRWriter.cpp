#include "generator/HRRWriter.hpp"
#include "generator/WriterBase.hpp"
#include "generator/Helpers.hpp"

HRRWriter::HRRWriter(const HRRBraKetStepList & hrrsteps, const QAMList & finalam) 
          : hrrsteps_(hrrsteps)
{
    // determine top bras/kets
    for(const auto & it : hrrsteps.first)
    {
        if(it.src1.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopbras_[it.src1.am()].insert(it.src1);
        if(it.src2.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopbras_[it.src2.am()].insert(it.src2);
    }   

    for(const auto & it : hrrsteps.second)
    {
        if(it.src1.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopkets_[it.src1.am()].insert(it.src1);
        if(it.src2.flags & DOUBLET_HRRTOPLEVEL)
            hrrtopkets_[it.src2.am()].insert(it.src2);
    }

    // we may need to add ( a s | as a top bra for AM quartets
    // that do not have HRR in the bra part
    // these might be ( X s |  or  | X s )
    if(hrrtopbras_.size() == 0)
    {   
        GaussianSet gs = AllGaussiansForAM(finalam[0]);
        for(const auto & it : gs)
            hrrtopbras_[it.am()].insert({DoubletType::BRA, it, {0,0,0}, DOUBLET_INITIAL | DOUBLET_HRRTOPLEVEL});
    }

    if(hrrtopkets_.size() == 0)
    {   
        GaussianSet gs = AllGaussiansForAM(finalam[2]);
        for(const auto & it : gs)
            hrrtopkets_[it.am()].insert({DoubletType::KET, it, {0,0,0}, DOUBLET_INITIAL | DOUBLET_HRRTOPLEVEL});
    }


    // create the top quartets
    for(const auto & it : hrrtopbras_)
    for(const auto & it2 : hrrtopkets_)
    {
        for(const auto & dit : it.second)
        for(const auto & dit2 : it2.second)
            hrrtopquartets_.insert({dit, dit2, 0, QUARTET_HRRTOPLEVEL});
    }
}



const DoubletSetMap & HRRWriter::TopBras(void) const
{
    return hrrtopbras_;
}



const DoubletSetMap & HRRWriter::TopKets(void) const
{
    return hrrtopkets_;
}



const QuartetSet & HRRWriter::TopQuartets(void) const
{
    return hrrtopquartets_;
}


void HRRWriter::WriteHRR(std::ostream & os, const WriterBase & base) const
{
    QAMList finalam = base.FinalAM();
    const int ncart_bra = NCART(finalam[0]) * NCART(finalam[1]);

    if(hrrsteps_.first.size() > 0)
    {
        os << "    //////////////////////////////////////////////\n";
        os << "    // Contracted integrals: Horizontal recurrance\n";
        os << "    // Bra part\n";
        os << "    // Steps: " << hrrsteps_.first.size() << "\n";
        os << "    //////////////////////////////////////////////\n";
        os << "\n";
        os << "    #pragma simd\n";
        os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << "    {\n";

        for(const auto & it : hrrtopkets_)
        {
            os << "        // form " << base.ArrVarName({finalam[0], finalam[1], it.first, 0}) << "\n";
            os << "        for(iket = 0; iket < " << it.second.size() << "; ++iket)\n";
            os << "        {\n";
            for(const auto & hit : hrrsteps_.first)
            {
                os << std::string(12, ' ') << "// " << hit << "\n";
                os << HRRBraStepString(hit, it.first, base) << "\n\n";
            }
            os << "        }\n";
            os << "\n";
        }

        os << "\n";
        os << "    }\n";
        os << "\n";
        os << "\n";
    }

    if(hrrsteps_.second.size() > 0)
    {
        os << "    //////////////////////////////////////////////\n";
        os << "    // Contracted integrals: Horizontal recurrance\n";
        os << "    // Ket part\n";
        os << "    // Steps: " << hrrsteps_.second.size() << "\n";
        os << "    //////////////////////////////////////////////\n";
        os << "\n";

        DAMList braam{finalam[0], finalam[1]};

        os << "    #pragma simd\n";
        os << "    for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << "    {\n";
        os << "        for(ibra = 0; ibra < " << ncart_bra << "; ++ibra)\n"; 
        os << "        {\n"; 

        for(const auto & hit : hrrsteps_.second)
        {
            os << std::string(12, ' ') << "// " << hit << "\n";
            os << HRRKetStepString(hit, braam, base) << "\n\n";
        }

        os << "        }\n"; 
        os << "    }\n";
    }

    os << "\n";
    os << "\n";
}



bool HRRWriter::HasBraHRR(void) const
{
    return (hrrsteps_.first.size() > 0);
}



bool HRRWriter::HasKetHRR(void) const
{
    return (hrrsteps_.second.size() > 0);
}



bool HRRWriter::HasHRR(void) const
{
    return ( HasBraHRR() || HasKetHRR() );
}



std::string HRRWriter::HRRBraStepArrVar(const Doublet & d, int ketam, bool istarget, const WriterBase & base) const
{
    if(base.IsContArray({d.left.am(), d.right.am(), ketam, 0}))
    {
        int ncart_ket = NCART(ketam); // all am is on the left part of the ket

        std::stringstream ss;
        ss << base.ArrVarName({d.left.am(), d.right.am(), ketam, 0})
           << "[" "abcd * " << d.ncart() * ncart_ket << " + " << d.idx() << " * " << ncart_ket
           << " + iket]"; 

        return ss.str();
    }
    else
    {
        std::stringstream ss;
        if(istarget)
            ss << "const double ";
        ss << "Q_" << base.amchar(d.left.ijk[0])  << "_" << base.amchar(d.left.ijk[1])  << "_" << base.amchar(d.left.ijk[2]) << "_"
                   << base.amchar(d.right.ijk[0]) << "_" << base.amchar(d.right.ijk[1]) << "_" << base.amchar(d.right.ijk[2]) << "_"
                   << base.amchar(ketam);
        return ss.str();
    }
}



std::string HRRWriter::HRRKetStepArrVar(const Doublet & d, const DAMList & braam, bool istarget, const WriterBase & base) const
{
    if(base.IsContArray({braam[0], braam[1], d.left.am(), d.right.am()}))
    {
        const int ncart_bra = NCART(braam[0]) * NCART(braam[1]);

        std::stringstream ss;
        ss << base.ArrVarName({braam[0], braam[1], d.left.am(), d.right.am()})
           << "[abcd * " << ncart_bra * d.ncart() << " + ibra * " << d.ncart() << " + " << d.idx() << "]"; 
        return ss.str();
    }
    else
    {
        std::stringstream ss;
        if(istarget)
            ss << "const double ";
        ss << "Q_" << base.amchar(d.left.ijk[0])  << "_" << base.amchar(d.left.ijk[1])  << "_" << base.amchar(d.left.ijk[2]) << "_"
                   << base.amchar(d.right.ijk[0]) << "_" << base.amchar(d.right.ijk[1]) << "_" << base.amchar(d.right.ijk[2]) << "_"
                   << base.amchar(braam[0]) << "_" << base.amchar(braam[1]);
        return ss.str();
    }
}



std::string HRRWriter::HRRBraStepString(const HRRDoubletStep & hrr, int ketam, const WriterBase & base) const
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "AB_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRBraStepArrVar(hrr.target, ketam, true, base);

    ss << " = ";
    ss << HRRBraStepArrVar(hrr.src1, ketam, false, base);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRBraStepArrVar(hrr.src2, ketam, false, base);
    ss << " );";

    return ss.str();
}



std::string HRRWriter::HRRKetStepString(const HRRDoubletStep & hrr, const DAMList & braam, const WriterBase & base) const
{
    // determine P,Q, etc, for AB_x, AB_y, AB_z
    const char * xyztype = "CD_";

    std::stringstream ss;
    ss << std::string(12, ' ');

    ss << HRRKetStepArrVar(hrr.target, braam, true, base);

    ss << " = ";
    ss << HRRKetStepArrVar(hrr.src1, braam, false, base);
    ss << " + ( " << xyztype << hrr.xyz << "[abcd] * ";
    ss << HRRKetStepArrVar(hrr.src2, braam, false, base);
    ss << " );";

    return ss.str();
}
