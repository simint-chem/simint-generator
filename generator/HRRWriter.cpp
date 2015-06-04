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


    // pointers to results for the HRR Bra step
    // I need the final bra, plus any where
    // the AM on the right is zero (created from VRR)
    brahrr_ptrs_.insert({finalam[0], finalam[1]});

    for(const auto & it : hrrsteps_.first)
    {
        if(it.src1.right.am() == 0)
            brahrr_ptrs_.insert({it.src1.left.am(), 0});
        if(it.src2.right.am() == 0)
            brahrr_ptrs_.insert({it.src2.left.am(), 0});
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



void HRRWriter::WriteHRRBraSteps_(std::ostream & os, const WriterBase & base, const std::string & ncart_ket) const
{
    os << "\n";
    os << "        for(iket = 0; iket < " << ncart_ket << "; ++iket)\n";
    os << "        {\n";

    for(const auto & hit : hrrsteps_.first)
    {
        os << std::string(12, ' ') << "// " << hit << "\n";
    
        // determine P,Q, etc, for AB_x, AB_y, AB_z
        const char * xyztype = "AB_";

        std::stringstream ss;
        os << std::string(12, ' ');

        os << HRRBraStepArrVar_(hit.target, ncart_ket, true, base);

        os << " = ";
        os << HRRBraStepArrVar_(hit.src1, ncart_ket, false, base);
        os << " + ( " << xyztype << hit.xyz << "[abcd] * ";
        os << HRRBraStepArrVar_(hit.src2, ncart_ket, false, base);
        os << " );";
        os << "\n\n";
    }

    os << "        }\n";
    os << "\n";
}



//void HRRWriter::WriteHRRKetSteps_(std::ostream & os, const WriterBase & base, int ketam) const
//{
//std::string HRRWriter::HRRKetStepString_(const HRRDoubletStep & hrr, const DAMList & braam, const WriterBase & base) const
//}



void HRRWriter::WriteHRRInline_(std::ostream & os, const WriterBase & base) const
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

        // Declare pointers for all contracted targets/srcs
        for(auto & it : brahrr_ptrs_)
        {
            if(it[0] == finalam[0] && it[1] == finalam[1])
                os << "        double * restrict";
            else
                os << "        double const * restrict";
            os << " BRA_" << amchar[it[0]] << "_" << amchar[it[1]] << ";\n";
        }

        os << "\n";

        for(const auto & it : hrrtopkets_)
        {
            // it.first is the AM for the ket part
            os << "        // form " << base.ArrVarName({base.FinalAM()[0], base.FinalAM()[1], it.first, 0}) << "\n";

            // set up the pointers needed
            for(const auto & itb : brahrr_ptrs_)
            {
                    os << "        BRA_" << amchar[itb[0]] << "_" << amchar[itb[1]]
                    << " = " << base.ArrVarName({itb[0], itb[1], it.first, 0}) 
                    << " + ( abcd * " << NCART(itb[0]) * NCART(itb[1]) * NCART(it.first)
                    << " );\n";
            }

            // ncart_ket in string form
            std::stringstream ss;
            ss << NCART(it.first);
    
            // actually write out the steps now
            WriteHRRBraSteps_(os, base, ss.str());
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

            // determine P,Q, etc, for AB_x, AB_y, AB_z
            const char * xyztype = "CD_";

            std::stringstream ss;
            os << std::string(12, ' ');
    
            os << HRRKetStepArrVar_(hit.target, braam, true, base);

            os << " = ";
            os << HRRKetStepArrVar_(hit.src1, braam, false, base);
            os << " + ( " << xyztype << hit.xyz << "[abcd] * ";
            os << HRRKetStepArrVar_(hit.src2, braam, false, base);
            os << " );";
            os << "\n\n";
        }

        os << "        }\n"; 
        os << "    }\n";
    }

    os << "\n";
    os << "\n";
}


void HRRWriter::WriteHRRExternal_(std::ostream & os, const WriterBase & base) const
{
    // todo
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



std::string HRRWriter::HRRBraStepArrVar_(const Doublet & d, const std::string & ncart_ket, bool istarget, const WriterBase & base) const
{
    std::stringstream ss;

    QAMList fam = base.FinalAM();
    if( (d.left.am() == fam[0] && d.right.am() == fam[1]) || d.right.am() == 0)
    {
        ss << "BRA_" << amchar[d.left.am()] << "_" << amchar[d.right.am()]
           << "[" << d.idx() << " * " << ncart_ket << " + iket]";
    }
    else
    {
        if(istarget)
            ss << "const double ";
        ss << "B_" << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2] << "_"
                   << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2];
    }

    return ss.str();
}



std::string HRRWriter::HRRKetStepArrVar_(const Doublet & d, const DAMList & braam, bool istarget, const WriterBase & base) const
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
        ss << "Q_" << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2] << "_"
                   << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2];
        return ss.str();
    }
}


void HRRWriter::WriteHRR(std::ostream & os, const WriterBase & base) const
{
    if(base.GetOption(OPTION_INLINEHRR) > 0)
        WriteHRRInline_(os, base);
    else
        WriteHRRExternal_(os, base);
}
        

void HRRWriter::WriteHRRFile(std::ostream & os, const WriterBase & base) const
{
}



void HRRWriter::WriteHRRHeaderFile(std::ostream & os, const WriterBase & base) const
{
}

