#include "generator/HRR_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/Helpers.hpp"
#include "generator/HRR_Algorithm_Base.hpp"
#include "generator/Ncart.hpp"


HRR_Writer::HRR_Writer(const HRR_Algorithm_Base & hrr_algo) 
    : hrr_algo_(hrr_algo)
{
}




void HRR_Writer::AddConstants(void) const
{
}



void HRR_Writer::WriteBraSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                                const std::string & ncart_ket, const std::string & ketstr) const
{
    os << indent4 << "for(iket = 0; iket < " << ncart_ket << "; ++iket)\n";
    os << indent4 << "{\n";

    for(const auto & it : steps)
    {
        os << std::string(20, ' ') << "// " << it << "\n";
    
        const char * xyztype = "hAB_";

        std::stringstream ss;
        os << std::string(20, ' ');

        os << HRRBraStepVar_(it.target, ncart_ket, ketstr);

        os << " = ";
        os << HRRBraStepVar_(it.src[0], ncart_ket, ketstr);
        os << " + ( " << xyztype << it.xyz << " * ";
        os << HRRBraStepVar_(it.src[1], ncart_ket, ketstr);
        os << " );";
        os << "\n\n";
    }

    os << indent4 << "}\n";
    os << "\n";
}



void HRR_Writer::WriteKetSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                                const std::string & ncart_bra, const std::string & brastr) const
{
        os << indent4 << "for(ibra = 0; ibra < " << ncart_bra << "; ++ibra)\n"; 
        os << indent4 << "{\n"; 
        for(const auto & it : steps)
        {
            os << std::string(20, ' ') << "// " << it << "\n";

            const char * xyztype = "hCD_";

            os << std::string(20, ' ');
    
            os << HRRKetStepVar_(it.target, brastr);

            os << " = ";
            os << HRRKetStepVar_(it.src[0], brastr);
            os << " + ( " << xyztype << it.xyz << " * ";
            os << HRRKetStepVar_(it.src[1], brastr);
            os << " );";
            os << "\n\n";
        }

        os << indent4 << "}\n"; 
}



std::string HRR_Writer::HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, 
                                       const std::string & ketstr) const
{
    std::stringstream ss;

    std::string arrname = WriterInfo::ArrVarName(d.left.am(), d.right.am(), ketstr);
    ss << "HRR_" << arrname << "[" << d.idx() << " * " << ncart_ket << " + iket]";

    return ss.str();
}



std::string HRR_Writer::HRRKetStepVar_(const Doublet & d, const std::string & brastr) const
{
    std::stringstream ss;

    std::string arrname = WriterInfo::ArrVarName(brastr, d.left.am(), d.right.am());
    ss << "HRR_" << arrname << "[ibra * " << d.ncart() << " + " << d.idx() << "]"; 

    return ss.str();
}




void HRR_Writer::WriteHRRInline_(std::ostream & os) const
{
    QAM finalam = WriterInfo::FinalAM();
    DAM finalbra{finalam[0], finalam[1]};

    if(WriterInfo::HasHRR())
    {
        os << indent3 << "//////////////////////////////////////////////\n";
        os << indent3 << "// Contracted integrals: Horizontal recurrance\n";
        os << indent3 << "//////////////////////////////////////////////\n";
        os << "\n";

        os << "\n";


        if(!WriterInfo::Scalar())
            os << indent3 << "#pragma omp simd linear(real_abcd)\n";
        os << indent3 << "for(abcd = 0; abcd < nshellbatch; ++abcd, ++real_abcd)\n";
        os << indent3 << "{\n";

        os << "\n";
        os << indent4 << "// set up HRR pointers\n";
        for(const auto & it : hrr_algo_.TopAM())
            os << indent4 << "double const * restrict HRR_" << WriterInfo::ArrVarName(it) << " = " << WriterInfo::ArrVarName(it) << " + abcd * " << NCART(it) << ";\n";

        // and also for the final integral
        os << indent4 << "double * restrict HRR_" << WriterInfo::ArrVarName(finalam) << " = " << WriterInfo::ArrVarName(finalam) << " + real_abcd * " << NCART(finalam) << ";\n";
        os << "\n";



        if(WriterInfo::HasBraHRR())
        {
            os << "\n";
            os << indent4 << "const double hAB_x = AB_x[abcd];\n";
            os << indent4 << "const double hAB_y = AB_y[abcd];\n";
            os << indent4 << "const double hAB_z = AB_z[abcd];\n";
            os << "\n";

            os << "\n";

            for(const auto & itk : hrr_algo_.TopKetAM()) // for all needed ket am
            for(const auto & itb : hrr_algo_.GetBraAMOrder()) // form these
            {
                os << indent4 << "// form " << WriterInfo::ArrVarName({itb[0], itb[1], itk[0], itk[1]}) << "\n";

                // allocate temporary if needed
                QAM am = {itb[0], itb[1], itk[0], itk[1]};
                if(!WriterInfo::IsContArray(am) && !WriterInfo::IsFinalAM(am))
                    os << indent4 << "double HRR_" << WriterInfo::ArrVarName(am) << "[" << NCART(am) << "];\n";

                // ncart_ket in string form
                std::stringstream ss;
                ss << NCART(itk[0], itk[1]);
    
                // the ket part in string form
                std::stringstream ssket;
                ssket << amchar[itk[0]] << "_" << amchar[itk[1]];
        
                // actually write out the steps now
                WriteBraSteps_(os, hrr_algo_.GetBraSteps(itb), ss.str(), ssket.str());
            }
        }

        os << "\n";
        os << "\n";

        if(WriterInfo::HasKetHRR())
        {
            os << "\n";
            os << indent4 << "const double hCD_x = CD_x[abcd];\n";
            os << indent4 << "const double hCD_y = CD_y[abcd];\n";
            os << indent4 << "const double hCD_z = CD_z[abcd];\n";
            os << "\n";

            for(const auto & itk : hrr_algo_.GetKetAMOrder())
            {
                os << indent4 << "// form " << WriterInfo::ArrVarName({finalbra[0], finalbra[1], itk[0], itk[1]}) << "\n";

                // allocate temporary if needed
                QAM am = {finalbra[0], finalbra[1], itk[0], itk[1]};
                if(!WriterInfo::IsContArray(am) && !WriterInfo::IsFinalAM(am))
                    os << indent4 << "double HRR_" << WriterInfo::ArrVarName(am) << "[" << NCART(am) << "];\n";

                // ncart_bra in string form
                std::stringstream ss;
                ss << NCART(finalbra); 

                // the bra part in string form
                std::stringstream ssbra;
                ssbra << amchar[finalbra[0]] << "_" << amchar[finalbra[1]];

                WriteKetSteps_(os, hrr_algo_.GetKetSteps(itk), ss.str(), ssbra.str());
            }
        }

        os << "\n";
        os << indent3 << "}  // close HRR loop\n";

    }

    os << "\n";
    os << "\n";
}


void HRR_Writer::WriteHRRExternal_(std::ostream & os) const
{
    QAM finalam = WriterInfo::FinalAM();
    DAM finalbra{finalam[0], finalam[1]};

    if(WriterInfo::HasHRR())
    {
        os << indent3 << "//////////////////////////////////////////////\n";
        os << indent3 << "// Contracted integrals: Horizontal recurrance\n";
        os << indent3 << "//////////////////////////////////////////////\n";
        os << "\n";

        os << "\n";


        if(!WriterInfo::Scalar())
            os << indent3 << "#pragma omp simd linear(real_abcd)\n";
        os << indent3 << "for(abcd = 0; abcd < nshellbatch; ++abcd, ++real_abcd)\n";
        os << indent3 << "{\n";

        os << "\n";
        os << indent4 << "// set up HRR pointers\n";
        for(const auto & it : hrr_algo_.TopAM())
            os << indent4 << "double const * restrict HRR_" << WriterInfo::ArrVarName(it) << " = " << WriterInfo::ArrVarName(it) << " + abcd * " << NCART(it) << ";\n";

        // and also for the final integral
        os << indent4 << "double * restrict HRR_" << WriterInfo::ArrVarName(finalam) << " = " << WriterInfo::ArrVarName(finalam) << " + real_abcd * " << NCART(finalam) << ";\n";
        os << "\n";



        if(WriterInfo::HasBraHRR())
        {
            os << "\n";
            os << indent4 << "const double hAB_x = AB_x[abcd];\n";
            os << indent4 << "const double hAB_y = AB_y[abcd];\n";
            os << indent4 << "const double hAB_z = AB_z[abcd];\n";
            os << "\n";

            os << "\n";

            for(const auto & itk : hrr_algo_.TopKetAM()) // for all needed ket am
            for(const auto & itb : hrr_algo_.GetBraAMOrder()) // form these
            {
                os << indent4 << "// form " << WriterInfo::ArrVarName({itb[0], itb[1], itk[0], itk[1]}) << "\n";

                // allocate temporary if needed
                QAM am = {itb[0], itb[1], itk[0], itk[1]};
                if(!WriterInfo::IsContArray(am) && !WriterInfo::IsFinalAM(am))
                    os << indent4 << "double HRR_" << WriterInfo::ArrVarName(am) << "[" << NCART(am) << "];\n";

                // call function
                os << indent4 << "HRR_BRA_" << amchar[itb[0]] << "_" << amchar[itb[1]] << "(\n";

                // pointer to result buffer
                os << indent5 << "HRR_" << WriterInfo::ArrVarName({itb[0], itb[1], itk[0], itk[1]}) << ",\n";

                // pointer to requirements
                for(const auto & it : hrr_algo_.GetBraAMReq(itb))
                    os << indent5 << "HRR_" << WriterInfo::ArrVarName({it[0], it[1], itk[0], itk[1]}) << ",\n";

                os << indent5 << "hAB_x, hAB_y, hAB_z, " << NCART(itk[0], itk[1]) << ");\n"; 
                os << "\n\n";
            }
        }

        os << "\n";
        os << "\n";

        if(WriterInfo::HasKetHRR())
        {
            os << "\n";
            os << indent4 << "const double hCD_x = CD_x[abcd];\n";
            os << indent4 << "const double hCD_y = CD_y[abcd];\n";
            os << indent4 << "const double hCD_z = CD_z[abcd];\n";
            os << "\n";

            for(const auto & itk : hrr_algo_.GetKetAMOrder())
            {
                os << indent4 << "// form " << WriterInfo::ArrVarName({finalbra[0], finalbra[1], itk[0], itk[1]}) << "\n";

                // allocate temporary if needed
                QAM am = {finalbra[0], finalbra[1], itk[0], itk[1]};
                if(!WriterInfo::IsContArray(am) && !WriterInfo::IsFinalAM(am))
                    os << indent4 << "double HRR_" << WriterInfo::ArrVarName(am) << "[" << NCART(am) << "];\n";

                os << indent4 << "HRR_KET_";
                os << amchar[itk[0]] << "_" << amchar[itk[1]] << "(\n";

                // pointer to result buffer
                os << indent5 << "HRR_" << WriterInfo::ArrVarName({finalbra[0], finalbra[1], itk[0], itk[1]}) << ",\n";

                // pointer to requirements
                for(const auto & it : hrr_algo_.GetKetAMReq(itk))
                    os << indent5 << "HRR_" << WriterInfo::ArrVarName({finalbra[0], finalbra[1], it[0], it[1]}) << ",\n";

                os << indent5 << "hCD_x, hCD_y, hCD_z, " << NCART(finalbra) << ");\n";
                os << "\n\n";

            }
        }

        os << "\n";
        os << indent3 << "}  // close HRR loop\n";

    }

    os << "\n";
    os << "\n";
}



void HRR_Writer::WriteHRR(std::ostream & os) const
{
    if(WriterInfo::GetOption(OPTION_INLINEHRR) > 0)
        WriteHRRInline_(os);
    else
        WriteHRRExternal_(os);
}



void HRR_Writer::WriteHRRFile(std::ostream & ofb, std::ostream & ofk, std::ostream & ofh) const
{
    QAM finalam = WriterInfo::FinalAM();
    DAM braam = {finalam[0], finalam[1]};
    DAM ketam = {finalam[2], finalam[3]};

    // only do final
    HRRDoubletStepList brasteps = hrr_algo_.GetBraSteps(braam);
    HRRDoubletStepList ketsteps = hrr_algo_.GetKetSteps(ketam);

    if(brasteps.size() > 0)
    {
        ofb << "//////////////////////////////////////////////\n";
        ofb << "// BRA: ( " << amchar[braam[0]] << " " << amchar[braam[1]] << " |\n";
        ofb << "// Steps: " << brasteps.size() << "\n";
        ofb << "//////////////////////////////////////////////\n";
        ofb << "\n";

        std::stringstream prototype;

        if(!WriterInfo::Scalar())
            prototype << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(ncart_ket)\n";
        prototype << "void HRR_BRA_";
        prototype << amchar[braam[0]] << "_" << amchar[braam[1]] << "(\n";

        // pointer to result buffer
        prototype << indent5 << "double * const restrict HRR_" << WriterInfo::ArrVarName(finalam[0], finalam[1], "X_X") << ",\n";

        // pointer to requirements
        for(const auto & it : hrr_algo_.GetBraAMReq(braam))
            prototype << indent5 << "double const * const restrict HRR_" << WriterInfo::ArrVarName(it[0], it[1], "X_X") << ",\n";

        prototype << indent5 << "const double hAB_x, const double hAB_y, const double hAB_z, const int ncart_ket)";



        ofb << prototype.str() << "\n";
        ofb << "{\n";
        ofb << indent1 << "int iket;\n";
        ofb << "\n";

        WriteBraSteps_(ofb, brasteps, "ncart_ket", "X_X"); 

        ofb << "\n";
        ofb << "}\n";
        ofb << "\n";
        ofb << "\n";


        // header
        ofh << prototype.str() << ";\n";

    }

    if(ketsteps.size() > 0)
    {
        ofk << "//////////////////////////////////////////////\n";
        ofk << "// KET: ( " << amchar[ketam[0]] << " " << amchar[ketam[1]] << " |\n";
        ofk << "// Steps: " << ketsteps.size() << "\n";
        ofk << "//////////////////////////////////////////////\n";
        ofk << "\n";

        std::stringstream prototype;

        if(!WriterInfo::Scalar())
            prototype << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(ncart_bra)\n";
        prototype << "void HRR_KET_";
        prototype << amchar[ketam[0]] << "_" << amchar[ketam[1]] << "(\n";

        // pointer to result buffer
        prototype << indent5 << "double * const restrict HRR_" << WriterInfo::ArrVarName("X_X", finalam[2], finalam[3]) << ",\n";

        // pointer to requirements
        for(const auto & it : hrr_algo_.GetKetAMReq(ketam))
            prototype << indent5 << "double const * const restrict HRR_" << WriterInfo::ArrVarName("X_X", it[0], it[1]) << ",\n";

        prototype << indent5 << "const double hCD_x, const double hCD_y, const double hCD_z, const int ncart_bra)\n";



        ofk << prototype.str() << "\n";
        ofk << "{\n";
        ofk << indent1 << "int ibra;\n";
        ofk << "\n";

        WriteKetSteps_(ofk, ketsteps, "ncart_bra", "X_X"); 

        ofk << "\n";
        ofk << "}\n";
        ofk << "\n";
        ofk << "\n";


        // header
        ofh << prototype.str() << ";\n";
    }
}


