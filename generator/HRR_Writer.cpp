#include "generator/HRR_Writer.hpp"
#include "generator/WriterInfo.hpp"
#include "generator/Helpers.hpp"
#include "generator/HRR_Algorithm_Base.hpp"
#include "generator/Ncart.hpp"


HRR_Writer::HRR_Writer(const HRR_Algorithm_Base & hrr_algo) 
{
    brakettop_ = hrr_algo.TopBraKet();
    brakettopam_ = hrr_algo.TopBraKetAM();
    hrrsteps_ = hrr_algo.DoubletStepLists();
    topquartetam_ = hrr_algo.TopQAM();
}



void HRR_Writer::WriteIncludes(std::ostream & os) const
{
}



void HRR_Writer::AddConstants(void) const
{
}



void HRR_Writer::WriteBraSteps_(std::ostream & os, const std::string & ncart_ket, const std::string & ketstr) const
{
    os << indent4 << "for(iket = 0; iket < " << ncart_ket << "; ++iket)\n";
    os << indent4 << "{\n";

    for(const auto & hit : hrrsteps_.first)
    {
        os << std::string(20, ' ') << "// " << hit << "\n";
    
        const char * xyztype = "AB_";

        std::stringstream ss;
        os << std::string(20, ' ');

        os << HRRBraStepVar_(hit.target, ncart_ket, ketstr, true);

        os << " = ";
        os << HRRBraStepVar_(hit.src[0], ncart_ket, ketstr, false);
        os << " + ( " << xyztype << hit.xyz << " * ";
        os << HRRBraStepVar_(hit.src[1], ncart_ket, ketstr, false);
        os << " );";
        os << "\n\n";
    }

    os << indent4 << "}\n";
    os << "\n";
}



void HRR_Writer::WriteKetSteps_(std::ostream & os, const std::string & ncart_bra, const std::string & brastr) const
{
        os << indent4 << "for(ibra = 0; ibra < " << ncart_bra << "; ++ibra)\n"; 
        os << indent4 << "{\n"; 
        for(const auto & hit : hrrsteps_.second)
        {
            os << std::string(20, ' ') << "// " << hit << "\n";

            const char * xyztype = "CD_";

            os << std::string(20, ' ');
    
            os << HRRKetStepVar_(hit.target, ncart_bra, brastr, true);

            os << " = ";
            os << HRRKetStepVar_(hit.src[0], ncart_bra, brastr, false);
            os << " + ( " << xyztype << hit.xyz << " * ";
            os << HRRKetStepVar_(hit.src[1], ncart_bra, brastr, false);
            os << " );";
            os << "\n\n";
        }

        os << indent4 << "}\n"; 
}



std::string HRR_Writer::HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, 
                                       const std::string & ketstr, bool istarget) const
{
    std::stringstream ss;

    QAM finalam = WriterInfo::FinalAM();

    // is final doublet?
    bool isfinald = (d.left.am() == finalam[0] && d.right.am() == finalam[1]);

    // is top level bra? (ie result from VRR/ET)
    bool istop = ( d.right.am() == 0);

    // actually the final quartet?
    bool isfinalq = isfinald && !WriterInfo::HasKetHRR();

    std::string arrname = WriterInfo::ArrVarName(d.left.am(), d.right.am(), ketstr);

    if(istop || isfinalq || isfinald)
        ss << "HRR_" << arrname << "[" << d.idx() << " * " << ncart_ket << " + iket]";
    else
    {
        if(istarget)
            ss << "const double ";
        ss << "B_" << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2] << "_"
                   << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2] << "_" << ketstr;
    }

    return ss.str();
}



std::string HRR_Writer::HRRKetStepVar_(const Doublet & d, const std::string & ncart_bra,
                                       const std::string & brastr, bool istarget) const
{
    std::stringstream ss;

    QAM finalam = WriterInfo::FinalAM();

    // is final doublet? This would be the final quartet also
    bool isfinal = (d.left.am() == finalam[2] && d.right.am() == finalam[3]);

    // is from VRR or ET
    // Since this is run last, the bra part is ( finalam[0] finalam[1] |
    // so if that is ( X s | and this is | X s ) it is from VRR/ET
    bool istop = ( d.right.am() == 0 && finalam[1] == 0 );

    // is top level ket? (ie result from bra, but not VRR, etc)
    bool istopbra = ( !istop && d.right.am() == 0 );

    std::string arrname = WriterInfo::ArrVarName(brastr, d.left.am(), d.right.am());

    if(istopbra || istop || isfinal)
        ss << "HRR_" << arrname << "[ibra * " << d.ncart() << " + " << d.idx() << "]"; 
    else
    {
        if(istarget)
            ss << "const double ";
        ss << "K_" << brastr << "_" << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2] << "_"
                   << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2];
    }
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
        for(const auto & it : topquartetam_)
            os << indent4 << "double const * restrict HRR_" << WriterInfo::ArrVarName(it) << " = " << WriterInfo::ArrVarName(it) << " + abcd * " << NCART(it) << ";\n";

        // and also for the final integral
        os << indent4 << "double * restrict HRR_" << WriterInfo::ArrVarName(finalam) << " = " << WriterInfo::ArrVarName(finalam) << " + real_abcd * " << NCART(finalam) << ";\n";
        os << "\n";



        if(WriterInfo::HasBraHRR())
        {
            os << "\n";
            os << indent4 << "const double AB_x = AB_x[abcd];\n";
            os << indent4 << "const double AB_y = AB_y[abcd];\n";
            os << indent4 << "const double AB_z = AB_z[abcd];\n";
            os << "\n";

            // allocate HRR Temporaries. These are the top ket (left) AM, with the final bra
            // but skip the final AM. we store directly there
            for(const auto & it : brakettopam_.second)
            {
                QAM qam{finalbra[0], finalbra[1], it[0], 0};
                if(!WriterInfo::IsFinalAM(qam))
                {
                    os << indent4 << "double HRR_" << WriterInfo::ArrVarName(qam)
                       << "[" << NCART(finalbra[0], finalbra[1], it[0]) << "];\n";
                }
            }

            os << "\n";

            for(const auto & it : brakettopam_.second)
            {
                // it.first is the AM for the ket part
                os << indent4 << "// form " << WriterInfo::ArrVarName({finalbra[0], finalbra[1], it[0], 0}) << "\n";

                // ncart_ket in string form
                std::stringstream ss;
                ss << NCART(it[0]);
    
                // the ket part in string form
                std::stringstream ssket;
                ssket << amchar[it[0]] << "_s";
        
                // actually write out the steps now
                WriteBraSteps_(os, ss.str(), ssket.str());
            }
        }

        os << "\n";
        os << "\n";

        if(WriterInfo::HasKetHRR())
        {
            os << "\n";
            os << indent4 << "const double CD_x = CD_x[abcd];\n";
            os << indent4 << "const double CD_y = CD_y[abcd];\n";
            os << indent4 << "const double CD_z = CD_z[abcd];\n";
            os << "\n";

            // ncart_bra in string form
            std::stringstream ss;
            ss << NCART(finalbra); 

            // the bra part in string form
            std::stringstream ssbra;
            ssbra << amchar[finalbra[0]] << "_" << amchar[finalbra[1]];

            WriteKetSteps_(os, ss.str(), ssbra.str());

    
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
        for(const auto & it : topquartetam_)
            os << indent4 << "double const * restrict HRR_" << WriterInfo::ArrVarName(it) << " = " << WriterInfo::ArrVarName(it) << " + abcd * " << NCART(it) << ";\n";

        // and also for the final integral
        os << indent4 << "double * restrict HRR_" << WriterInfo::ArrVarName(finalam) << " = " << WriterInfo::ArrVarName(finalam) << " + real_abcd * " << NCART(finalam) << ";\n";
        os << "\n";

        if(WriterInfo::HasBraHRR())
        {
            // allocate HRR Temporaries. These are the top ket (left) AM, with the final bra
            // but skip the final AM. we store directly there
            for(const auto & it : brakettopam_.second)
            {
                QAM qam{finalbra[0], finalbra[1], it[0], 0};
                if(!WriterInfo::IsFinalAM(qam))
                {
                    os << indent4 << "double HRR_" << WriterInfo::ArrVarName(qam)
                       << "[" << NCART(finalbra[0], finalbra[1], it[0]) << "];\n";
                }
            }
            os << "\n";

            for(const auto & it : brakettopam_.second)
            {
                // it.first is the AM for the ket part

                // ncart_ket in string form
                // ( should be |X s) )
                std::stringstream ss;
                ss << NCART(it[0]);
    

                os << indent4 << "// form " << WriterInfo::ArrVarName({finalbra[0], finalbra[1], it[0], 0}) << "\n";
                os << indent4 << "HRR_BRA_" << amchar[finalbra[0]] << "_" << amchar[finalbra[1]] << "(\n";
                for(const auto & it2 : brakettopam_.first)
                   os << indent5 << "HRR_" << WriterInfo::ArrVarName({it2[0], it2[1], it[0], it[1]}) << ",\n";

                os << indent5 << "HRR_" << WriterInfo::ArrVarName({finalbra[0], finalbra[1], it[0], it[1]}) << ",\n";
                os << indent5 << "AB_x[abcd], AB_y[abcd], AB_z[abcd], " << ss.str() << ");\n";
                os << "\n";
            }
        }

        os << "\n";
        os << "\n";

        if(WriterInfo::HasKetHRR())
        {
            // ncart_bra in string form
            std::stringstream ss;
            ss << NCART(finalbra);

            os << indent4 << "// form " << WriterInfo::ArrVarName(finalam) << "\n";
            os << indent4 << "HRR_KET_" << amchar[finalam[2]] << "_" << amchar[finalam[3]] << "(\n";
            for(const auto & it2 : brakettopam_.second)
               os << indent5 << "HRR_" << WriterInfo::ArrVarName({finalbra[0], finalbra[1], it2[0], it2[1]}) << ",\n";

            os << indent5 << "HRR_" << WriterInfo::ArrVarName(finalam) << ",\n";
            os << indent5 << "CD_x[abcd], CD_y[abcd], CD_z[abcd], " << ss.str() << ");\n";
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



void HRR_Writer::WriteHRRFile(std::ostream & ofb, std::ostream & ofk) const
{
    QAM finalam = WriterInfo::FinalAM();


    if(hrrsteps_.first.size() > 0)
    {
        ofb << "//////////////////////////////////////////////\n";
        ofb << "// BRA: ( " << amchar[finalam[0]] << " " << amchar[finalam[1]] << " |\n";
        ofb << "// Steps: " << hrrsteps_.first.size() << "\n";
        ofb << "//////////////////////////////////////////////\n";
        ofb << "\n";

        // it.first is the AM for the ket part
        if(!WriterInfo::Scalar())
            ofb << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(ncart_ket)\n";
        ofb << "void HRR_BRA_";
        ofb << amchar[finalam[0]] << "_" << amchar[finalam[1]] << "(\n";

        // pointers to buffers (top bras)
        for(const auto & it : brakettopam_.first)
            ofb << indent5 << "double const * const restrict HRR_" << WriterInfo::ArrVarName(it[0], it[1], "X_X") << ",\n";

        // pointer to result buffer
        ofb << indent5 << "double * const restrict HRR_" << WriterInfo::ArrVarName(finalam[0], finalam[1], "X_X") << ",\n";

        ofb << indent5 << "const double AB_x, const double AB_y, const double AB_z, const int ncart_ket)\n";

        ofb << "{\n";
        ofb << indent1 << "int iket;\n";
        ofb << "\n";

        WriteBraSteps_(ofb, "ncart_ket", "X_X"); 

        ofb << "\n";
        ofb << "}\n";
        ofb << "\n";
        ofb << "\n";
    }

    if(hrrsteps_.second.size() > 0)
    {
        ofk << "//////////////////////////////////////////////\n";
        ofk << "// KET: ( " << amchar[finalam[0]] << " " << amchar[finalam[1]] << " |\n";
        ofk << "// Steps: " << hrrsteps_.second.size() << "\n";
        ofk << "//////////////////////////////////////////////\n";
        ofk << "\n";

        if(!WriterInfo::Scalar())
            ofk << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(ncart_bra)\n";
        ofk << "void HRR_KET_";
        ofk << amchar[finalam[2]] << "_" << amchar[finalam[3]] << "(\n";

        // pointers to buffers (top bras)
        for(const auto & it : brakettopam_.second)
            ofk << indent5 << "double const * const restrict HRR_" << WriterInfo::ArrVarName("X_X", it[0], it[1])  << ",\n";

        // pointer to result buffer
        ofk << indent5 << "double * const restrict HRR_" << WriterInfo::ArrVarName("X_X", finalam[2], finalam[3]) << ",\n";

        ofk << indent5 << "const double CD_x, const double CD_y, const double CD_z, const int ncart_bra)\n";
        ofk << "{\n";
        ofk << indent1 << "int ibra;\n";
        ofk << "\n";

        WriteKetSteps_(ofk, "ncart_bra", "X_X"); 

        ofk << "\n";
        ofk << "}\n";
        ofk << "\n";
        ofk << "\n";
    }
}



void HRR_Writer::WriteHRRHeaderFile(std::ostream & os) const
{
    QAM finalam = WriterInfo::FinalAM();

    if(hrrsteps_.first.size() > 0)
    {
        os << indent1 << "//////////////////////////////////////////////\n";
        os << indent1 << "// BRA: ( " << amchar[finalam[0]] << " " << amchar[finalam[1]] << " |\n";
        os << indent1 << "// Steps: " << hrrsteps_.first.size() << "\n";
        os << indent1 << "//////////////////////////////////////////////\n";

        // it.first is the AM for the ket part
        if(!WriterInfo::Scalar())
            os << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(ncart_ket)\n";

        os << "void HRR_BRA_";
        os << amchar[finalam[0]] << "_" << amchar[finalam[1]] << "(\n";

        // pointers to buffers (top bras)
        for(const auto & it : brakettopam_.first)
            os << indent5 << "double const * const restrict HRR_" << WriterInfo::ArrVarName(it[0], it[1], "X_X") << ",\n";

        // pointer to result buffer
        os << indent5 << "double * const restrict HRR_" << WriterInfo::ArrVarName(finalam[0], finalam[1], "X_X") << ",\n";

        os << indent5 << "const double AB_x, const double AB_y, const double AB_z, const int ncart_ket);\n\n\n";

    }

    if(hrrsteps_.second.size() > 0)
    {
        os << indent1 << "//////////////////////////////////////////////\n";
        os << indent1 << "// KET: ( " << amchar[finalam[0]] << " " << amchar[finalam[1]] << " |\n";
        os << indent1 << "// Steps: " << hrrsteps_.second.size() << "\n";
        os << indent1 << "//////////////////////////////////////////////\n";

        if(!WriterInfo::Scalar())
            os << "#pragma omp declare simd simdlen(SIMINT_SIMD_LEN) uniform(ncart_bra)\n";

        os << "void HRR_KET_";
        os << amchar[finalam[2]] << "_" << amchar[finalam[3]] << "(\n";

        // pointers to buffers (top bras)
        for(const auto & it : brakettopam_.second)
            os << indent5 << "double const * const restrict HRR_" << WriterInfo::ArrVarName("X_X", it[0], it[1])  << ",\n";

        // pointer to result buffer
        os << indent5 << "double * const restrict HRR_" << WriterInfo::ArrVarName("X_X", finalam[2], finalam[3]) << ",\n";

        os << indent5 << "const double CD_x, const double CD_y, const double CD_z, const int ncart_bra);\n\n\n";
    }
}

