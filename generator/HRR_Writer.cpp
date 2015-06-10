#include "generator/HRR_Writer.hpp"
#include "generator/WriterBase.hpp"
#include "generator/Helpers.hpp"
#include "generator/HRR_Algorithm_Base.hpp"


HRR_Writer::HRR_Writer(const HRR_Algorithm_Base & hrr_algo) 
{
    brakettop_ = hrr_algo.TopBraKet();
    brakettopam_ = hrr_algo.TopBraKetAM();
    hrrsteps_ = hrr_algo.DoubletStepLists();
    topquartetam_ = hrr_algo.TopQAM();
}



void HRR_Writer::WriteIncludes(std::ostream & os, const WriterBase & base) const
{
    if(base.GetOption(OPTION_INLINEHRR) == 0)
        os << "#include \"eri/hrr.gen/hrr.h\"\n";
}



void HRR_Writer::WriteBraSteps_(std::ostream & os, const WriterBase & base, const std::string & ncart_ket, const std::string & ketstr) const
{
    os << "\n";
    os << indent4 << "for(iket = 0; iket < " << ncart_ket << "; ++iket)\n";
    os << indent4 << "{\n";

    for(const auto & hit : hrrsteps_.first)
    {
        os << std::string(20, ' ') << "// " << hit << "\n";
    
        const char * xyztype = "AB_";

        std::stringstream ss;
        os << std::string(20, ' ');

        os << HRRBraStepVar_(hit.target, ncart_ket, ketstr, true, base);

        os << " = ";
        os << HRRBraStepVar_(hit.src1, ncart_ket, ketstr, false, base);
        os << " + ( b" << xyztype << hit.xyz << " * ";
        os << HRRBraStepVar_(hit.src2, ncart_ket, ketstr, false, base);
        os << " );";
        os << "\n\n";
    }

    os << indent4 << "}\n";
    os << "\n";
}



void HRR_Writer::WriteKetSteps_(std::ostream & os, const WriterBase & base, const std::string & ncart_bra, const std::string & brastr) const
{
        os << indent4 << "for(ibra = 0; ibra < " << ncart_bra << "; ++ibra)\n"; 
        os << indent4 << "{\n"; 
        for(const auto & hit : hrrsteps_.second)
        {
            os << std::string(20, ' ') << "// " << hit << "\n";

            const char * xyztype = "CD_";

            os << std::string(20, ' ');
    
            os << HRRKetStepVar_(hit.target, ncart_bra, brastr, true, base);

            os << " = ";
            os << HRRKetStepVar_(hit.src1, ncart_bra, brastr, false, base);
            os << " + ( k" << xyztype << hit.xyz << " * ";
            os << HRRKetStepVar_(hit.src2, ncart_bra, brastr, false, base);
            os << " );";
            os << "\n\n";
        }

        os << indent4 << "}\n"; 
}



std::string HRR_Writer::HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, 
                                       const std::string & ketstr, bool istarget, const WriterBase & base) const
{
    std::stringstream ss;

    QAM finalam = base.FinalAM();

    // is final doublet?
    bool isfinald = (d.left.am() == finalam[0] && d.right.am() == finalam[1]);

    // is top level bra? (ie result from VRR/ET)
    bool istop = ( d.right.am() == 0);

    // actually the final quartet?
    bool isfinalq = isfinald && !base.HasKetHRR();

    std::string arrname = base.ArrVarName(d.left.am(), d.right.am(), ketstr);

    if(istop || isfinalq)
        ss << arrname << "[idx_" << arrname << " + " << d.idx() << " * " << ncart_ket << " + iket]";
    else if(isfinald)
        ss << arrname << "[" << d.idx() << " * " << ncart_ket << " + iket]";
    else
    {
        if(istarget)
            ss << "const double ";
        ss << "B_" << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2] << "_"
                   << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2] << ketstr;
    }

    return ss.str();
}



std::string HRR_Writer::HRRKetStepVar_(const Doublet & d, const std::string & ncart_bra,
                                       const std::string & brastr, bool istarget, const WriterBase & base) const
{
    std::stringstream ss;

    QAM finalam = base.FinalAM();

    // is final doublet? This would be the final quartet also
    bool isfinal = (d.left.am() == finalam[2] && d.right.am() == finalam[3]);

    // is from VRR or ET
    // Since this is run last, the bra part is ( finalam[0] finalam[1] |
    // so if that is ( X s | and this is | X s ) it is from VRR/ET
    bool istop = ( d.right.am() == 0 && finalam[1] == 0 );

    // is top level ket? (ie result from bra, but not VRR, etc)
    bool istopbra = ( !istop && d.right.am() == 0 );

    std::string arrname = base.ArrVarName(brastr, d.left.am(), d.right.am());

    if(istopbra)
        ss << arrname << "[ibra * " << d.ncart() << " + " << d.idx() << "]"; 
    else if(istop || isfinal)
        ss << arrname << "[idx_" << arrname << " + ibra * " << d.ncart() << " + " << d.idx() << "]"; 
    else
    {
        if(istarget)
            ss << "const double ";
        ss << "K_" << brastr << "_" << d.left.ijk[0]  << "_" << d.left.ijk[1]  << "_" << d.left.ijk[2] << "_"
                   << d.right.ijk[0] << "_" << d.right.ijk[1] << "_" << d.right.ijk[2];
    }
    return ss.str();
}




void HRR_Writer::WriteHRRInline_(std::ostream & os, const WriterBase & base) const
{
    QAM finalam = base.FinalAM();
    DAM finalbra{finalam[0], finalam[1]};

    if(base.HasHRR())
    {
        os << indent3 << "//////////////////////////////////////////////\n";
        os << indent3 << "// Contracted integrals: Horizontal recurrance\n";
        os << indent3 << "//////////////////////////////////////////////\n";
        os << "\n";
        os << "            #pragma omp simd linear(real_abcd)\n";
        os << indent3 << "for(abcd = 0; abcd < nshell1234; ++abcd, ++real_abcd)\n";
        os << indent3 << "{\n";

        // variables
        os << "\n";
        if(base.HasBraHRR())
        {
            os << indent4 << "const double bAB_x = AB_x[abcd];\n";
            os << indent4 << "const double bAB_y = AB_y[abcd];\n";
            os << indent4 << "const double bAB_z = AB_z[abcd];\n";
        }
        if(base.HasKetHRR())
        {
            os << indent4 << "const double kCD_x = CD_x[abcd];\n";
            os << indent4 << "const double kCD_y = CD_y[abcd];\n";
            os << indent4 << "const double kCD_z = CD_z[abcd];\n";
        }
        os << "\n";


        // we may be able to do some index math up front
        for(const auto & it : topquartetam_)
            os << indent4 << "const int idx_" << base.ArrVarName(it) << " = abcd * " << NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]) << ";\n";

        // and also for the final integral
        os << indent4 << "const int idx_" << base.ArrVarName(finalam) << " = real_abcd * " << NCART(finalam[0]) * NCART(finalam[1]) * NCART(finalam[2]) * NCART(finalam[3]) << ";\n";
        os << "\n";

        if(base.HasBraHRR())
        {
            // allocate HRR Temporaries. These are the top ket (left) AM, with the final bra
            // but skip the final AM if we aren't permuting. They are stored directly
            for(const auto & it : brakettopam_.second)
            {
                QAM qam{finalbra[0], finalbra[1], it[0], 0};
                if(!base.IsFinalAM(qam))
                {
                    os << indent4 << "double " << base.ArrVarName(qam)
                       << "[" << NCART(finalbra[0]) * NCART(finalbra[1]) * NCART(it[0]) << "];\n";
                }
            }

            os << "\n";

            for(const auto & it : brakettopam_.second)
            {
                // it.first is the AM for the ket part
                os << indent4 << "// form " << base.ArrVarName({finalbra[0], finalbra[1], it[0], 0}) << "\n";

                // ncart_ket in string form
                std::stringstream ss;
                ss << NCART(it[0]);
    
                // the ket part in string form
                std::stringstream ssket;
                ssket << amchar[it[0]] << "_s";
        
                // actually write out the steps now
                WriteBraSteps_(os, base, ss.str(), ssket.str());
            }
        }

        os << "\n";
        os << "\n";

        if(base.HasKetHRR())
        {
            // ncart_bra in string form
            DAM braam{finalam[0], finalam[1]};
            std::stringstream ss;
            ss << (NCART(braam[0]) * NCART(braam[1])); 

            // the bra part in string form
            std::stringstream ssbra;
            ssbra << amchar[finalbra[0]] << "_" << amchar[finalbra[1]];

            WriteKetSteps_(os, base, ss.str(), ssbra.str());

    
        }

        os << indent3 << "}  // close HRR loop\n";

    }

    os << "\n";
    os << "\n";
}


void HRR_Writer::WriteHRRExternal_(std::ostream & os, const WriterBase & base) const
{
/*
    QAM finalam = base.FinalAM();

    if(hrrsteps_.first.size() > 0)
    {
        os << indent1 << "//////////////////////////////////////////////\n";
        os << indent1 << "// Contracted integrals: Horizontal recurrance\n";
        os << indent1 << "// Bra part\n";
        os << indent1 << "// Steps: " << hrrsteps_.first.size() << "\n";
        os << indent1 << "//////////////////////////////////////////////\n";
        os << "\n";
        os << "    #pragma simd\n";
        os << indent1 << "for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << indent1 << "{\n";


        os << "\n";

        for(const auto & it : hrrtopkets_)
        {
            // it.first is the AM for the ket part
            QAM thisam{finalam[0], finalam[1], it.first, 0};
            os << indent2 << "// form " << base.ArrVarName(thisam) << "\n";
            os << indent2 << "HRR_BRA_" << amchar[finalam[0]] << "_" << amchar[finalam[1]] << "(\n";

            // Pass the pointers
            for(const auto & itb : brahrr_ptrs_)
            {
                 os << "               "
                    << base.ArrVarName({itb[0], itb[1], it.first, 0}) 
                    << " + ( abcd * " << NCART(itb[0]) * NCART(itb[1]) * NCART(it.first) << " ),\n";
            }

            os << "               AB_x[abcd], AB_y[abcd], AB_z[abcd], " << NCART(it.first) << ");\n";
        }

        os << "\n";
        os << indent1 << "}\n";
        os << "\n";
        os << "\n";
    }

    if(hrrsteps_.second.size() > 0)
    {
        os << indent1 << "//////////////////////////////////////////////\n";
        os << indent1 << "// Contracted integrals: Horizontal recurrance\n";
        os << indent1 << "// Ket part\n";
        os << indent1 << "// Steps: " << hrrsteps_.second.size() << "\n";
        os << indent1 << "//////////////////////////////////////////////\n";
        os << "\n";

        DAM braam{finalam[0], finalam[1]};
        std::stringstream ss;
        ss << (NCART(braam[0]) * NCART(braam[1])); 

        os << "    #pragma simd\n";
        os << indent1 << "for(abcd = 0; abcd < nshell1234; ++abcd)\n";
        os << indent1 << "{\n";

        os << indent2 << "// form " << base.ArrVarName(finalam) << "\n";
        os << indent2 << "HRR_KET_" << amchar[finalam[2]] << "_" << amchar[finalam[3]] << "(\n";

        // Pass the pointes
        for(auto & it : kethrr_ptrs_)
        {
            os << "               "
               << base.ArrVarName({finalam[0], finalam[1], it[0], it[1]})
               << " + ( abcd * " << NCART(braam[0]) * NCART(braam[1]) * NCART(it[0]) * NCART(it[1]) << " ),\n"; 
        }
        os << "               CD_x[abcd], CD_y[abcd], CD_z[abcd], " << (NCART(braam[0]) * NCART(braam[1])) << ");\n";
        os << "\n";
        os << indent1 << "}\n";
        os << "\n";
        os << "\n";
    }
*/
os << "TODO\n";
}




void HRR_Writer::WriteHRR(std::ostream & os, const WriterBase & base) const
{
    if(base.GetOption(OPTION_INLINEHRR) > 0)
        WriteHRRInline_(os, base);
    else
        WriteHRRExternal_(os, base);
}



void HRR_Writer::WriteHRRFile(std::ostream & ofb, std::ostream & ofk, const WriterBase & base) const
{
    QAM finalam = base.FinalAM();


    if(hrrsteps_.first.size() > 0)
    {
        ofb << "//////////////////////////////////////////////\n";
        ofb << "// BRA: ( " << amchar[finalam[0]] << " " << amchar[finalam[1]] << " |\n";
        ofb << "// Steps: " << hrrsteps_.first.size() << "\n";
        ofb << "//////////////////////////////////////////////\n";
        ofb << "\n";

        // it.first is the AM for the ket part
        ofb << "#pragma omp declare simd simdlen(SIMD_LEN) uniform(ncart_ket)\n";
        ofb << "void HRR_BRA_";
        ofb << amchar[finalam[0]] << "_" << amchar[finalam[1]] << "(\n";

        // pointers to buffers (top bras)
        for(const auto & it : brakettopam_.first)
            ofb << indent5 << "double const * const restrict BRA_" << amchar[it[0]] << "_" << amchar[it[1]] << ",\n";

        // pointer to result buffer
        ofb << indent5 << "double * const restrict BRA_" << amchar[finalam[0]] << "_" << amchar[finalam[1]] << ",\n";

        ofb << indent5 << "const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket)\n";

        ofb << "{\n";
        ofb << indent1 << "int iket;\n";
        ofb << "\n";

        WriteBraSteps_(ofb, base, "ncart_ket", ""); 

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

        ofk << "#pragma omp declare simd simdlen(SIMD_LEN) uniform(ncart_bra)\n";
        ofk << "void HRR_KET_";
        ofk << amchar[finalam[2]] << "_" << amchar[finalam[3]] << "(\n";

        // pointers to buffers (top bras)
        for(const auto & it : brakettopam_.second)
            ofk << indent5 << "double const * const restrict KET_" << amchar[it[0]] << "_" << amchar[it[1]] << ",\n";

        // pointer to result buffer
        ofk << indent5 << "double * const restrict KET_" << amchar[finalam[2]] << "_" << amchar[finalam[3]] << ",\n";

        ofk << indent5 << "const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra)\n";
        ofk << "{\n";
        ofk << indent1 << "int ibra;\n";
        ofk << "\n";

        WriteKetSteps_(ofk, base, "ncart_bra", ""); 

        ofk << "\n";
        ofk << "}\n";
        ofk << "\n";
        ofk << "\n";
    }
}



void HRR_Writer::WriteHRRHeaderFile(std::ostream & os, const WriterBase & base) const
{
    QAM finalam = base.FinalAM();

    if(hrrsteps_.first.size() > 0)
    {
        os << indent1 << "//////////////////////////////////////////////\n";
        os << indent1 << "// BRA: ( " << amchar[finalam[0]] << " " << amchar[finalam[1]] << " |\n";
        os << indent1 << "// Steps: " << hrrsteps_.first.size() << "\n";
        os << indent1 << "//////////////////////////////////////////////\n";
        os << "\n";

        // it.first is the AM for the ket part
        os << "#pragma omp declare simd simdlen(SIMD_LEN) uniform(ncart_ket)\n";
        os << "void HRR_BRA_" << amchar[finalam[0]] << "_" << amchar[finalam[1]] << "(\n";
/*
        for(const auto & itb : brahrr_ptrs_)
            os << "                   double * const restrict BRA_" << amchar[itb[0]] << "_" << amchar[itb[1]] << ",\n";
*/
        os << indent5 << "const double bAB_x, const double bAB_y, const double bAB_z, const int ncart_ket);\n";
        os << "\n";
        os << "\n";
    }

    if(hrrsteps_.second.size() > 0)
    {
        os << indent1 << "//////////////////////////////////////////////\n";
        os << indent1 << "// KET: ( " << amchar[finalam[0]] << " " << amchar[finalam[1]] << " |\n";
        os << indent1 << "// Steps: " << hrrsteps_.second.size() << "\n";
        os << indent1 << "//////////////////////////////////////////////\n";
        os << "\n";

        os << "#pragma omp declare simd simdlen(SIMD_LEN) uniform(ncart_bra)\n";
        os << "void HRR_KET_" << amchar[finalam[2]] << "_" << amchar[finalam[3]] << "(\n";

        // Pass the pointes
/*
        for(auto & it : kethrr_ptrs_)
            os << "                  double * const restrict KET_" << amchar[it[0]] << "_" << amchar[it[1]] << ",\n";
*/

        os << indent5  << "const double kCD_x, const double kCD_y, const double kCD_z, const int ncart_bra)\n";
        os << "\n";
        os << "\n";
    }
}

