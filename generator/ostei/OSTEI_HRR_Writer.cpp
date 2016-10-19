#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"


////////////////////////
// Base HRR Writer
////////////////////////
OSTEI_HRR_Writer::OSTEI_HRR_Writer(const OSTEI_HRR_Algorithm_Base & hrr_algo, const OSTEI_GeneratorInfo & info, 
                                   int start_external, int start_general)
    : hrr_algo_(hrr_algo), info_(info), vinfo_(info.GetVectorInfo()),
      start_external_(start_external), start_general_(start_general)
{ }


bool OSTEI_HRR_Writer::HasHRR(void) const
{
    return hrr_algo_.HasHRR();
}

bool OSTEI_HRR_Writer::HasBraHRR(void) const
{
    return hrr_algo_.HasBraHRR();
}

bool OSTEI_HRR_Writer::HasKetHRR(void) const
{
    return hrr_algo_.HasKetHRR();
}

ConstantMap OSTEI_HRR_Writer::GetConstants(void) const
{
    return ConstantMap();
}


void OSTEI_HRR_Writer::WriteBraSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                                      const std::string & ncart_ket, const std::string & ketstr) const
{
    os << indent4 << "for(iket = 0; iket < " << ncart_ket << "; ++iket)\n";
    os << indent4 << "{\n";

    for(const auto & it : steps)
    {
        //os << std::string(20, ' ') << "// " << it << "\n";
    
        const char * sign = " + ";
        if(it.type == RRStepType::I) // moving from J->I
            sign = " - ";

        os << std::string(20, ' ');

        os << HRRBraStepVar_(it.target, ncart_ket, ketstr);

        os << " = ";
        os << HRRBraStepVar_(it.src[0], ncart_ket, ketstr);
        os << sign << "( hAB[" << static_cast<int>(it.xyz) << "] * ";
        os << HRRBraStepVar_(it.src[1], ncart_ket, ketstr);
        os << " );";
        os << "\n\n";
    }

    os << indent4 << "}\n";
    os << "\n";
}


void OSTEI_HRR_Writer::WriteKetSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                                      const std::string & ncart_bra, const std::string & brastr) const
{
        os << indent4 << "for(ibra = 0; ibra < " << ncart_bra << "; ++ibra)\n"; 
        os << indent4 << "{\n"; 
        for(const auto & it : steps)
        {
            //os << std::string(20, ' ') << "// " << it << "\n";

            const char * sign = " + ";
            if(it.type == RRStepType::K) // Moving from L->K
                sign = " - ";

            os << std::string(20, ' ');
    
            os << HRRKetStepVar_(it.target, brastr);

            os << " = ";
            os << HRRKetStepVar_(it.src[0], brastr);
            os << sign << "( hCD[" << static_cast<int>(it.xyz) << "] * ";
            os << HRRKetStepVar_(it.src[1], brastr);
            os << " );";
            os << "\n\n";
        }

        os << indent4 << "}\n"; 
}


std::string OSTEI_HRR_Writer::HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, 
                                           const std::string & ketstr) const
{
    std::string arrname = ArrVarName(d.left.am(), d.right.am(), ketstr);
    return StringBuilder("HRR_", arrname, "[", d.index(), " * ", ncart_ket, " + iket]"); 
}


std::string OSTEI_HRR_Writer::HRRKetStepVar_(const Doublet & d, const std::string & brastr) const
{
    std::string arrname = ArrVarName(brastr, d.left.am(), d.right.am());
    return StringBuilder("HRR_", arrname, "[ibra * ", d.ncart(), " + ", d.index(), "]"); 
}


void OSTEI_HRR_Writer::WriteHRR(std::ostream & os) const
{
    QAM finalam = info_.FinalAM();
    DAM finalbra{finalam[0], finalam[1]};

    os << indent3 << "//////////////////////////////////////////////\n";
    os << indent3 << "// Contracted integrals: Horizontal recurrance\n";
    os << indent3 << "//////////////////////////////////////////////\n";
    os << "\n";

    if(HasBraHRR())
    {
        os << "\n";
        os << indent4 << "const double hAB[3] = { P.AB_x[ab], P.AB_y[ab], P.AB_z[ab] };\n"; 
        os << "\n";
    }

    os << indent3 << "for(abcd = 0; abcd < nshellbatch; ++abcd, ++real_abcd)\n";
    os << indent3 << "{\n";

    os << "\n";
    os << indent4 << "// set up HRR pointers\n";
    for(const auto & it : hrr_algo_.TopAM())
        os << indent4 << "double const * restrict " << HRRVarName(it) << " = " << ArrVarName(it) << " + abcd * " << NCART(it) << ";\n";

    // and also for the final integral
    os << indent4 << "double * restrict " << HRRVarName(finalam) << " = " << ArrVarName(finalam) << " + real_abcd * " << NCART(finalam) << ";\n";
    os << "\n";
    
    if(HasBraHRR())
    {
        for(const auto & itk : hrr_algo_.TopKetAM()) // for all needed ket am
        for(const auto & itb : hrr_algo_.GetBraAMOrder()) // form these
        {
            QAM am = {itb[0], itb[1], itk[0], itk[1]};
            os << indent4 << "// form " << ArrVarName(am) << "\n";
            // allocate temporary if needed
            if(!info_.IsContQ(am) && !info_.IsFinalAM(am))
                os << indent4 << "double " << HRRVarName(am) << "[" << NCART(am) << "];\n";

            int L = am[0] + am[1];

            if(L < start_external_)
                WriteHRR_Bra_Inline_(os, am);
            else if(L < start_general_)
                WriteHRR_Bra_External_(os, am);
            else
                WriteHRR_Bra_General_(os, am);
            os << "\n";
        }
    }

    os << "\n";
    os << "\n";

    if(HasKetHRR())
    {
        os << "\n";
        os << indent4 << "const double hCD[3] = { Q.AB_x[cd+abcd], Q.AB_y[cd+abcd], Q.AB_z[cd+abcd] };\n";
        os << "\n";

        for(const auto & itk : hrr_algo_.GetKetAMOrder())
        {
            QAM am = {finalbra[0], finalbra[1], itk[0], itk[1]};
            os << indent4 << "// form " << ArrVarName(am) << "\n";
            // allocate temporary if needed
            if(!info_.IsContQ(am) && !info_.IsFinalAM(am))
                os << indent4 << "double " << HRRVarName(am) << "[" << NCART(am) << "];\n";

            int L = am[2] + am[3];

            if(L < start_external_)
                WriteHRR_Ket_Inline_(os, am);
            else if(L < start_general_)
                WriteHRR_Ket_External_(os, am);
            else
                WriteHRR_Ket_General_(os, am);
            os << "\n";
        }
    }


    os << "\n";
    os << indent3 << "}  // close HRR loop\n";
    os << "\n";
    os << "\n";
}

void OSTEI_HRR_Writer::WriteHRR_Bra_Inline_(std::ostream & os, QAM am) const
{
    // ncart_ket in string form
    std::string ncart_ket_str = StringBuilder(NCART(am[2], am[3]));

    // the ket part in string form
    std::string ket_str = StringBuilder(amchar[am[2]], "_", amchar[am[3]]);

    // actually write out the steps now
    WriteBraSteps_(os, hrr_algo_.GetBraSteps({am[0], am[1]}), ncart_ket_str, ket_str);
}

void OSTEI_HRR_Writer::WriteHRR_Ket_Inline_(std::ostream & os, QAM am) const
{
    // ncart_bra in string form
    std::string ncart_bra_str = StringBuilder(NCART(DAM{am[0], am[1]}));

    // the bra part in string form
    std::string bra_str = StringBuilder(amchar[am[0]], "_", amchar[am[1]]);

    WriteKetSteps_(os, hrr_algo_.GetKetSteps(DAM{am[2], am[3]}), ncart_bra_str, bra_str);
}

void OSTEI_HRR_Writer::WriteHRR_Bra_External_(std::ostream & os, QAM am) const
{
    // call function
    RRStepType rrstep = hrr_algo_.GetBraRRStep({am[0], am[1]});
    os << indent4 << "HRR_BRA_" << RRStepTypeToStr(rrstep) << "_"
       << amchar[am[0]] << "_" << amchar[am[1]] << "(\n";

    // pointer to result buffer
    os << indent5 << "" << HRRVarName(am) << ",\n";

    // pointer to requirements
    for(const auto & it : hrr_algo_.GetBraAMReq({am[0], am[1]}))
        os << indent5 << "" << HRRVarName({it[0], it[1], am[2], am[3]}) << ",\n";

    os << indent5 << "hAB, " << NCART(am[2], am[3]) << ");\n"; 
}

void OSTEI_HRR_Writer::WriteHRR_Ket_External_(std::ostream & os, QAM am) const
{
    RRStepType rrstep = hrr_algo_.GetKetRRStep({am[2], am[3]});
    os << indent4 << "HRR_KET_" << RRStepTypeToStr(rrstep) << "_"
       << amchar[am[2]] << "_" << amchar[am[3]] << "(\n";

    // pointer to result buffer
    os << indent5 << "" << HRRVarName(am) << ",\n";

    // pointer to requirements
    for(const auto & it : hrr_algo_.GetKetAMReq({am[2], am[3]}))
        os << indent5 << "" << HRRVarName({am[0], am[1], it[0], it[1]}) << ",\n";

    os << indent5 << "hCD, " << NCART(DAM{am[0], am[1]}) << ");\n";
}

void OSTEI_HRR_Writer::WriteHRR_Bra_General_(std::ostream & os, QAM am) const
{
    // call function
    RRStepType rrstep = hrr_algo_.GetBraRRStep({am[0], am[1]});
    os << indent4 << "ostei_general_hrr_" << RRStepTypeToStr(rrstep)
       << "(" << am[0] << ", " << am[1] << ", " << am[2] << ", " << am[3] << ", hAB, ";

    for(const auto & it : hrr_algo_.GenerateAMReq(am, rrstep))
        os << HRRVarName(it) << ", ";

    os << HRRVarName(am) << ");\n";
}

void OSTEI_HRR_Writer::WriteHRR_Ket_General_(std::ostream & os, QAM am) const
{
    // call function
    RRStepType rrstep = hrr_algo_.GetKetRRStep({am[2], am[3]});
    os << indent4 << "ostei_general_hrr_" << RRStepTypeToStr(rrstep)
       << "(" << am[0] << ", " << am[1] << ", " << am[2] << ", " << am[3] << ", hCD, ";

    for(const auto & it : hrr_algo_.GenerateAMReq(am, rrstep))
        os << HRRVarName(it) << ", ";

    os << HRRVarName(am) << ");\n";
}

void OSTEI_HRR_Writer::WriteHRRFile(std::ostream & of, std::ostream & ofh) const
{
    QAM finalam = info_.FinalAM();
    DAM braam = {finalam[0], finalam[1]};
    DAM ketam = {finalam[2], finalam[3]};

    // only call the function for the last doublets
    HRRDoubletStepList brasteps = hrr_algo_.GetBraSteps(braam);
    HRRDoubletStepList ketsteps = hrr_algo_.GetKetSteps(ketam);

    if(brasteps.size() > 0)
    {
        RRStepType brasteptype = hrr_algo_.GetBraRRStep(braam);
        std::string steptypestr = RRStepTypeToStr(brasteptype);

        of << "//////////////////////////////////////////////\n";
        of << "// BRA: ( " << amchar[braam[0]] << " " << amchar[braam[1]] << " |\n";
        of << "// Steps: " << brasteps.size() << "\n";
        of << "//////////////////////////////////////////////\n";
        of << "\n";

        std::stringstream prototype;

        prototype << "void HRR_BRA_" << steptypestr << "_";
        prototype << amchar[braam[0]] << "_" << amchar[braam[1]] << "(";

        // pointer to result buffer
        prototype << "double * const restrict " << HRRVarName(finalam[0], finalam[1], "X_X") << ",\n";

        // pointer to requirements
        for(const auto & it : hrr_algo_.GetBraAMReq(braam))
            prototype << indent5 << "double const * const restrict " << HRRVarName(it[0], it[1], "X_X") << ",\n";

        prototype << indent5 << "const double hAB[3], const int ncart_ket)";



        of << prototype.str() << "\n";
        of << "{\n";
        of << indent1 << "int iket;\n";
        of << "\n";

        WriteBraSteps_(of, brasteps, "ncart_ket", "X_X"); 

        of << "\n";
        of << "}\n";
        of << "\n";
        of << "\n";


        // header
        ofh << prototype.str() << ";\n\n";

    }

    if(ketsteps.size() > 0)
    {
        RRStepType ketsteptype = hrr_algo_.GetKetRRStep(ketam);
        std::string steptypestr = RRStepTypeToStr(ketsteptype);

        of << "//////////////////////////////////////////////\n";
        of << "// KET: | " << amchar[ketam[0]] << " " << amchar[ketam[1]] << " )\n";
        of << "// Steps: " << ketsteps.size() << "\n";
        of << "//////////////////////////////////////////////\n";
        of << "\n";

        std::stringstream prototype;

        prototype << "void HRR_KET_" << steptypestr << "_";
        prototype << amchar[ketam[0]] << "_" << amchar[ketam[1]] << "(";

        // pointer to result buffer
        prototype << "double * const restrict " << HRRVarName("X_X", finalam[2], finalam[3]) << ",\n";

        // pointer to requirements
        for(const auto & it : hrr_algo_.GetKetAMReq(ketam))
            prototype << indent5 << "double const * const restrict " << HRRVarName("X_X", it[0], it[1]) << ",\n";

        prototype << indent5 << "const double hCD[3], const int ncart_bra)";



        of << prototype.str() << "\n";
        of << "{\n";
        of << indent1 << "int ibra;\n";
        of << "\n";

        WriteKetSteps_(of, ketsteps, "ncart_bra", "X_X"); 

        of << "\n";
        of << "}\n";
        of << "\n";
        of << "\n";


        // header
        ofh << prototype.str() << ";\n\n";
    }
}


