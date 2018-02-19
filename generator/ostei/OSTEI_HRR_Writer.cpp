#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_HRR_Writer.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"


////////////////////////
// Base HRR Writer
////////////////////////
OSTEI_HRR_Writer::OSTEI_HRR_Writer(const OSTEI_HRR_Algorithm_Base & hrr_algo, const OSTEI_GeneratorInfo & info, 
                                   int start_external, int start_general)
    : hrr_algo_(hrr_algo), info_(info),
      start_external_(start_external), start_general_(start_general)
{ }


const OSTEI_HRR_Algorithm_Base & OSTEI_HRR_Writer::Algo(void) const
{
    return hrr_algo_;
}

void OSTEI_HRR_Writer::WriteBraSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                                      const std::string & tag,
                                      const std::string & ncart_ket, const std::string & ketstr) const
{
    os << indent4 << "for(iket = 0; iket < " << ncart_ket << "; ++iket)\n";
    os << indent4 << "{\n";

    for(const auto & it : steps)
    {
        //os << std::string(20, ' ') << "// " << it << "\n";

        // add the appropriate integral tags
        Doublet target(it.target);
        Doublet src0(it.src[0]);
        Doublet src1(it.src[1]);
        target.tag = tag;
        src0.tag = tag;
        src1.tag = tag;
    
        const char * sign = " + ";
        if(it.type == RRStepType::I) // moving from J->I
            sign = " - ";

        os << std::string(20, ' ');

        os << HRRBraStepVar_(target, ncart_ket, ketstr);

        os << " = ";
        os << HRRBraStepVar_(src0, ncart_ket, ketstr);
        os << sign << "( hAB[" << static_cast<int>(it.xyz) << "] * ";
        os << HRRBraStepVar_(src1, ncart_ket, ketstr);
        os << " );";
        os << "\n\n";
    }

    os << indent4 << "}\n";
    os << "\n";
}


void OSTEI_HRR_Writer::WriteKetSteps_(std::ostream & os, const HRRDoubletStepList & steps,
                                      const std::string & tag,
                                      const std::string & ncart_bra, const std::string & brastr) const
{
    //if(info_.Vectorized())
    //    os << indent4 << "#pragma omp simd simdlen(SIMINT_SIMD_LEN)\n";

    os << indent4 << "for(ibra = 0; ibra < " << ncart_bra << "; ++ibra)\n"; 
    os << indent4 << "{\n"; 
    for(const auto & it : steps)
    {
        //os << std::string(20, ' ') << "// " << it << "\n";

        // add the appropriate integral tags
        Doublet target(it.target);
        Doublet src0(it.src[0]);
        Doublet src1(it.src[1]);
        target.tag = tag;
        src0.tag = tag;
        src1.tag = tag;

        const char * sign = " + ";
        if(it.type == RRStepType::K) // Moving from L->K
            sign = " - ";

        os << std::string(20, ' ');
    
        os << HRRKetStepVar_(target, brastr);

        os << " = ";

        os << HRRKetStepVar_(src0, brastr);
        os << sign << "( hCD[" << static_cast<int>(it.xyz) << "] * ";
        os << HRRKetStepVar_(src1, brastr);
        os << " );";
        os << "\n\n";
    }

    os << indent4 << "}\n"; 
}


std::string OSTEI_HRR_Writer::HRRBraStepVar_(const Doublet & d, const std::string & ncart_ket, 
                                             const std::string & ketstr) const
{
    std::string arrname = ArrVarName(d, ketstr);
    return StringBuilder("HRR_", arrname, "[", d.index(), " * ", ncart_ket, " + iket]"); 
}


std::string OSTEI_HRR_Writer::HRRKetStepVar_(const Doublet & d, const std::string & brastr) const
{
    std::string arrname = ArrVarName(brastr, d);
    return StringBuilder("HRR_", arrname, "[ibra * ", d.ncart(), " + ", d.index(), "]"); 
}


void OSTEI_HRR_Writer::WriteHRR(std::ostream & os) const
{
    QAM finalam = info_.FinalAM();

    os << indent3 << "//////////////////////////////////////////////\n";
    os << indent3 << "// Contracted integrals: Horizontal recurrance\n";
    os << indent3 << "//////////////////////////////////////////////\n";
    os << "\n\n";

    if(hrr_algo_.HasBraHRR())
        os << indent3 << "const double hAB[3] = { P.AB_x[ab], P.AB_y[ab], P.AB_z[ab] };\n"; 
    os << "\n\n";

    os << indent3 << "for(abcd = 0; abcd < nshellbatch; ++abcd, ++real_abcd)\n";
    os << indent3 << "{\n";

    if(hrr_algo_.HasKetHRR())
        os << indent4 << "const double hCD[3] = { Q.AB_x[cd+abcd], Q.AB_y[cd+abcd], Q.AB_z[cd+abcd] };\n";
    os << "\n";
    os << indent4 << "// set up HRR pointers\n";
    for(const auto & it : hrr_algo_.TopAM())
        os << indent4 << "double const * restrict " << HRRVarName(it) << " = " << ArrVarName(it) << " + abcd * " << NCART(it) << ";\n";

    // and also for the final integral
    if(info_.Deriv() == 0)
        os << indent4 << "double * restrict " << HRRVarName(finalam) << " = " << ArrVarName(finalam) << " + real_abcd * " << NCART(finalam) << ";\n";
    else if(info_.Deriv() == 1)
        os << indent4 << "double * restrict " << HRRVarName(finalam) << " = " << ArrVarName(finalam) << " + real_abcd * " << NCART(finalam) << " * 12;\n";
    os << "\n";

    for(auto am : hrr_algo_.GetAMOrder())
    {
        DoubletType steptype = hrr_algo_.GetDoubletStep(am);

        if(steptype == DoubletType::BRA)
        {
            os << indent4 << "// form " << ArrVarName(am) << "\n";
            int L = am[0] + am[1];

            if(L < start_external_)
                WriteHRR_Bra_Inline_(os, am);
            else if(L < start_general_)
                WriteHRR_Bra_External_(os, am);
            else
                WriteHRR_Bra_General_(os, am);
            os << "\n";
        }
        else
        {
            os << indent4 << "// form " << ArrVarName(am) << "\n";
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

}

void OSTEI_HRR_Writer::WriteHRR_Bra_Inline_(std::ostream & os, QAM am) const
{
    // ncart_ket in string form
    std::string ncart_ket_str = StringBuilder(NCART(am[2], am[3]));

    // the ket part in string form
    std::string ket_str = StringBuilder(amchar[am[2]], "_", amchar[am[3]]);

    // actually write out the steps now
    auto brasteps = hrr_algo_.GetBraSteps(DAM{am[0], am[1]});
    WriteBraSteps_(os, brasteps, am.tag, ncart_ket_str, ket_str);
}

void OSTEI_HRR_Writer::WriteHRR_Ket_Inline_(std::ostream & os, QAM am) const
{
    // ncart_bra in string form
    std::string ncart_bra_str = StringBuilder(NCART(DAM{am[0], am[1]}));

    // the bra part in string form
    std::string bra_str = StringBuilder(amchar[am[0]], "_", amchar[am[1]]);

    auto ketsteps = hrr_algo_.GetKetSteps(DAM{am[2], am[3]});
    WriteKetSteps_(os, ketsteps, am.tag, ncart_bra_str, bra_str);
}

void OSTEI_HRR_Writer::WriteHRR_Bra_External_(std::ostream & os, QAM am) const
{
    // call function
    RRStepType rrstep = hrr_algo_.GetBraRRStep({am[0], am[1]});
    os << indent4 << "HRR_" << RRStepTypeToStr(rrstep) << "_"
       << amchar[am[0]] << "_" << amchar[am[1]] << "(\n";

    // pointer to result buffer
    os << indent5 << "" << HRRVarName(am) << ",\n";

    // pointer to requirements
    for(const auto & it : hrr_algo_.GetBraAMReq({am[0], am[1]}))
        os << indent5 << "" << HRRVarName({it[0], it[1], am[2], am[3], am.tag}) << ",\n";

    os << indent5 << "hAB, " << NCART(am[2], am[3]) << ");\n"; 

    // Mark this as required in the log file
    std::cout << "SIMINT EXTERNAL HRR " << RRStepTypeToStr(rrstep)
              << " " << am[0] << " " << am[1] << "\n";
}

void OSTEI_HRR_Writer::WriteHRR_Ket_External_(std::ostream & os, QAM am) const
{
    RRStepType rrstep = hrr_algo_.GetKetRRStep({am[2], am[3]});
    os << indent4 << "HRR_" << RRStepTypeToStr(rrstep) << "_"
       << amchar[am[2]] << "_" << amchar[am[3]] << "(\n";

    // pointer to result buffer
    os << indent5 << "" << HRRVarName(am) << ",\n";

    // pointer to requirements
    for(const auto & it : hrr_algo_.GetKetAMReq({am[2], am[3]}))
        os << indent5 << "" << HRRVarName({am[0], am[1], it[0], it[1], am.tag}) << ",\n";

    os << indent5 << "hCD, " << NCART(DAM{am[0], am[1]}) << ");\n";

    // Mark this as required in the log file
    std::cout << "SIMINT EXTERNAL HRR " << RRStepTypeToStr(rrstep)
              << " " << am[2] << " " << am[3] << "\n";
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
        

        prototype << "void HRR_" << steptypestr << "_";
        prototype << amchar[braam[0]] << "_" << amchar[braam[1]] << "(";

        // pointer to result buffer
        prototype << "double * const restrict " << HRRVarName(finalam[0], finalam[1], "X_X") << ",\n";

        // pointer to requirements
        for(const auto & it : hrr_algo_.GetBraAMReq(braam))
            prototype << indent5 << "double const * const restrict " << HRRVarName(it[0], it[1], "X_X") << ",\n";

        prototype << indent5 << "const double hAB[3], const int ncart_ket)";

        // Create different functions for different ncart_ket values
        const int ncart_ket_vals[15] = {1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120};
        for (int ival = 0; ival < 15; ival++)
        {
            std::stringstream prototype_ncartket;
            prototype_ncartket << "void HRR_" << steptypestr << "_";
            prototype_ncartket << amchar[braam[0]] << "_" << amchar[braam[1]] << "_" << ncart_ket_vals[ival] << "(";

            // pointer to result buffer
            prototype_ncartket << "double * const restrict " << HRRVarName(finalam[0], finalam[1], "X_X") << ",\n";

            // pointer to requirements
            for(const auto & it : hrr_algo_.GetBraAMReq(braam))
                prototype_ncartket << indent5 << "double const * const restrict " << HRRVarName(it[0], it[1], "X_X") << ",\n";

            prototype_ncartket << indent5 << "const double hAB[3])";
            
            of << prototype_ncartket.str() << "\n";
            of << "{\n";
            of << indent1 << "int iket;\n";
            of << "\n";
            
            WriteBraSteps_(of, brasteps, "", std::to_string(ncart_ket_vals[ival]), "X_X"); 
            
            of << "\n";
            of << "}\n";
            of << "\n";
            of << "\n";
        }
        
        
        of << prototype.str() << "\n{\n";  // Original function header
        
        // Wrapper for different ncart_ket value functions
        of << indent1 << "int done_brastep = 0;\n";
        for (int ival = 0; ival < 15; ival++)
        {
            of << indent1 << "if (ncart_ket == " << ncart_ket_vals[ival] << ")\n";
            of << indent1 << "{\n";
            
            // function name
            of << indent1 << indent1 << "HRR_" << steptypestr << "_" << amchar[braam[0]] << "_";
            of << amchar[braam[1]] << "_" << ncart_ket_vals[ival] << "(";
            
            // pointer to result buffer
            of << HRRVarName(finalam[0], finalam[1], "X_X") << ", ";
            
            // pointer to requirements
            for(const auto & it : hrr_algo_.GetBraAMReq(braam))
                of << HRRVarName(it[0], it[1], "X_X") << ", ";
            of << "hAB);\n";
            
            of << indent1 << indent1 << "done_brastep = 1;\n";
            
            of << indent1 << "}\n\n";
        }
        
        
        // General situtaion
        of << indent1 << "if (done_brastep == 0) \n";
        of << indent1 << "{\n";
        
        of << indent1 << indent1 << "int iket;\n";
        of << "\n";
        
        WriteBraSteps_(of, brasteps, "", "ncart_ket", "X_X"); 
        
        of << indent1 << "}\n";

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

        prototype << "void HRR_" << steptypestr << "_";
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

        WriteKetSteps_(of, ketsteps, "", "ncart_bra", "X_X"); 

        of << "\n";
        of << "}\n";
        of << "\n";
        of << "\n";


        // header
        ofh << prototype.str() << ";\n\n";
    }
}


