#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"

#include <algorithm> // std::count


static bool IsPointer(const std::string & var)
{
    static const std::set<std::string> ptrvars{
                "P_PA", "P_PB", "Q_PA", "Q_PB", "PQ", "Pxyz",
                "aop_PQ", "aoq_PQ" };
    return ptrvars.count(var);
}


///////////////////////////////
// Base VRR Writer
///////////////////////////////
OSTEI_VRR_Writer::OSTEI_VRR_Writer(const OSTEI_VRR_Algorithm_Base & vrr_algo, const OSTEI_GeneratorInfo & info,
                                   int start_external, int start_general)
    : vrr_algo_(vrr_algo), info_(info),
      start_external_(start_external), start_general_(start_general)
{ 
}

const OSTEI_VRR_Algorithm_Base & OSTEI_VRR_Writer::Algo(void) const
{
    return vrr_algo_;
}

void OSTEI_VRR_Writer::WriteVRRSteps_(std::ostream & os, QAM am, const VRR_StepSet & vs, const std::string & num_n) const
{
    os << "\n";
    os << indent5 << "// Forming " << PrimVarName(am) << "[" << num_n << " * " << NCART(am) << "];\n";

    os << indent5 << "for(n = 0; n < " << num_n << "; ++n)  // loop over orders of auxiliary function\n";
    os << indent5 << "{\n";

    os << "\n";

    // iterate over the requirements
    // should be in order since it's a set
    for(const auto & it : vs)
    {
        // Get the stepping
        XYZStep step = it.xyz;
        std::string stepdir = StringBuilder("[", static_cast<int>(step), "]");

        std::string primname = StringBuilder(PrimVarName(am), "[n * ", NCART(am), " + ", it.target.index(), "]");
        std::string srcname[8];

        if(it.src[0])
            srcname[0] = StringBuilder(PrimVarName(it.src[0].amlist()), "[n * ", it.src[0].ncart(), " + ", it.src[0].index(), "]");
        if(it.src[1])
            srcname[1] = StringBuilder( PrimVarName(it.src[1].amlist()), "[(n+1) * ", it.src[1].ncart(), " + ", it.src[1].index(), "]");
        if(it.src[2])
            srcname[2] = StringBuilder(PrimVarName(it.src[2].amlist()), "[n * ", it.src[2].ncart(), " + ", it.src[2].index(), "]");
        if(it.src[3])
            srcname[3] = StringBuilder(PrimVarName(it.src[3].amlist()), "[(n+1) * ", it.src[3].ncart(), " + ", it.src[3].index(), "]");
        if(it.src[4])
            srcname[4] = StringBuilder(PrimVarName(it.src[4].amlist()), "[n * ", it.src[4].ncart(), " + ", it.src[4].index(), "]");
        if(it.src[5])
            srcname[5] = StringBuilder(PrimVarName(it.src[5].amlist()), "[(n+1) * ", it.src[5].ncart(), " + ", it.src[5].index(), "]");
        if(it.src[6])
            srcname[6] = StringBuilder(PrimVarName(it.src[6].amlist()), "[(n+1) * ", it.src[6].ncart(), " + ", it.src[6].index(), "]");
        if(it.src[7])
            srcname[7] = StringBuilder(PrimVarName(it.src[7].amlist()), "[(n+1) * ", it.src[7].ncart(), " + ", it.src[7].index(), "]");

        std::string aoppq, aover;
        std::string vrr_const0, vrr_const1, vrr_const2, vrr_const3;

        // Note - the signs on a_over_p, a_over_q, etc, are taken care of in FileWriter.cpp
        if(it.type == RRStepType::I || it.type == RRStepType::J)
        {
            aoppq = StringBuilder("aop_PQ", stepdir);
            aover = "a_over_p";
            vrr_const0 = StringBuilder("vrr_const_", it.ijkl[0], "_over_2p");
            vrr_const1 = StringBuilder("vrr_const_", it.ijkl[1], "_over_2p");
            vrr_const2 = StringBuilder("vrr_const_", it.ijkl[2], "_over_2pq");
            vrr_const3 = StringBuilder("vrr_const_", it.ijkl[3], "_over_2pq");
        }
        else
        {
            aoppq = StringBuilder("aoq_PQ", stepdir);
            aover = "a_over_q";
            vrr_const0 = StringBuilder("vrr_const_", it.ijkl[2], "_over_2q");
            vrr_const1 = StringBuilder("vrr_const_", it.ijkl[3], "_over_2q");
            vrr_const2 = StringBuilder("vrr_const_", it.ijkl[0], "_over_2pq");
            vrr_const3 = StringBuilder("vrr_const_", it.ijkl[1], "_over_2pq");
        }


        if(it.type == RRStepType::I || it.type == RRStepType::J)
        {
            std::string ppa = StringBuilder("P_PA", stepdir);

            if(it.type != RRStepType::I)
                ppa = StringBuilder("P_PB", stepdir);

            os << indent6 << primname << " = SIMINT_MUL(" << ppa << ", " << srcname[0] << ");\n";
        }
        else
        {
            std::string qpa = StringBuilder("Q_PA", stepdir);

            if(it.type != RRStepType::K)
                qpa = StringBuilder("Q_PB", stepdir);

            os << indent6 << primname << " = SIMINT_MUL(" << qpa << ", " << srcname[0] << ");\n";
        }

        if(it.src[1])
            os << indent6 << primname << " = SIMINT_FMADD( " << aoppq << ", " << srcname[1] << ", " << primname << ");\n";
        if(it.src[2] && it.src[3])
            os << indent6 << primname << " = SIMINT_FMADD( " << vrr_const0 << ", SIMINT_FMADD(" << aover << ", " << srcname[3] << ", " << srcname[2] << "), " << primname << ");\n";
        if(it.src[4] && it.src[5])
            os << indent6 << primname << " = SIMINT_FMADD( " << vrr_const1 << ", SIMINT_FMADD(" << aover << ", " << srcname[5] << ", " << srcname[4] << "), " << primname << ");\n";
        if(it.src[6])
            os << indent6 << primname << " = SIMINT_FMADD( " << vrr_const2 << ", " << srcname[6] << ", " << primname << ");\n";
        if(it.src[7])
            os << indent6 << primname << " = SIMINT_FMADD( " << vrr_const3 << ", " << srcname[7] << ", " << primname << ");\n";
        
        os << "\n";

    }

    os << indent5 << "}\n";
}


ConstantMap OSTEI_VRR_Writer::GetConstants(void) const
{
    ConstantMap cm;
    for(int i = 1; i <= vrr_algo_.GetMaxInt(); i++)
        cm.emplace(StringBuilder("const_", i), StringBuilder(i));
    return cm;
}


void OSTEI_VRR_Writer::WriteVRR(std::ostream & os) const
{
    os << "\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << indent5 << "// Primitive integrals: Vertical recurrance\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << "\n";

    // constants
    for(const auto & it : vrr_algo_.GetAllInt_2p())
    {
        if(it == 1)
            os << indent5 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2p = one_over_2p;\n"; 
        else
            os << indent5 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2p = SIMINT_MUL(const_" << it << ", one_over_2p);\n";
    }

    for(const auto & it : vrr_algo_.GetAllInt_2q())
    {
        if(it == 1)
            os << indent5 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2q = one_over_2q;\n"; 
        else
            os << indent5 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2q = SIMINT_MUL(const_" << it << ", one_over_2q);\n";
    }

    for(const auto & it : vrr_algo_.GetAllInt_2pq())
    {
        if(it == 1)
            os << indent5 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2pq = one_over_2pq;\n"; 
        else
            os << indent5 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2pq = SIMINT_MUL(const_" << it << ", one_over_2pq);\n";
    }

    os << "\n\n";

    // actually write out
    for(const auto & am : vrr_algo_.GetAMOrder())
    {
        if(am == QAM{0,0,0,0})
            continue;

        int L = am[0] + am[1] + am[2] + am[3];
        if(L < start_external_)
            WriteVRR_Inline_(os, am);
        else if(L < start_general_)
            WriteVRR_External_(os, am);
        else
            WriteVRR_General_(os, am);

        os << "\n\n";
    }
    
}


void OSTEI_VRR_Writer::WriteVRR_Inline_(std::ostream & os, QAM am) const
{
    WriteVRRSteps_(os, am, vrr_algo_.GetSteps(am), 
                   std::to_string(vrr_algo_.GetMReq(am)+1));
}


void OSTEI_VRR_Writer::WriteVRR_External_(std::ostream & os, QAM am) const
{

    // iterate over increasing am
    RRStepType rrstep = vrr_algo_.GetRRStep(am);

    os << indent5 << "VRR_" << RRStepTypeToStr(rrstep) << "_" 
                  << amchar[am[0]] << "_" << amchar[am[1]] << "_"
                  << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

    os << indent7 << PrimVarName(am) << ",\n";

    for(const auto & it : vrr_algo_.GetAMReq(am))
        os << indent7 << PrimVarName(it) << ",\n";
    
    for(const auto & it : vrr_algo_.GetVarReq(am))
        os << indent7 << it << ",\n";

    os << indent7 << (vrr_algo_.GetMReq(am)+1) << ");\n";

    // Mark this as required in the log file
    std::cout << "SIMINT EXTERNAL VRR " << RRStepTypeToStr(rrstep)
              << " " << am[0] << " " << am[1] << " " << am[2] << " " << am[3] << "\n";
}

void OSTEI_VRR_Writer::WriteVRR_General_(std::ostream & os, QAM am) const
{
    // is this a special case (incrementing with all other AM == 0)?
    RRStepType rrstep = vrr_algo_.GetRRStep(am);
    std::string stepstr = RRStepTypeToStr(rrstep);
    auto varreq = vrr_algo_.GetVarReq(am);

    auto amreq = vrr_algo_.GenerateAMReq(am, rrstep, true);
    std::vector<std::string> amreq_str;
    for(const auto &it : amreq)
    {
        if(!ValidQAM(it))
            amreq_str.push_back("NULL");
        else
            amreq_str.push_back(PrimVarName(it));
    }
    

    if(std::count(am.begin(), am.end(), 0) == 3)
    {
        os << indent5 << "ostei_general_vrr1_" << stepstr << "("
                      << am[static_cast<int>(rrstep)] << ", " << (vrr_algo_.GetMReq(am)+1) << ",\n"
           << indent7;

        if(rrstep == RRStepType::I || rrstep == RRStepType::J)
            os << "one_over_2p, a_over_p, aop_PQ, ";
        else
            os << "one_over_2q, a_over_q, aoq_PQ, ";

        switch(rrstep)
        {
            case RRStepType::I:
                os << "P_PA,\n";
                os << indent7 << amreq_str[0] << ", "
                              << amreq_str[1] << ", ";
                break;
        
            case RRStepType::J:
                os << "P_PB,\n";
                os << indent7 << amreq_str[0] << ", "
                              << amreq_str[2] << ", ";
                break;

            case RRStepType::K:
                os << "Q_PA,\n";
                os << indent7 << amreq_str[0] << ", "
                              << amreq_str[1] << ", ";
                break;

            case RRStepType::L:
                os << "Q_PB,\n";
                os << indent7 << amreq_str[0] << ", "
                              << amreq_str[2] << ", ";
                break;
        }

        os << PrimVarName(am) << ");\n"; // output array

    }
    else
    {
        os << indent5 << "ostei_general_vrr_" << stepstr << "("
                      << am[0] << ", " << am[1] << ", " << am[2] << ", " << am[3] << ", " << (vrr_algo_.GetMReq(am)+1) << ",\n"
           << indent7;

        if(rrstep == RRStepType::I || rrstep == RRStepType::J)
            os << "one_over_2p, a_over_p, one_over_2pq, aop_PQ, ";
        else
            os << "one_over_2q, a_over_q, one_over_2pq, aoq_PQ, ";


        if(rrstep == RRStepType::I)
            os << "P_PA,\n";
        if(rrstep == RRStepType::J)
            os << "P_PB,\n";
        if(rrstep == RRStepType::K)
            os << "Q_PA,\n";
        if(rrstep == RRStepType::L)
            os << "Q_PB,\n";

        os << indent7;
        for(const auto &it : amreq_str)
            os << it << ", ";

        os << PrimVarName(am) << ");\n";
    }
}


void OSTEI_VRR_Writer::WriteVRRFile(std::ostream & os, std::ostream & osh) const
{
    QAM am = info_.FinalAM();

    os << "//////////////////////////////////////////////\n";
    os << "// VRR: ( " << amchar[am[0]] << " " << amchar[am[1]] << " | " << amchar[am[2]] << " " << amchar[am[3]] << " )\n";
    os << "//////////////////////////////////////////////\n";

    RRStepType rrstep = vrr_algo_.GetRRStep(am);

    std::stringstream prototype;

    std::string stepchar = RRStepTypeToStr(rrstep);

    prototype << "void VRR_" << stepchar << "_"
              << amchar[am[0]] << "_" << amchar[am[1]] << "_"
              << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

    // final target
    prototype << indent3 << "SIMINT_DBLTYPE * const restrict " << PrimVarName(am) << ",\n";

    for(const auto & it : vrr_algo_.GetAMReq(am))
        prototype << indent3 << "const SIMINT_DBLTYPE * const restrict " << PrimVarName(it) << ",\n";
    
    for(const auto & it : vrr_algo_.GetVarReq(am))
    {
        if(IsPointer(it))
            prototype << indent3 << "const SIMINT_DBLTYPE * " << it << ",\n";
        else
            prototype << indent3 << "const SIMINT_DBLTYPE " << it << ",\n";
    }

    prototype << indent3 << "const int num_n)";

    os << prototype.str() << "\n";
    os << "{\n";


    // is this in standard order, or one of the permutations?
    if(am[1] == 0 && am[3] == 0)
    {
        // standard ( X s | Y s )
        // Make this one

        os << "    int n = 0;\n";

        for(const auto & it : vrr_algo_.GetIntReq_2p(am))
            os << indent1 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2p = SIMINT_MUL(SIMINT_DBLSET1(" << it << "), one_over_2p);\n"; 

        for(const auto & it : vrr_algo_.GetIntReq_2q(am))
            os << indent1 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2q = SIMINT_MUL(SIMINT_DBLSET1(" << it << "), one_over_2q);\n"; 

        for(const auto & it : vrr_algo_.GetIntReq_2pq(am))
            os << indent1 << "const SIMINT_DBLTYPE vrr_const_" << it << "_over_2pq = SIMINT_MUL(SIMINT_DBLSET1(" << it << "), one_over_2pq);\n"; 

        // Write out the steps
        WriteVRRSteps_(os, am, vrr_algo_.GetSteps(am), "num_n");

    }
    else
    {
        // permutation. Call one of the others
        QAM tocall = am;
        if(am[0] == 0)
            std::swap(tocall[0], tocall[1]);
        if(am[2] == 0)
            std::swap(tocall[2], tocall[3]);

        if(rrstep == RRStepType::J)
            stepchar = 'I';
        else if(rrstep == RRStepType::L)
            stepchar = 'K';

        os << indent1 << "// Routines are identical except for swapping of\n";
        os << indent1 << "// PA with PB and QC with QD\n";
        os << "\n";
        os << indent1 << "VRR_" << stepchar << "_"
                      << amchar[tocall[0]] << "_" << amchar[tocall[1]] << "_"
                      << amchar[tocall[2]] << "_" << amchar[tocall[3]]  << "(\n";
    
        // final target
        os << indent3 << PrimVarName(am) << ",\n";
        
        for(const auto & it : vrr_algo_.GetAMReq(am))
            os << indent3 << PrimVarName(it) << ",\n";
    
        for(const auto & it : vrr_algo_.GetVarReq(am))
            os << indent3 << it << ",\n";

        os << indent3 << "num_n);";
    }

    os << "\n";
    os << "}\n";
    os << "\n\n";

    // header
    osh << prototype.str() << ";\n\n";
}

