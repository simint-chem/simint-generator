#include "generator/ostei/OSTEI_GeneratorInfo.hpp"
#include "generator/ostei/OSTEI_VRR_Writer.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"


// map from RRStepType to a character, for creating function names
static const std::map<RRStepType, char> rrstep_char{ {RRStepType::I, 'I'},
                                                     {RRStepType::J, 'J'},
                                                     {RRStepType::K, 'K'},
                                                     {RRStepType::L, 'L'} };

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
OSTEI_VRR_Writer::OSTEI_VRR_Writer(const OSTEI_VRR_Algorithm_Base & vrr_algo, const OSTEI_GeneratorInfo & info)
    : vrr_algo_(vrr_algo), info_(info), vinfo_(info_.GetVectorInfo())
{ 
}

bool OSTEI_VRR_Writer::HasVRR(void) const
{
    return vrr_algo_.HasVRR();
}

bool OSTEI_VRR_Writer::HasBraVRR(void) const
{
    return vrr_algo_.HasBraVRR();
}

bool OSTEI_VRR_Writer::HasKetVRR(void) const
{
    return vrr_algo_.HasKetVRR();
}

bool OSTEI_VRR_Writer::HasVRR_I(void) const
{
    return vrr_algo_.HasVRR_I();
}

bool OSTEI_VRR_Writer::HasVRR_J(void) const
{
    return vrr_algo_.HasVRR_J();
}

bool OSTEI_VRR_Writer::HasVRR_K(void) const
{
    return vrr_algo_.HasVRR_K();
}

bool OSTEI_VRR_Writer::HasVRR_L(void) const
{
    return vrr_algo_.HasVRR_L();
}

ConstantMap OSTEI_VRR_Writer::GetConstants(void) const
{
    // by default, return empty
    return ConstantMap();
}

void OSTEI_VRR_Writer::DeclarePrimArrays(std::ostream & os) const
{
    os << indent5 << "// Holds the auxiliary integrals ( i 0 | k 0 )^m in the primitive basis\n";
    os << indent5 << "// with m as the slowest index\n";


    for(const auto & qam : vrr_algo_.GetAllAM()) 
    {
        // add +1 fromm required m values to account for 0
        os << indent5 << "// AM = (" << qam[0] << " " << qam[1] << " | " << qam[2] << " " << qam[3] << " )\n";
        os << indent5 << vinfo_.DoubleType() << " " << PrimVarName(qam)
           << "[" << (vrr_algo_.GetMReq(qam)+1) << " * " 
           << NCART(qam) << "] SIMINT_ALIGN_ARRAY_DBL;\n";
    }

    os << "\n\n";
}

void OSTEI_VRR_Writer::DeclarePrimPointers(std::ostream & os) const
{
    for(const auto & qam : vrr_algo_.GetAllAM()) 
    {
        if(info_.IsContQ(qam))
            os << indent4 << "double * restrict " << PrimPtrName(qam)
               << " = " << ArrVarName(qam) << " + abcd * " << NCART(qam) << ";\n";
    }

    os << "\n\n";
}

void OSTEI_VRR_Writer::WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepSet & vs, const std::string & num_n) const
{
    os << "\n";
    os << indent5 << "// Forming " << PrimVarName(qam) << "[" << num_n << " * " << NCART(qam) << "];\n";

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

        std::string primname = StringBuilder(PrimVarName(qam), "[n * ", NCART(qam), " + ", it.target.index(), "]");
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


        if(info_.HasFMA())
        {
            if(it.type == RRStepType::I || it.type == RRStepType::J)
            {
                os << indent6 << primname << " = ";
                if(it.type == RRStepType::I)
                    os << "P_PA" << stepdir; 
                else
                    os << "P_PB" << stepdir;

                os << " * " << srcname[0] << ";\n";
            }
            else
            {
                os << indent6 << primname;

                if(it.type == RRStepType::K)
                    os << " = Q_PA" << stepdir;
                else
                    os << " = Q_PB" << stepdir;

                os << " * " << srcname[0] << ";\n";
            }

            if(it.src[1])
                os << indent6 << primname << " = " << vinfo_.FMAdd(aoppq, srcname[1], primname) << ";\n"; 
            if(it.src[2] && it.src[3])
                os << indent6 << primname << " = " << vinfo_.FMAdd(vrr_const0, vinfo_.FMAdd(aover, srcname[3], srcname[2]), primname) << ";\n"; 
            if(it.src[4] && it.src[5])
                os << indent6 << primname << " = " << vinfo_.FMAdd(vrr_const1, vinfo_.FMAdd(aover, srcname[5], srcname[4]), primname) << ";\n"; 
            if(it.src[6])
                os << indent6 << primname << " = " << vinfo_.FMAdd(vrr_const2, srcname[6], primname) << ";\n";
            if(it.src[7])
                os << indent6 << primname << " = " << vinfo_.FMAdd(vrr_const3, srcname[7], primname) << ";\n";
        }
        else
        {
            if(it.type == RRStepType::I || it.type == RRStepType::J)
            {
                os << indent6 << primname;

                if(it.type == RRStepType::I)
                    os << " = P_PA" << stepdir;
                else                            
                    os << " = P_PB" << stepdir;

                os << " * " << srcname[0];
            }
            else
            {
                os << indent6 << primname;

                if(it.type == RRStepType::K)
                    os << " = Q_PA" << stepdir;
                else
                    os << " = Q_PB" << stepdir;

                os << " * " << srcname[0];
            }

            if(it.src[1])
                os << " + " << aoppq << " * " << srcname[1]; 
            if(it.src[2] && it.src[3])
                os << " + " << vrr_const0 << " * (" << srcname[2] << " + " << aover << " * " << srcname[3] << ")"; 
            if(it.src[4] && it.src[5])
                os << " + " << vrr_const1 << " * (" << srcname[4] << " + " << aover << " * " << srcname[5] << ")"; 
            if(it.src[6])
                os << " + " << vrr_const2 << " * " << srcname[6];
            if(it.src[7])
                os << " + " << vrr_const3 << " * " << srcname[7];
            os << ";\n";
        }
        
        os << "\n";

    }

    os << indent5 << "}\n";
}





///////////////////////////////////////
// Inline VRR Writer
///////////////////////////////////////
ConstantMap OSTEI_VRR_Writer_Inline::GetConstants(void) const
{
    ConstantMap cm;
    for(int i = 1; i <= vrr_algo_.GetMaxInt(); i++)
        cm.emplace(StringBuilder("const_", i), StringBuilder(i));
    return cm;
}


void OSTEI_VRR_Writer_Inline::WriteVRR(std::ostream & os) const
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
            os << indent5 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2p = one_over_2p;\n"; 
        else
            os << indent5 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2p = " << vinfo_.IntConstant(it) << " * one_over_2p;\n"; 
    }

    for(const auto & it : vrr_algo_.GetAllInt_2q())
    {
        if(it == 1)
            os << indent5 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2q = one_over_2q;\n"; 
        else
            os << indent5 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2q = " << vinfo_.IntConstant(it) << " * one_over_2q;\n"; 
    }

    for(const auto & it : vrr_algo_.GetAllInt_2pq())
    {
        if(it == 1)
            os << indent5 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2pq = one_over_2pq;\n"; 
        else
            os << indent5 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2pq = " << vinfo_.IntConstant(it) << " * one_over_2pq;\n"; 
    }


    // iterate over increasing am
    for(const auto & qam : vrr_algo_.GetAMOrder())
    {
        // Write out the steps
        if(qam != QAM{0,0,0,0})
            WriteVRRSteps_(os, qam, vrr_algo_.GetSteps(qam), 
                           std::to_string(vrr_algo_.GetMReq(qam)+1));

        os << "\n\n";
    }

    os << "\n\n";
}


void OSTEI_VRR_Writer_Inline::WriteVRRFile(std::ostream & os, std::ostream & osh) const
{ }



///////////////////////////////////////
// External VRR Writer
///////////////////////////////////////
void OSTEI_VRR_Writer_External::WriteVRR(std::ostream & os) const
{
    os << "\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << indent5 << "// Primitive integrals: Vertical recurrance\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << "\n";

    // iterate over increasing am
    for(const auto & am : vrr_algo_.GetAMOrder())
    {
        if(am != QAM{0,0,0,0})
        {
            RRStepType rrstep = vrr_algo_.GetRRStep(am);

            os << indent5 << "VRR_" << rrstep_char.at(rrstep) << "_" 
                          << amchar[am[0]] << "_" << amchar[am[1]] << "_"
                          << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

            os << indent7 << PrimVarName(am) << ",\n";

            for(const auto & it : vrr_algo_.GetAMReq(am))
                os << indent7 << PrimVarName(it) << ",\n";
            
            for(const auto & it : vrr_algo_.GetVarReq(am))
                os << indent7 << it << ",\n";

            os << indent7 << (vrr_algo_.GetMReq(am)+1) << ");\n";
            os << "\n";
        }
    }

    os << "\n\n";
}


void OSTEI_VRR_Writer_External::WriteVRRFile(std::ostream & os, std::ostream & osh) const
{
    QAM am = info_.FinalAM();

    os << "//////////////////////////////////////////////\n";
    os << "// VRR: ( " << amchar[am[0]] << " " << amchar[am[1]] << " | " << amchar[am[2]] << " " << amchar[am[3]] << " )\n";
    os << "//////////////////////////////////////////////\n";

    RRStepType rrstep = vrr_algo_.GetRRStep(am);

    std::stringstream prototype;

    char stepchar = rrstep_char.at(rrstep);

    prototype << "void VRR_" << stepchar << "_"
              << amchar[am[0]] << "_" << amchar[am[1]] << "_"
              << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

    // final target
    prototype << indent3 << vinfo_.DoubleType() << " * const restrict " << PrimVarName(am) << ",\n";

    for(const auto & it : vrr_algo_.GetAMReq(am))
        prototype << indent3 << vinfo_.ConstDoubleType() << " * const restrict " << PrimVarName(it) << ",\n";
    
    for(const auto & it : vrr_algo_.GetVarReq(am))
    {
        if(IsPointer(it))
            prototype << indent3 << vinfo_.ConstDoubleType() << " * " << it << ",\n";
        else
            prototype << indent3 << vinfo_.ConstDoubleType() << " " << it << ",\n";
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
            os << indent1 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2p = " 
                          << vinfo_.DoubleSet1(std::to_string(it)) << " * one_over_2p;\n"; 

        for(const auto & it : vrr_algo_.GetIntReq_2q(am))
            os << indent1 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2q = "
                          << vinfo_.DoubleSet1(std::to_string(it)) << " * one_over_2q;\n"; 

        for(const auto & it : vrr_algo_.GetIntReq_2pq(am))
            os << indent1 << vinfo_.ConstDoubleType() << " vrr_const_" << it << "_over_2pq = "
               << vinfo_.DoubleSet1(std::to_string(it)) << " * one_over_2pq;\n"; 

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
