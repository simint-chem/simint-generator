#include "generator/ETWriter.hpp"
#include "generator/WriterBase.hpp"


ETWriter::ETWriter(const ETStepList & etsl) 
          : etsl_(etsl)
{
    for(const auto & it : etsl)
    {
        // should only add the first element of the bra
        if(it.src1 && (it.src1.flags & QUARTET_ETTOPLEVEL))
            etrm_[it.src1.bra.left.am()].insert(it.src1.bra.left);
        if(it.src2 && (it.src2.flags & QUARTET_ETTOPLEVEL))
            etrm_[it.src2.bra.left.am()].insert(it.src2.bra.left);
        if(it.src3 && (it.src3.flags & QUARTET_ETTOPLEVEL))
            etrm_[it.src3.bra.left.am()].insert(it.src3.bra.left);
        if(it.src4 && (it.src4.flags & QUARTET_ETTOPLEVEL))
            etrm_[it.src4.bra.left.am()].insert(it.src4.bra.left);
            
        // add all integrals for this step to the set
        if(it.src1)
            etint_.insert(it.src1.amlist());
        if(it.src2)
            etint_.insert(it.src2.amlist());
        if(it.src3)
            etint_.insert(it.src3.amlist());
        if(it.src3)
            etint_.insert(it.src4.amlist());
        if(it.target)
            etint_.insert(it.target.amlist());
    }       
}



QAMSet ETWriter::ETInt(void) const
{
    return etint_;
}



GaussianMap ETWriter::ETRMap(void) const
{
    return etrm_;
}



bool ETWriter::HasET(void) const
{
    return (etsl_.size() > 0);
}



void ETWriter::DeclarePointers(std::ostream & os, const WriterBase & base) const
{
    if(etint_.size() > 0)
    {
        os << "            // set up pointers to the contracted integrals - Electron Transfer\n";

        for(const auto & it : etint_)
        {
            if(base.IsContArray(it))
                os << "            double * const restrict " << base.PrimPtrName(it) << " = "
                   << base.ArrVarName(it) << " + (abcd * " << (NCART(it[0])*NCART(it[1])*NCART(it[2])*NCART(it[3])) << ");\n";
        }
    }
}



void ETWriter::DeclarePrimArrays(std::ostream & os, const WriterBase & base) const
{
    if(etint_.size())
    {
        os << "                    // Holds temporary integrals for electron transfer\n";

        for(const auto & it : etint_)
        {
            // only if these aren't from vrr
            if(it[1] > 0 || it[2] > 0 || it[3] > 0)
                os << "                    double " << base.PrimVarName(it) << "[" << NCART(it[0]) * NCART(it[2]) << "];\n";
        } 

        os << "\n\n";

    }
}



void ETWriter::WriteETInline(std::ostream & os, const WriterBase & base) const
{
    os << "\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "                    // Primitive integrals: Electron transfer\n";
    os << "                    //////////////////////////////////////////////\n";
    os << "\n";

    std::string indent1(20, ' ');
    std::string indent2(24, ' ');

    if(etsl_.size() == 0)
        os << indent1 << "//...nothing to do...\n";
    else
    {
        for(const auto & it : etsl_)
        {
            os << indent1 << "// " << it << "\n";
            os << ETStepString(it, base);
            os << "\n";
        }
    }

    // add to needed contracted integrals
    for(const auto & it : etint_)
    {
        if(base.IsContArray(it))
        {
            int ncart = NCART(it[0])*NCART(it[2]);

            os << "\n";
            os << indent1 << "// Accumulating in contracted workspace\n";
            os << indent1 << "for(n = 0; n < " << ncart << "; n++)\n";

            if(base.IsFinalAM(it))
                os << indent2 << base.ArrVarName(it) << "[n] += " << base.PrimVarName(it) << "[n];\n";
            else
                os << indent2 << base.ArrVarName(it) << "[n * nshell1234 + abcd] += " << base.PrimVarName(it) << "[n];\n";
        }
    }

}



std::string ETWriter::ETStepVar(const Quartet & q, const WriterBase & base)
{
    std::stringstream ss; 
    ss << base.PrimVarName(q.amlist()) << "[" << q.idx() << "]";
    return ss.str();
}



std::string ETWriter::ETStepString(const ETStep & et, const WriterBase & base)
{
    int stepidx = XYZStepToIdx(et.xyz);
    int ival = et.target.bra.left.ijk[stepidx];
    int kval = et.target.ket.left.ijk[stepidx]-1;

    std::stringstream ss;
    ss << std::string(20, ' ');
    ss << ETStepVar(et.target, base);

    ss << " = ";

    ss << "etfac[" << stepidx << "] * " << ETStepVar(et.src1, base);

    if(et.src2.bra.left && et.src2.ket.left)
        ss << " + " << ival << " * one_over_2q * " << ETStepVar(et.src2, base);
    if(et.src3.bra.left && et.src3.ket.left)
        ss << " + " << kval << " * one_over_2q * " << ETStepVar(et.src3, base);
    if(et.src4.bra.left && et.src4.ket.left)
        ss << " - p_over_q * " << ETStepVar(et.src4, base);
    ss << ";\n";

    return ss.str();
}

