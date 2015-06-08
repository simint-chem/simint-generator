#include "generator/WriterBase.hpp"


WriterBase::WriterBase(const OptionsMap & options, const QAM & finalam)
     : options_(options), finalam_(finalam)
{ }



void WriterBase::SetContQ(const QAMSet & topquartets)
{
    contq_.clear();

    // add hrr top level stuff
    for(const auto & it : topquartets)
        contq_.insert(it);

    // calculate the memory
    memory_ = 0;
    for(const auto & it : contq_)
    {
        if(it != finalam_)
            memory_ += (sizeof(double) * NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]));
    }
}



size_t WriterBase::MemoryReq(void) const
{
    return memory_;
}



int WriterBase::L(void) const
{
    return finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3];
}



int WriterBase::GetOption(int option) const
{
    if(options_.count(option) > 0)
        return options_.at(option);
    else
        return 0;
}



bool WriterBase::IsContArray(const QAM & am) const
{
    return contq_.count(am);
}



bool WriterBase::IsFinalAM(const QAM & am) const
{
    return am == finalam_;
}


QAM WriterBase::FinalAM(void) const
{
    return finalam_;
}



std::string WriterBase::ArrVarName(const QAM & am, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]];
    return ss.str();
}

std::string WriterBase::ArrVarName(int am1, int am2, const std::string & ketstr, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << amchar[am1] << "_" << amchar[am2] << "_" << ketstr;
    return ss.str();
}

std::string WriterBase::ArrVarName(const std::string & brastr, int am3, int am4, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << brastr << "_" << amchar[am3] << "_" << amchar[am4];
    return ss.str();
}


std::string WriterBase::HRRVarName(const QAM & am)
{
    return ArrVarName(am, "HRR");
}


std::string WriterBase::ArrVarName(int am1, int am2, const std::string & ketstr)
{
    return ArrVarName(am1, am2, ketstr, "HRR");
}

std::string WriterBase::ArrVarName(const std::string & brastr, int am3, int am4)
{
    return ArrVarName(brastr, am3, am4, "HRR");
}


std::string WriterBase::PrimVarName(const QAM & am)
{
    return ArrVarName(am, "PRIM");
}


std::string WriterBase::PrimPtrName(const QAM & am)
{
    return ArrVarName(am, "PRIM_PTR");
}
        

void WriterBase::PermuteResult(std::ostream & os, const std::string & src) const
{
    int ncart = NCART(finalam_[0]) * NCART(finalam_[1]) * NCART(finalam_[2]) * NCART(finalam_[3]);
    os << "\n";
    os << "        //Permute the result\n";
    os << "        for(ir = 0; ir < " << ncart << "; ir++)\n";
    os << "            result[abcd * " << ncart << " + ir] = " << src << "[ir];\n"; 
    os << "\n";
}


void WriterBase::DeclareContwork(std::ostream & os) const
{
    if(memory_ > 0)
    {
        os << "    // Workspace for contracted integrals\n";
        os << "    double * const contwork = malloc(nshell1234 * " << memory_ << ");\n";
        os << "    memset(contwork, 0, nshell1234 * " << memory_ << ");\n";
        os << "\n";
        os << "    // partition workspace into shells\n";

        size_t ptidx = 0;
        for(const auto & it : contq_)
        {
            if(it != finalam_)
            {
                os << "    double * const " << ArrVarName(it) << " = contwork + (nshell1234 * " << ptidx << ");\n";
                ptidx += NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]);
            }
        
        }
    }
}

void WriterBase::FreeContwork(std::ostream & os) const
{
    if(memory_ > 0)
    {
        os << "    // Free contracted work space\n";
        os << "    free(contwork);\n";
        os << "\n";
    }
}

bool WriterBase::HasVRR(void) const
{
    return ((finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3]) > 0);
}

bool WriterBase::HasET(void) const
{
    return (finalam_[2]+finalam_[3] > 0);
}

bool WriterBase::HasHRR(void) const
{
    return (HasBraHRR() || HasKetHRR());
}

bool WriterBase::HasBraHRR(void) const
{
    return (finalam_[1] > 0);
}

bool WriterBase::HasKetHRR(void) const
{
    return (finalam_[3] > 0);
}

bool WriterBase::Permute(void) const
{
    return (GetOption(OPTION_PERMUTE) > 0 && HasHRR());
}
