#include "generator/WriterBase.hpp"


WriterBase::WriterBase(const OptionsMap & options, const QAMList & finalam)
     : options_(options), finalam_(finalam)
{ }



void WriterBase::SetContQ(const QuartetSet & topquartets, const DoubletSetMap & topkets) 
{
    // add the final am
    contq_.insert(finalam_);

    // add hrr top level stuff
    for(const auto & it : topquartets)
        contq_.insert(it.amlist());

    // also add final bra combined with all the kets
    for(const auto & it : topkets)
    {
        QAMList qam{finalam_[0], finalam_[1], it.first, 0};
        contq_.insert(qam);
    }

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



bool WriterBase::IsContArray(const QAMList & am) const
{
    return contq_.count(am);
}



QAMList WriterBase::FinalAM(void) const
{
    return finalam_;
}



std::string WriterBase::ArrVarName(const QAMList & am)
{
    std::stringstream ss;
    ss << "INT__"  << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]];
    return ss.str();
}



std::string WriterBase::AuxName(int i)
{
    return std::string("AUX_") + ArrVarName({i, 0, 0, 0});
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

