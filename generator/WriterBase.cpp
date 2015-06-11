#include <stdexcept>
#include <fstream>
#include <cctype> // for tolower
#include "generator/WriterBase.hpp"


WriterBase::WriterBase(const OptionsMap & options, const std::string & prefix, const QAM & finalam)
     : options_(options), prefix_(prefix), finalam_(finalam)
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

    // add memory for AB_{xyz} and CD_{xyz}
    if(HasBraHRR())
        memory_ += 3 * sizeof(double);
    if(HasKetHRR())
        memory_ += 3 * sizeof(double);
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
    {
        std::stringstream ss;
        ss << "No value for option: " << option;
        throw std::runtime_error(ss.str());
    }
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
        

/*
void WriterBase::PermuteResult(std::ostream & os, const std::string & src) const
{
    int ncart = NCART(finalam_[0]) * NCART(finalam_[1]) * NCART(finalam_[2]) * NCART(finalam_[3]);
    os << "\n";
    os << "        //Permute the result\n";
    os << "        for(ir = 0; ir < " << ncart << "; ir++)\n";
    os << "            result[abcd * " << ncart << " + ir] = " << src << "[ir];\n"; 
    os << "\n";
}
*/

void WriterBase::DeclareContwork(std::ostream & os) const
{

    if(memory_ > 0)
    {
        os << "    // Workspace for contracted integrals\n";

        if(memory_ > GetOption(OPTION_STACKMEM))
            os << "    double * const contwork = malloc(SIMINT_NSHELL_SIMD * " << memory_ << ");\n";
        else
            os << "    double contwork[SIMINT_NSHELL_SIMD * " << memory_ << "];\n";
        os << "\n";

        os << "    // partition workspace into shells\n";

        size_t ptidx = 0;
        for(const auto & it : contq_)
        {
            if(it != finalam_)
            {
                os << "    double * const " << ArrVarName(it) << " = contwork + (SIMINT_NSHELL_SIMD * " << ptidx << ");\n";
                ptidx += NCART(it[0]) * NCART(it[1]) * NCART(it[2]) * NCART(it[3]);
            }
        
        }
        os << "\n";

        // AB_ and CD_ arrays
        if(HasBraHRR())
        {
            os << "    // Holds AB_{xyz} in a flattened fashion for later\n";
            os << "    double * const restrict AB_x = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict AB_y = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict AB_z = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "\n";
        }

        if(HasKetHRR())
        {
            os << "    // Holds CD_{xyz} in a flattened fashion for later\n";
            os << "    double * const restrict CD_x = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict CD_y = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "    double * const restrict CD_z = contwork + (SIMINT_NSHELL_SIMD * " << ptidx++ << ");\n";
            os << "\n";
        }
    }
}


void WriterBase::ZeroContWork(std::ostream & os, const std::string & nshell) const
{
    if(memory_ > 0)
        os << "            memset(contwork, 0, " << nshell << " * " << memory_ << ");\n";
}


void WriterBase::FreeContwork(std::ostream & os) const
{
    if((memory_ > 0) && (memory_ > GetOption(OPTION_STACKMEM)))
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

const std::string & WriterBase::Prefix(void) const
{
    return prefix_;
}


//bool WriterBase::Permute(void) const
//{
//    return (GetOption(OPTION_PERMUTE) > 0 && HasHRR());
//}


void WriterBase::ReadCPUFlags(const std::string & file)
{
    std::set<std::string> cpuflags_tmp;

    std::ifstream cpufile(file.c_str());
    if(!cpufile.is_open())
        throw std::runtime_error("Error - cannot open cpu flags file!");

    while(cpufile)
    {
        std::string s;
        cpufile >> s;

        // convert to lower case
        for(auto & it : s)
            it = tolower(it);

        if(s.size())
            cpuflags_tmp.insert(s);
    }    

    /*
    std::cout << "Read " << cpuflags_tmp.size() << " flags\n";
    for(const auto & it : cpuflags_tmp)
        std::cout << "  \"" << it << "\"\n";
    std::cout << "\n";
    */

    // add equivalent ones where underscore is replaced by a period
    cpuflags_ = cpuflags_tmp;

    for(auto it : cpuflags_tmp)
    {
        for(auto & it2 : it)
        {
            if(it2 == '_')
                it2 = '.';
        }

        // its a set, so inserting duplicates is ok
        cpuflags_.insert(it);
    }

    /*
    std::cout << "Processed into " << cpuflags_.size() << " flags\n";
    for(const auto & it : cpuflags_)
        std::cout << "  \"" << it << "\"\n";
    std::cout << "\n";
    */
}
        
bool WriterBase::HasCPUFlag(const std::string & flag) const
{
    return (cpuflags_.count(flag) > 0);
}
