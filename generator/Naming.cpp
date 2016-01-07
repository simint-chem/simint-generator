#include "generator/Naming.hpp"

std::string ArrVarName(const QAM & am, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]];
    return ss.str();
}

std::string ArrVarName(int am1, int am2, const std::string & ketstr, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << amchar[am1] << "_" << amchar[am2] << "_" << ketstr;
    return ss.str();
}

std::string ArrVarName(const std::string & brastr, int am3, int am4, const std::string & prefix)
{
    std::stringstream ss;

    if(prefix.size())
        ss << prefix << "_";

    ss << "INT__"  << brastr << "_" << amchar[am3] << "_" << amchar[am4];
    return ss.str();
}


std::string HRRVarName(const QAM & am)
{
    return ArrVarName(am, "HRR");
}


std::string HRRVarName(int am1, int am2, const std::string & ketstr)
{
    return ArrVarName(am1, am2, ketstr, "HRR");
}

std::string HRRVarName(const std::string & brastr, int am3, int am4)
{
    return ArrVarName(brastr, am3, am4, "HRR");
}


std::string PrimVarName(const QAM & am)
{
    return ArrVarName(am, "PRIM");
}


std::string PrimPtrName(const QAM & am)
{
    return ArrVarName(am, "PRIM_PTR");
}
