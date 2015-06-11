#include <fstream>
#include <stdexcept>

#include "generator/Boys.hpp"
#include "generator/WriterBase.hpp"
#include "generator/Helpers.hpp"
        
////////////////////////////////////
// Base class
////////////////////////////////////
std::vector<std::string> BoysGen::Includes(void) const
{
    return std::vector<std::string>();
}


void BoysGen::WriteIncludes(std::ostream & os) const
{
    auto boysinc = Includes();
    for(const auto & it : boysinc)
        os << "#include \"" << it << "\"\n";
}


////////////////////////////////////
// FO fitting
////////////////////////////////////
BoysFO::BoysFit::BoysFit(const std::string & filepath)
{
    std::ifstream f(filepath.c_str());
    if(!f.is_open())
        throw std::runtime_error(std::string("Error - cannot open ") + filepath);

    f.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);

    f >> v >> n >> m;

    a.resize(n+1);
    b.resize(m+2);

    std::string s;

    // n is the order, so there is n+1 terms with
    // the leading a0 term
    for(int i = 0; i < n+1; i++)
        f >> a[i];

    // note that the leading 1 term is not included
    // in the file, but the trailing bp is
    b[0] = "1.0";
    for(int i = 1; i < m+2; i++)
        f >> b[i];

    f.close();
}

void BoysFO::WriteBoysSingle_(std::ostream & os, const WriterBase & base, int m, bool prefac) const
{
    const BoysFit & bf = bfmap_.at(m); 

    std::stringstream ssvar, sspow;
    ssvar << base.PrimVarName({0,0,0,0}) << "[" << bf.v << "]";
    sspow << bf.v << ".0+0.5";

    os << indent6 << ssvar.str() << " =\n";
    os << indent6 << "              (\n";
    os << indent6 << "                 (\n";
    os << indent6 << "                   (\n";
    os << indent6 << "                               " << base.DoubleSet(bf.a[0]) << "\n";
    for(int i = 1; i < bf.a.size(); i++)
        os << indent6 << "                     + F_x * ( " << base.DoubleSet(bf.a[i]) << "\n";
    os << indent6 << "                             " << std::string(bf.a.size()-1, ')') << "\n";  // prints a bunch of close paren
    os << indent6 << "                   )\n";
    os << indent6 << "                   /\n";
    os << indent6 << "                   (\n";
    os << indent6 << "                               " << base.DoubleSet(bf.b[0]) << "\n";
    for(int i = 1; i < bf.b.size(); i++)
        os << indent6 << "                     + F_x * ( " << base.DoubleSet(bf.b[i]) << "\n";
    os << indent6 << "                             " << std::string(bf.b.size()-1, ')') << "\n";  // prints a bunch of close paren
    os << indent6 << "                   )\n";
    os << indent6 << "                 )\n";
    os << indent6 << "              );\n";
    os << "\n";

    // apply prefac and power
    os << indent6 << ssvar.str() << " = ";
    if(prefac)
        os << "allprefac * ";
    os << base.Power(ssvar.str(), base.DoubleSet(sspow.str())) << ";\n";
}

void BoysFO::WriteBoys(std::ostream & os, const WriterBase & base) const
{
    std::string primname = base.PrimVarName({0,0,0,0});

    // just calculate them all if L is small
    // value of 3 found by testing
    if(base.L() < 3)
    {
        for(int m = 0; m <= base.L(); m++)
            WriteBoysSingle_(os, base, m, true);
    }
    else
    {
        // calculate highest, and recurse down
        WriteBoysSingle_(os, base, base.L(), false);

        // calculate the downward recursion factors
        os << indent6 << base.ConstDoubleType() << " x2 = " << base.DoubleSet("2.0") << " * F_x;\n";
        os << indent6 << base.ConstDoubleType() << " ex = " << base.Exp("-F_x") << ";\n";

        for(int m = base.L()-1; m >= 0; m--)
        {
            std::stringstream ss;
            ss.precision(18);
            ss << (1.0/(2.0*m+1.0));
            os << indent6 << primname << "[" << m << "] = (x2 * " << primname << "[" << (m+1) << "] + ex) * " << ss.str() << ";\n";
        }

        // add prefac now
        os << "\n";
        os << indent6 << "for(n = 0; n <= " << base.L() << "; ++n)\n";
        os << indent6 << "    " << primname << "[n] *= allprefac;\n";
        os << "\n";
    }
}


BoysFO::BoysFO(std::string dir)
{
    if(dir.back() == '/')
        dir = dir.substr(0, dir.size()-1);

    std::string flist = dir + "/filelist";

    std::ifstream f(flist.c_str());
    if(!f.is_open())
        throw std::runtime_error(std::string("Error - cannot open ") + flist);

    f.exceptions(std::ifstream::badbit);

    std::string fline;

    while(std::getline(f, fline).good())
    {
        BoysFit bf(dir + "/" + fline);
        bfmap_[bf.v] = bf;
    }
}


////////////////////////////////////
// Boys Split
////////////////////////////////////

void BoysSplit::WriteBoys(std::ostream & os, const WriterBase & base) const
{
    os << indent6 << "Boys_F_split(" << base.PrimVarName({0,0,0,0}) << ", " << base.L() << ", F_x);\n";
    for(int i = 0; i <= base.L(); i++)
        os << indent6 << base.PrimVarName({0,0,0,0}) << "[" << i << "] *= allprefac;\n";
}

std::vector<std::string> BoysSplit::Includes(void) const
{
    std::vector<std::string> v{"boys/boys_split.h"};
    return v;
}


////////////////////////////////////
// Boys Valeev Reference
////////////////////////////////////

void BoysVRef::WriteBoys(std::ostream & os, const WriterBase & base) const
{
    os << indent6 << "Boys_F_VRef(" << base.PrimVarName({0,0,0,0}) << ", " << base.L() << ", F_x);\n";
    for(int i = 0; i <= base.L(); i++)
        os << indent6 << base.PrimVarName({0,0,0,0}) << "[" << i << "] *= allprefac;\n";
}

std::vector<std::string> BoysVRef::Includes(void) const
{
    std::vector<std::string> v{"boys/boys_vref.h"};
    return v;
}

