#include <fstream>
#include <stdexcept>

#include "generator/Boys.hpp"
#include "generator/WriterInfo.hpp"
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



void BoysGen::AddConstants(void) const
{
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

std::string BoysFO::GetFOConstant(std::string ab, int m, int i) const
{
    std::stringstream cname;
    cname << "FO_" << m << "_" << ab << "_" << i;
    return WriterInfo::NamedConstant(cname.str());
}



void BoysFO::WriteBoysSingle_(std::ostream & os, int m, bool prefac) const
{
    const BoysFit & bf = bfmap_.at(m); 

    std::stringstream ssvar, sspow;
    ssvar << WriterInfo::PrimVarName({0,0,0,0}) << "[" << bf.v << "]";
    sspow << bf.v << ".0+0.5";


    if(WriterInfo::HasFMA())
    {
        // Horner's rule
        // Numerator
        os << indent5 << "num = " << WriterInfo::FMAdd("F_x", GetFOConstant("a", m, bf.a.size()-1), GetFOConstant("a", m, bf.a.size()-2)) << ";\n";
        for(int i = bf.a.size()-3; i >= 0; i--)
            os << indent5 << "num = " << WriterInfo::FMAdd("F_x", "num", GetFOConstant("a", m, i)) << ";\n";
        os << "\n";

        // Denominator
        os << indent5 << "den = " << WriterInfo::FMAdd("F_x", GetFOConstant("b", m, bf.b.size()-1), GetFOConstant("b", m, bf.b.size()-2)) << ";\n"; 
        for(int i = bf.b.size()-3; i >= 1; i--)
            os << indent5 << "den = " << WriterInfo::FMAdd("F_x", "den", GetFOConstant("b", m, i)) << ";\n";
        os << indent5 << "den = " << WriterInfo::FMAdd("F_x", "den", WriterInfo::IntConstant(1)) << ";\n";
        os << "\n";
        os << indent5 << ssvar.str() << " = num / den;\n";
    }
    else
    {
        os << indent5 << ssvar.str() << " =\n";
        os << indent5 << "              (\n";
        os << indent5 << "                 (\n";
        os << indent5 << "                   (\n";
        os << indent5 << "                               " << GetFOConstant("a", m, 0) << "\n";
        for(size_t i = 1; i < bf.a.size(); i++)
            os << indent5 << "                     + F_x * ( " << GetFOConstant("a", m, i) << "\n";

        os << indent5 << "                             " << std::string(bf.a.size()-1, ')') << "\n";  // prints a bunch of close paren
        os << indent5 << "                   )\n";
        os << indent5 << "                   /\n";
        os << indent5 << "                   (\n";
        os << indent5 << "                               " << WriterInfo::IntConstant(1) << "\n";   //<< GetFOConstant("b", 0, m) << "\n";
        for(size_t i = 1; i < bf.b.size(); i++)
            os << indent5 << "                     + F_x * ( " << GetFOConstant("b", m, i) << "\n";

        os << indent5 << "                             " << std::string(bf.b.size()-1, ')') << "\n";  // prints a bunch of close paren
        os << indent5 << "                   )\n";
        os << indent5 << "                 )\n";
        os << indent5 << "              );\n";
        os << "\n";
    }


    // apply prefac and power
    // calculate the prefactor if this is the first time it's needed
    if(prefac && m == 0)
        os << indent5 << WriterInfo::ConstDoubleType() << " prefac = P_prefac * " << WriterInfo::DoubleLoad("Q.prefac", "j") << ";\n";



    os << indent5 << ssvar.str() << " = ";
    
    if(prefac)
        os << "prefac * " << WriterInfo::Sqrt(ssvar.str() + " * one_over_PQalpha_sum");
    else
        os << WriterInfo::Sqrt(ssvar.str());

    for(int i = 0; i < bf.v; i++)
        os << " * " << ssvar.str();
    os << ";\n";
    os << "\n";
}



void BoysFO::AddConstants(void) const
{
    const int L = WriterInfo::L();

    // for b0
    WriterInfo::AddIntConstant(1);

    if(L < 3)
    {
        for(int m = 0; m <= L; m++)
        {
            const BoysFit & bf = bfmap_.at(m);
            for(size_t i = 0; i < bf.a.size(); i++)
            { 
                std::stringstream cname;
                cname << "FO_" << m << "_a_" << i;
                WriterInfo::AddNamedConstant(cname.str(), bf.a[i]);
            }        
            for(size_t i = 1; i < bf.b.size(); i++) // skip b0 = 1.0
            { 
                std::stringstream cname;
                cname << "FO_" << m << "_b_" << i;
                WriterInfo::AddNamedConstant(cname.str(), bf.b[i]);
            }        
        }

    } 
    else
    {
        // add the factor of 2.0
        WriterInfo::AddIntConstant(2);

        const BoysFit & bf = bfmap_.at(L);
        for(size_t i = 0; i < bf.a.size(); i++)
        { 
            std::stringstream cname;
            cname << "FO_" << L << "_a_" << i;
            WriterInfo::AddNamedConstant(cname.str(), bf.a[i]);
        }        
        for(size_t i = 1; i < bf.b.size(); i++) // skip b0 = 1.0
        { 
            std::stringstream cname;
            cname << "FO_" << L << "_b_" << i;
            WriterInfo::AddNamedConstant(cname.str(), bf.b[i]);
        }        

        // constants for the recursion
        // skipping i = 0 since that is 1.0
        for(int i = L-1; i > 0; i--)
        {
            std::stringstream cname, ssval;
            cname << "FO_RECUR_" << i;
            ssval.precision(18);
            ssval << (1.0/(2.0*i+1.0));
            WriterInfo::AddNamedConstant(cname.str(), ssval.str());
        }        
    }
       
}

void BoysFO::WriteBoys(std::ostream & os) const
{
    std::string primname = WriterInfo::PrimVarName({0,0,0,0});

    // if using FMA, we need to declare the numerator and denominator separately
    if(WriterInfo::HasFMA())
        os << indent5 << WriterInfo::DoubleType() << " num, den;\n\n";

    // just calculate them all if L is small
    // value of 3 found by testing
    if(WriterInfo::L() < 3)
    {
        for(int m = 0; m <= WriterInfo::L(); m++)
            WriteBoysSingle_(os, m, true);
    }
    else
    {
        // calculate highest, and recurse down
        WriteBoysSingle_(os, WriterInfo::L(), false);

        // calculate the downward recursion factors
        os << indent5 << WriterInfo::ConstDoubleType() << " x2 = " << WriterInfo::IntConstant(2) << " * F_x;\n";
        os << indent5 << WriterInfo::ConstDoubleType() << " ex = " << WriterInfo::Exp("-F_x") << ";\n";

        for(int m = WriterInfo::L()-1; m > 0; m--)
        {
            std::stringstream cname;
            cname << "FO_RECUR_" << m;
            os << indent5 << primname << "[" << m << "] = (x2 * " << primname << "[" << (m+1) << "] + ex) * " << WriterInfo::NamedConstant(cname.str()) << ";\n";
        }

        // do m = 0
        os << indent5 << primname << "[0] = (x2 * " << primname << "[1] + ex);\n"; // times 1.0


        // add prefac now
        os << indent5 << WriterInfo::ConstDoubleType() << " prefac = " << WriterInfo::Sqrt("one_over_PQalpha_sum") << " * P_prefac * " << WriterInfo::DoubleLoad("Q.prefac", "j") << ";\n";
 
        os << "\n";
        os << indent5 << "for(n = 0; n <= " << WriterInfo::L() << "; ++n)\n";
        os << indent6 << primname << "[n] *= prefac;\n";
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

void BoysSplit::WriteBoys(std::ostream & os) const
{
    os << indent5 << "Boys_F_split_simd((double *)" << WriterInfo::PrimVarName({0,0,0,0}) << ", " 
                  << WriterInfo::L() << ", (const double *)(&F_x));\n";

    os << indent5 << WriterInfo::ConstDoubleType() << " prefac = " << WriterInfo::Sqrt("one_over_PQalpha_sum") << " * P_prefac * " << WriterInfo::DoubleLoad("Q.prefac", "j") << ";\n";

    for(int i = 0; i <= WriterInfo::L(); i++)
        os << indent5 << WriterInfo::PrimVarName({0,0,0,0}) << "[" << i << "] *= prefac;\n";
}

std::vector<std::string> BoysSplit::Includes(void) const
{
    std::vector<std::string> v{"boys/boys_split.h"};
    return v;
}


////////////////////////////////////
// Boys Valeev Reference
////////////////////////////////////

void BoysVRef::WriteBoys(std::ostream & os) const
{
    os << indent5 << "Boys_F_VRef(" << WriterInfo::PrimVarName({0,0,0,0}) << ", " << WriterInfo::L() << ", F_x);\n";
    for(int i = 0; i <= WriterInfo::L(); i++)
        os << indent5 << WriterInfo::PrimVarName({0,0,0,0}) << "[" << i << "] *= allprefac;\n";
}

std::vector<std::string> BoysVRef::Includes(void) const
{
    std::vector<std::string> v{"boys/boys_vref.h"};
    return v;
}

