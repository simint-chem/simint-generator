/*! \file
 *
 * \brief Class for generating the Boys function via Modified Fructl-Otto
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#include <fstream>
#include <stdexcept>

#include "generator/BoysGenerator.hpp"
#include "generator/GeneratorInfoBase.hpp"
#include "generator/Naming.hpp"
#include "generator/Printing.hpp"
        

////////////////////////////////////
// FO fitting
////////////////////////////////////
ConstantMap BoysFO::GetConstants(void) const
{
    ConstantMap cm;

    const int L = info_.L();

    // for b0
    cm.emplace("const_1", "1");

    if(L < BOYS_FO_RECUR)
    {
        for(int m = 0; m <= L; m++)
        {
            const BoysFit & bf = bfmap_.at(m);
            for(size_t i = 0; i < bf.a.size(); i++)
                cm.emplace(StringBuilder("FO_", m, "_a_", i), bf.a[i]);
            for(size_t i = 1; i < bf.b.size(); i++) // skip b0 = 1.0
                cm.emplace(StringBuilder("FO_", m, "_b_", i), bf.b[i]);
        }

    } 
    else
    {
        // add the factor of 2.0
        cm.emplace("const_2", "2");

        const BoysFit & bf = bfmap_.at(L);
        for(size_t i = 0; i < bf.a.size(); i++)
            cm.emplace(StringBuilder("FO_", L, "_a_", i), bf.a[i]);
        for(size_t i = 1; i < bf.b.size(); i++) // skip b0 = 1.0
            cm.emplace(StringBuilder("FO_", L, "_b_", i), bf.b[i]);


        // constants for the recursion
        // skipping i = 0 since that is 1.0
        for(int i = L-1; i > 0; i--)
            cm.emplace(StringBuilder("FO_RECUR_", i),
                       StringBuilder((1.0/(2.0*i+1.0))));
    }
       
    return cm;

}

IncludeSet BoysFO::GetIncludes(void) const
{
    return {"\"simint/boys/boys_FO.h\""};
}



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

std::string BoysFO::GetFOConstant_(std::string ab, int m, int i) const
{
    return StringBuilder("FO_", m, "_", ab, "_", i);
}



void BoysFO::WriteBoysSingle_(std::ostream & os, int m, bool prefac) const
{
    const BoysFit & bf = bfmap_.at(m); 

    std::string varstr = StringBuilder(PrimVarName({0,0,0,0}), "[", bf.v, "]");
    std::string powstr = StringBuilder(bf.v, ".0+0.5");


    if(info_.HasFMA())
    {
        // Horner's rule
        // Numerator
        os << indent5 << "num = " << vinfo_.FMAdd("F_x", GetFOConstant_("a", m, bf.a.size()-1), GetFOConstant_("a", m, bf.a.size()-2)) << ";\n";
        for(int i = bf.a.size()-3; i >= 0; i--)
            os << indent5 << "num = " << vinfo_.FMAdd("F_x", "num", GetFOConstant_("a", m, i)) << ";\n";
        os << "\n";

        // Denominator
        os << indent5 << "den = " << vinfo_.FMAdd("F_x", GetFOConstant_("b", m, bf.b.size()-1), GetFOConstant_("b", m, bf.b.size()-2)) << ";\n"; 
        for(int i = bf.b.size()-3; i >= 1; i--)
            os << indent5 << "den = " << vinfo_.FMAdd("F_x", "den", GetFOConstant_("b", m, i)) << ";\n";
        os << indent5 << "den = " << vinfo_.FMAdd("F_x", "den", vinfo_.IntConstant(1)) << ";\n";
        os << "\n";
        os << indent5 << varstr << " = num / den;\n";
    }
    else
    {
        os << indent5 << varstr << " =\n";
        os << indent5 << "              (\n";
        os << indent5 << "                 (\n";
        os << indent5 << "                   (\n";
        os << indent5 << "                               " << GetFOConstant_("a", m, 0) << "\n";
        for(size_t i = 1; i < bf.a.size(); i++)
            os << indent5 << "                     + F_x * ( " << GetFOConstant_("a", m, i) << "\n";

        os << indent5 << "                             " << std::string(bf.a.size()-1, ')') << "\n";  // prints a bunch of close paren
        os << indent5 << "                   )\n";
        os << indent5 << "                   /\n";
        os << indent5 << "                   (\n";
        os << indent5 << "                               " << vinfo_.IntConstant(1) << "\n";   //<< GetFOConstant_("b", 0, m) << "\n";
        for(size_t i = 1; i < bf.b.size(); i++)
            os << indent5 << "                     + F_x * ( " << GetFOConstant_("b", m, i) << "\n";

        os << indent5 << "                             " << std::string(bf.b.size()-1, ')') << "\n";  // prints a bunch of close paren
        os << indent5 << "                   )\n";
        os << indent5 << "                 )\n";
        os << indent5 << "              );\n";
        os << "\n";
    }


    // apply prefac and power
    // calculate the prefactor if this is the first time it's needed
    if(prefac && m == 0)
        os << indent5 << vinfo_.ConstDoubleType() << " prefac = P_prefac * Q_prefac;\n";



    os << indent5 << varstr << " = ";
    
    if(prefac)
        os << "prefac * " << vinfo_.Sqrt(varstr + " * one_over_PQalpha_sum");
    else
        os << vinfo_.Sqrt(varstr);

    for(int i = 0; i < bf.v; i++)
        os << " * " << varstr;
    os << ";\n";
    os << "\n";
}



void BoysFO::WriteBoys(std::ostream & os) const
{
    const int L = info_.L();

    std::string primname = PrimVarName({0,0,0,0});

    // if using FMA, we need to declare the numerator and denominator separately
    if(info_.HasFMA())
        os << indent5 << vinfo_.DoubleType() << " num, den;\n\n";

    // just calculate them all if L is small
    // value of 3 found by testing
    if(L < BOYS_FO_RECUR)
    {
        for(int m = 0; m <= L; m++)
            WriteBoysSingle_(os, m, true);
    }
    else
    {
        // calculate highest, and recurse down
        WriteBoysSingle_(os, L, false);

        // calculate the downward recursion factors
        os << indent5 << vinfo_.ConstDoubleType() << " x2 = " << vinfo_.IntConstant(2) << " * F_x;\n";
        os << indent5 << vinfo_.ConstDoubleType() << " ex = " << vinfo_.Exp("-F_x") << ";\n";

        for(int m = L-1; m > 0; m--)
        {
            std::string cname = StringBuilder("FO_RECUR_", m);
            os << indent5 << primname << "[" << m << "] = (x2 * " << primname << "[" << (m+1) << "] + ex) * " << cname << ";\n";
        }

        // do m = 0
        os << indent5 << primname << "[0] = (x2 * " << primname << "[1] + ex);\n"; // times 1.0


        // add prefac now
        os << indent5 << vinfo_.ConstDoubleType() << " prefac = " << vinfo_.Sqrt("one_over_PQalpha_sum") << " * P_prefac * Q_prefac;\n";
 
        os << "\n";
        os << indent5 << "for(n = 0; n <= " << L << "; ++n)\n";
        os << indent6 << primname << "[n] *= prefac;\n";
        os << "\n";
    }
}


BoysFO::BoysFO(const GeneratorInfoBase & info, std::string dir)
    : BoysGenerator(info)
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

