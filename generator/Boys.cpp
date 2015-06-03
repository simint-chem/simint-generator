#include <fstream>
#include <stdexcept>

#include "generator/Boys.hpp"
#include "generator/WriterBase.hpp"
        
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

void BoysFO::WriteBoys(std::ostream & os, const WriterBase & base) const
{
    const std::string indent(20, ' ');

    for(int i = 0; i <= base.L(); i++)
    {
        const BoysFit & bf = bfmap_.at(i); 

        os << indent << base.AuxName(0) << "[" << bf.v << "] = allprefac\n"; 
        os << indent << "         * pow(\n";
        os << indent << "                 (\n";
        os << indent << "                   (\n";
        os << indent << "                               " << bf.a[0] << "\n";
        for(int i = 1; i < bf.a.size(); i++)
            os << indent << "                     + F_x * ( " << bf.a[i] << "\n";
        os << indent << "                             " << std::string(bf.a.size()-1, ')') << "\n";
        os << indent << "                   )\n";
        os << indent << "                   /\n";
        os << indent << "                   (\n";
        os << indent << "                               " << bf.b[0] << "\n";
        for(int i = 1; i < bf.b.size(); i++)
            os << indent << "                     + F_x * ( " << bf.b[i] << "\n";
        os << indent << "                             " << std::string(bf.b.size()-1, ')') << "\n";
        os << indent << "                   )\n";
        os << indent << "                 ), " << bf.v << ".0+0.5);\n";
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
    const std::string indent(20, ' ');
    os << indent << "Boys_F_split(" << base.AuxName(0) << ", " << base.L() << ", F_x);\n";
    for(int i = 0; i <= base.L(); i++)
        os << indent << base.AuxName(0) << "[" << i << "] *= allprefac;\n";
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
    const std::string indent(20, ' ');
    os << indent << "Boys_F_VRef(" << base.AuxName(0) << ", " << base.L() << ", F_x);\n";
    for(int i = 0; i <= base.L(); i++)
        os << indent << base.AuxName(0) << "[" << i << "] *= allprefac;\n";
}

std::vector<std::string> BoysVRef::Includes(void) const
{
    std::vector<std::string> v{"boys/boys_vref.h"};
    return v;
}

