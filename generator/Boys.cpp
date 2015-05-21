#include <fstream>
#include <stdexcept>
#include <sstream>

#include "generator/Boys.hpp"
        
////////////////////////////////////
// Base class
////////////////////////////////////
std::vector<std::string> BoysGen::includes(void) const
{
    return std::vector<std::string>();
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

std::string BoysFO::all_code_lines(int maxam) const
{
    std::stringstream ss;
    for(int i = 0; i <= maxam; i++)
        ss << bfmap_.at(i).code_line() << "\n";
    return ss.str();
}

std::string BoysFO::BoysFit::code_line(void) const
{
    const std::string indent(20, ' ');
    std::stringstream ss;
    ss << indent << "AUX_INT__s_s_s_s[" << v << "] = allprefac\n"; 
    ss << indent << "         * pow(\n";
    ss << indent << "                 (\n";
    ss << indent << "                   (\n";
    ss << indent << "                               " << a[0] << "\n";
    for(int i = 1; i < a.size(); i++)
        ss << indent << "                     + F_x * ( " << a[i] << "\n";
    ss << indent << "                             " << std::string(a.size()-1, ')') << "\n";
    ss << indent << "                   )\n";
    ss << indent << "                   /\n";
    ss << indent << "                   (\n";
    ss << indent << "                               " << b[0] << "\n";
    for(int i = 1; i < b.size(); i++)
        ss << indent << "                     + F_x * ( " << b[i] << "\n";
    ss << indent << "                             " << std::string(b.size()-1, ')') << "\n";
    ss << indent << "                   )\n";
    ss << indent << "                 ), " << v << ".0+0.5);\n";
    return ss.str();  
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

std::string BoysSplit::all_code_lines(int maxam) const
{
    const std::string indent(20, ' ');
    std::stringstream ss;
    ss << indent << "Boys_F_split(AUX_INT__s_s_s_s, " << maxam << ", F_x);\n";
    for(int i = 0; i <= maxam; i++)
        ss << indent << "AUX_INT__s_s_s_s[" << i << "] *= allprefac;\n";

    return ss.str();
}

std::vector<std::string> BoysSplit::includes(void) const
{
    std::vector<std::string> v{"boys/boys_split.h"};
    return v;
}


////////////////////////////////////
// Boys Valeev Reference
////////////////////////////////////

std::string BoysVRef::all_code_lines(int maxam) const
{
    const std::string indent(20, ' ');
    std::stringstream ss;
    ss << indent << "Boys_F_VRef(AUX_INT__s_s_s_s, " << maxam << ", F_x);\n";
    for(int i = 0; i <= maxam; i++)
        ss << indent << "AUX_INT__s_s_s_s[" << i << "] *= allprefac;\n";

    return ss.str();
}

std::vector<std::string> BoysVRef::includes(void) const
{
    std::vector<std::string> v{"boys/boys_vref.h"};
    return v;
}

