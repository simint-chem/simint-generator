#include <fstream>
#include <stdexcept>
#include <sstream>

#include "generator/Boys.hpp"

BoysFit::BoysFit(const std::string & filepath)
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

std::string BoysFit::code_line(void) const
{
    const std::string indent(20, ' ');
    std::stringstream ss;
    ss << indent << "AUX_S_0_0_0_0[" << v << "] = allprefac\n"; 
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

BoysMap ReadBoysFitInfo(std::string dir)
{
    if(dir.back() == '/')
        dir = dir.substr(0, dir.size()-1);

    std::string flist = dir + "/filelist";

    std::ifstream f(flist.c_str());
    if(!f.is_open())
        throw std::runtime_error(std::string("Error - cannot open ") + flist);

    f.exceptions(std::ifstream::badbit);

    BoysMap bm;
    std::string fline;

    while(std::getline(f, fline).good())
    {
        BoysGenPtr bf(new BoysFit(dir + "/" + fline));
        bm[bf->v] = bf;
    }

    return bm;

}

