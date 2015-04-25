#ifndef BOYS_H
#define BOYS_H

#include <string>
#include <vector>
#include <map>

struct BoysFit
{
    int v; // order of the boys function
    int n; // order of the numerator polynomial
    int m; // order of the denominator polynomial

    std::vector<std::string> a; // numerator coefficients
    std::vector<std::string> b; // denominator coefficients

    BoysFit(const std::string & filepath);

    // let the compiler generate these
    BoysFit(const BoysFit & rhs) = default;
    BoysFit() = default;

    std::string code_line(void) const;
};


typedef std::map<int, BoysFit> BoysMap;


BoysMap ReadBoysInfo(std::string dir);


#endif
