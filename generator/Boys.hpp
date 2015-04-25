#ifndef BOYS_HPP
#define BOYS_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>

struct BoysGen
{
    int v; // order of the boys function

    virtual ~BoysGen() { };
    virtual std::string code_line(void) const = 0;
};

struct BoysFit : public BoysGen
{
    int n; // order of the numerator polynomial
    int m; // order of the denominator polynomial

    std::vector<std::string> a; // numerator coefficients
    std::vector<std::string> b; // denominator coefficients

    BoysFit(const std::string & filepath);

    // let the compiler generate these
    BoysFit(const BoysFit & rhs) = default;
    BoysFit() = default;

    virtual std::string code_line(void) const;
};


typedef std::shared_ptr<BoysGen> BoysGenPtr;
typedef std::map<int, BoysGenPtr> BoysMap;


BoysMap ReadBoysFitInfo(std::string dir);


#endif
