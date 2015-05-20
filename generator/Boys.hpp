#ifndef BOYS_HPP
#define BOYS_HPP

#include <memory>
#include <string>
#include <vector>
#include <map>

class BoysGen
{
    public:
        virtual std::string all_code_lines(int maxam) const = 0;

        virtual std::vector<std::string> includes(void) const;
        virtual ~BoysGen() { };
};

class BoysFO : public BoysGen
{
    public:
        BoysFO(std::string dir); // read from directory

        virtual std::string all_code_lines(int maxam) const;

    private:
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

        std::map<int, BoysFit> bfmap_;
};


struct BoysSplit : public BoysGen
{
    public:
        // default constructors ok

        virtual std::string all_code_lines(int maxam) const;
        virtual std::vector<std::string> includes(void) const;
};




#endif
