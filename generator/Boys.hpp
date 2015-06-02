#ifndef BOYS_HPP
#define BOYS_HPP

#include <string>
#include <vector>
#include <map>

class WriterBase;

class BoysGen
{
    public:
        void WriteIncludes(std::ostream & os) const;

        virtual void WriteBoys(std::ostream & os, const WriterBase & base) const = 0;

        virtual std::vector<std::string> includes(void) const;
        virtual ~BoysGen() { };
};

class BoysFO : public BoysGen
{
    public:
        BoysFO(std::string dir); // read from directory

        virtual void WriteBoys(std::ostream & os, const WriterBase & base) const;

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
        };

        std::map<int, BoysFit> bfmap_;
};


struct BoysSplit : public BoysGen
{
    public:
        // default constructors ok

        virtual void WriteBoys(std::ostream & os, const WriterBase & base) const;
        virtual std::vector<std::string> includes(void) const;
};


struct BoysVRef : public BoysGen
{
    public:
        // default constructors ok

        virtual void WriteBoys(std::ostream & os, const WriterBase & base) const;
        virtual std::vector<std::string> includes(void) const;
};


#endif
