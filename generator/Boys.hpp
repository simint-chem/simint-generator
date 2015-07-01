#ifndef BOYS_HPP
#define BOYS_HPP

#include <string>
#include <vector>
#include <map>



class BoysGen
{
    public:
        virtual void WriteBoys(std::ostream & os) const = 0;

        void WriteIncludes(std::ostream & os) const;
        void WriteConstants(std::ostream & os) const;

        virtual ~BoysGen() { };

    protected:
        virtual std::vector<std::string> Includes(void) const;
};

class BoysFO : public BoysGen
{
    public:
        BoysFO(std::string dir); // read from directory

        virtual void WriteBoys(std::ostream & os) const;

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

        void WriteBoysSingle_(std::ostream & os, int m, bool prefac) const;
};


struct BoysSplit : public BoysGen
{
    public:
        // default constructors ok
        virtual void WriteBoys(std::ostream & os) const;

    protected:
        virtual std::vector<std::string> Includes(void) const;
};


struct BoysVRef : public BoysGen
{
    public:
        // default constructors ok
        virtual void WriteBoys(std::ostream & os) const;

    protected:
        virtual std::vector<std::string> Includes(void) const;
};


#endif
