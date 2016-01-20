/*! \file
 *
 * \brief Classes for generating the Boys function
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef BOYS_HPP
#define BOYS_HPP

#include <vector>

#include "generator/WriterBase.hpp" // for ConstantMap and IncludeSet

// Forward declarations
class ERIGeneratorInfo;
class VectorInfo;


/////////////////////////////////////////////////////////////
// Where to switch from calculating all the boys parameters
// to calculating only the highest and then using the
// downard recurrance
/////////////////////////////////////////////////////////////
#define BOYS_FO_RECUR 3
#define BOYS_SPLIT_RECUR 500 // effectively disabled




/*! \brief Base class for boys function generation
 */
class BoysGen
{
    public:
        BoysGen(const ERIGeneratorInfo & info);

        virtual ConstantMap GetConstants(void) const;
        virtual IncludeSet GetIncludes(void) const = 0;
        virtual void WriteBoys(std::ostream & os) const = 0;

    protected:
        const ERIGeneratorInfo & info_;
        const VectorInfo & vinfo_;
};




class BoysFO : public BoysGen
{
    public:
        BoysFO(const ERIGeneratorInfo & info, std::string dir); // read from directory

        virtual ConstantMap GetConstants(void) const;
        virtual IncludeSet GetIncludes(void) const;
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

        std::string GetFOConstant_(std::string ab, int m, int i) const;

        void WriteBoysSingle_(std::ostream & os, int m, bool prefac) const;
};


class BoysSplit : public BoysGen
{
    public:
        using BoysGen::BoysGen;

        virtual ConstantMap GetConstants(void) const;
        virtual IncludeSet GetIncludes(void) const;
        virtual void WriteBoys(std::ostream & os) const;
};




#endif
