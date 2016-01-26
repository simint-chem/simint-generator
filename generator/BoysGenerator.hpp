/*! \file
 *
 * \brief Classes for generating the Boys function
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__BOYSGENERATOR_HPP_
#define SIMINT_GUARD_GENERATOR__BOYSGENERATOR_HPP_

#include <vector>

#include "generator/Types.hpp" // for ConstantMap and IncludeSet

// Forward declarations
class GeneratorInfoBase;
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
class BoysGenerator
{
    public:
        BoysGenerator(const GeneratorInfoBase & info);

        virtual ~BoysGenerator() = default;

        virtual ConstantMap GetConstants(void) const;
        virtual IncludeSet GetIncludes(void) const = 0;
        virtual void WriteBoys(std::ostream & os) const = 0;

    protected:
        const GeneratorInfoBase & info_;
        const VectorInfo & vinfo_;
};




class BoysFO : public BoysGenerator
{
    public:
        BoysFO(const GeneratorInfoBase & info, std::string dir); // read from directory

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


class BoysSplit : public BoysGenerator
{
    public:
        using BoysGenerator::BoysGenerator;

        virtual ConstantMap GetConstants(void) const;
        virtual IncludeSet GetIncludes(void) const;
        virtual void WriteBoys(std::ostream & os) const;
};




#endif
