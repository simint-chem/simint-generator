/*! \file
 *
 * \brief Classes for generating the Boys function
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef BOYS_HPP
#define BOYS_HPP

#include <string>
#include <vector>
#include <map>

// Forward declarations
class ERIGeneratorInfo;


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
        virtual void AddConstants(ERIGeneratorInfo & info) const = 0;
        virtual void AddIncludes(ERIGeneratorInfo & info) const = 0;
        virtual void WriteBoys(std::ostream & os, const ERIGeneratorInfo & info) const = 0;
};




class BoysFO : public BoysGen
{
    public:
        BoysFO(std::string dir); // read from directory

        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void AddIncludes(ERIGeneratorInfo & info) const;
        virtual void WriteBoys(std::ostream & os, const ERIGeneratorInfo & info) const;

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

        void WriteBoysSingle_(std::ostream & os, int m, bool prefac, const ERIGeneratorInfo & info) const;
};


struct BoysSplit : public BoysGen
{
    public:
        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void AddIncludes(ERIGeneratorInfo & info) const;
        virtual void WriteBoys(std::ostream & os, const ERIGeneratorInfo & info) const;
};




/*
struct BoysVRef : public BoysGen
{
    public:
        // default constructors ok
        virtual void WriteBoys(std::ostream & os) const;

        virtual void AddIncludes(ERIGeneratorInfo & info) const;

};
*/

#endif
