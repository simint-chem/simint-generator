#include <stdexcept>

#include "generator/BoysGenerator.hpp"
#include "generator/ERIGeneratorInfo.hpp"
#include "generator/Naming.hpp"
#include "generator/Printing.hpp"
        

////////////////////////////////////
// Boys Split
////////////////////////////////////

ConstantMap BoysSplit::GetConstants(void) const
{
    ConstantMap cm;

    const int L = info_.L();

    cm.emplace("const_1", "1");

    if(L >= BOYS_SPLIT_RECUR)
    {
        // add the factor of 2.0
        cm.emplace("const_2", "2");

        // constants for the recursion
        // skipping i = 0 since that is 1.0
        for(int i = L-1; i > 0; i--)
            cm.emplace(StringBuilder("FO_RECUR_", i),
                       StringBuilder((1.0/(2.0*i+1.0))));
    }

    return cm;
       
}

IncludeSet BoysSplit::GetIncludes(void) const
{
    return {"\"boys/boys_split.h\""};
}

void BoysSplit::WriteBoys(std::ostream & os) const
{
    const int L = info_.L();

    std::string primname = PrimVarName({0,0,0,0});

    if(L < BOYS_SPLIT_RECUR)
    {
        if(info_.Scalar())
            os << indent5 << "Boys_F_split((double *)" << primname << ", " 
                          << L << ", F_x);\n";
        else
            os << indent5 << "Boys_F_split_simd((double *)" << primname << ", " 
                          << L << ", (const double *)(&F_x));\n";

        os << indent5 << vinfo_.ConstDoubleType() << " prefac = " << vinfo_.Sqrt("one_over_PQalpha_sum") << " * P_prefac * " << vinfo_.DoubleLoad("Q.prefac", "j") << ";\n";

        for(int i = 0; i <= L; i++)
            os << indent5 << primname << "[" << i << "] *= prefac;\n";
    }
    else
    {
        // calculate highest, and recurse down
        os << indent5 << "Boys_F_split_single((double *)(" << primname << "+" << L << "), " << L << ", (const double *)(&F_x));\n";


        // calculate the downward recursion factors
        os << indent5 << vinfo_.ConstDoubleType() << " x2 = " << vinfo_.IntConstant(2) << " * F_x;\n";
        os << indent5 << vinfo_.ConstDoubleType() << " ex = " << vinfo_.Exp("-F_x") << ";\n";

        for(int m = L-1; m > 0; m--)
        {
            std::string cname = StringBuilder("FO_RECUR_", m);
            os << indent5 << primname << "[" << m << "] = (x2 * " << primname << "[" << (m+1) << "] + ex) * " << cname << ";\n";
        }

        // do m = 0
        os << indent5 << primname << "[0] = (x2 * " << primname << "[1] + ex);\n"; // times 1.0


        // add prefac now
        os << indent5 << vinfo_.ConstDoubleType() << " prefac = " << vinfo_.Sqrt("one_over_PQalpha_sum") << " * P_prefac * " << vinfo_.DoubleLoad("Q.prefac", "j") << ";\n";
 
        os << "\n";
        os << indent5 << "for(n = 0; n <= " << L << "; ++n)\n";
        os << indent6 << primname << "[n] *= prefac;\n";
        os << "\n";
    }

}

