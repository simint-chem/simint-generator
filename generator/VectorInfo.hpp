#ifndef SIMINT_GUARD_GENERATOR__VECTORINFO_HPP_
#define SIMINT_GUARD_GENERATOR__VECTORINFO_HPP_

#include "generator/StringBuilder.hpp"


class VectorInfo
{
public:
    // all default constructors, etc, ok

    // To be implemented by derived classes
    virtual std::string DoubleType(void) const = 0;
    virtual std::string UnionType(void) const = 0;
    virtual std::string DoubleSet1(const std::string & dbl) const = 0;
    virtual std::string DoubleLoad(const std::string & ptr, const std::string & idx) const = 0;
    virtual std::string DoubleStore(const std::string & var, const std::string & ptr, const std::string & idx) const = 0;

    virtual std::string FMAdd(const std::string & a, const std::string & b, const std::string & c) const = 0;
    virtual std::string FMSub(const std::string & a, const std::string & b, const std::string & c) const = 0;

    virtual std::string Sqrt(const std::string & val) const = 0;
    virtual std::string RSqrt(const std::string & val) const = 0;
    virtual std::string Power(const std::string & base, const std::string & exp) const = 0;
    virtual std::string Exp(const std::string & exp) const = 0;

    virtual std::string IntConstant(int i) const = 0;


    // Uses the pure virtual functions
    virtual std::string ConstDoubleType(void) const { return StringBuilder("const ", DoubleType()); }

    virtual std::string NewDoubleSet1(const std::string & var, const std::string & val) const
    {
        return StringBuilder(DoubleType(), " ", var, " = ", DoubleSet1(val));
    }


    virtual std::string NewDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx) const
    {
        return StringBuilder(DoubleType(), " ", var, " = ", DoubleLoad(ptr, idx));
    }

    virtual std::string NewConstDoubleSet1(const std::string & var, const std::string & val) const
    {
        return StringBuilder("const ", NewDoubleSet1(var, val));
    }


    virtual std::string NewConstDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx) const
    {
        return StringBuilder("const ", NewDoubleLoad(var, ptr, idx));
    }

    virtual std::string ConstUnionType(void) const
    {
        return StringBuilder("const ", UnionType());
    }
};




class ScalarVector : public VectorInfo
{
public:
    virtual std::string DoubleType(void) const { return "double"; }

    virtual std::string UnionType(void) const { return ""; }

    virtual std::string DoubleSet1(const std::string & dbl) const { return dbl; }

    virtual std::string DoubleLoad(const std::string & ptr, const std::string & idx) const { return ptr + "[" + idx + "]"; }

    virtual std::string DoubleStore(const std::string & var, const std::string & ptr, const std::string & idx) const { return ptr + "[" + idx + "] = " + var; }

    virtual std::string FMAdd(const std::string & a, const std::string & b, const std::string & c) const
    {
        return StringBuilder("fma(", a, ", ", b, ", ", c, ")");
    }

    virtual std::string FMSub(const std::string & a, const std::string & b, const std::string & c) const
    {
        return StringBuilder("fma(", a, ", ", b, ", -", c, ")");
    }

    virtual std::string Sqrt(const std::string & val) const
    {
        return StringBuilder("sqrt(", val, ")");
    }

    virtual std::string RSqrt(const std::string & val) const
    {
        return StringBuilder("(1.0 / sqrt(", val, "))");
    }

    virtual std::string Power(const std::string & base, const std::string & exp) const
    {
        return StringBuilder("pow(", base, ",", exp, ")");
    }

    virtual std::string Exp(const std::string & exp) const
    {
        return StringBuilder("exp(", exp, ")");
    }

    virtual std::string IntConstant(int i) const
    {
        return std::to_string(i);
    }

};



class BasicIntelSIMDVector : public VectorInfo
{
public:
    BasicIntelSIMDVector(int width)
        : width_(width), dwidth_(width/64)
    {}


    virtual std::string DoubleType(void) const { return StringBuilder("__m", width_, "d"); }

    virtual std::string UnionType(void) const { return StringBuilder("union double", dwidth_); }

    virtual std::string DoubleSet1(const std::string & dbl) const
    {
        return StringBuilder("_mm", width_, "_set1_pd(" + dbl + ")");
    }

    virtual std::string DoubleLoad(const std::string & ptr, const std::string & idx) const
    {
        std::string ptrstr(ptr);
        if(idx.size())
            ptrstr += " + " + idx;
        return StringBuilder("_mm", width_, "_load_pd(", ptrstr, ")");
    }

    virtual std::string DoubleStore(const std::string & var, const std::string & ptr, const std::string & idx) const
    {
        std::string ptrstr(ptr);
        if(idx.size())
            ptrstr += " + " + idx;
        return StringBuilder("_mm", width_, "_store_pd(", ptrstr, ", ", var, ")");
    }

    virtual std::string FMAdd(const std::string & a, const std::string & b, const std::string & c) const
    {
        return StringBuilder("_mm", width_, "_fmadd_pd(", a, ", ", b, ", ", c, ")");
    }

    virtual std::string FMSub(const std::string & a, const std::string & b, const std::string & c) const
    {
        return StringBuilder("_mm", width_, "_fmsub_pd(", a, ", ", b, ", ", c, ")");
    }

    virtual std::string Sqrt(const std::string & val) const
    {
        return StringBuilder("_mm", width_, "_sqrt_pd(", val, ")");
    }

    virtual std::string RSqrt(const std::string & val) const
    {
        return StringBuilder("(1.0/", Sqrt(val), ")");
    }

    virtual std::string Power(const std::string & base, const std::string & exp) const
    {
        return StringBuilder("_mm", width_, "_pow_pd(", base, ",", exp, ")");
    }

    virtual std::string Exp(const std::string & exp) const
    {
        return StringBuilder("_mm", width_, "_exp_pd(", exp, ")");
    }

    virtual std::string IntConstant(int i) const
    {
        return StringBuilder("const_", i);
    }

private:
    int width_;
    int dwidth_;
};






#endif
