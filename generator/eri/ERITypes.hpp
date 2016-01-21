#ifndef SIMINT_GUARD_ERI__ERITYPES_HPP_
#define SIMINT_GUARD_ERI__ERITYPES_HPP_

#include "generator/Types.hpp"

struct HRRDoubletStep
{
    RRStepType steptype;
    Doublet target;
    std::array<Doublet, 2> src;
    XYZStep xyz;    

    std::string str(void) const
    {
        const char * xyztype = (target.type == DoubletType::BRA ? "ab" : "cd");
        std::stringstream ss;
        ss << target << " = " << src[0] << " + " << xyz << "_" << xyztype << " * " << src[1];
        return ss.str();
    }

    bool operator==(const HRRDoubletStep & rhs) const
    {
        return (target == rhs.target &&
                src == rhs.src &&
                xyz == rhs.xyz);
    }
    
    bool operator<(const HRRDoubletStep & rhs) const
    {
        return this->target < rhs.target;
    }
};



inline std::ostream & operator<<(std::ostream & os, const HRRDoubletStep & hrr)
{
    os << hrr.str();
    return os;
}


struct ETStep
{
    Quartet target;
    std::array<Quartet, 4> src;
    std::array<int, 2> ik;
    XYZStep xyz;

    std::string str(void) const
    {
        std::stringstream ss;
        ss << target << " = " << xyz << " * " << src[0];
        if(src[1].bra.left && src[1].ket.left)
            ss << " + " << src[1];
        if(src[2].bra.left && src[2].ket.left)
            ss << " + " << src[2];
        if(src[3].bra.left && src[3].ket.left)
            ss << " - " << src[3];
 
        return ss.str();
    }

    bool operator==(const ETStep & rhs) const
    {
        return (target == rhs.target &&
                src == rhs.src &&
                xyz == rhs.xyz &&
                ik == rhs.ik);
    }

    bool operator<(const ETStep & rhs) const
    {
        return this->target < rhs.target;
    }
};



inline std::ostream & operator<<(std::ostream & os, const ETStep & et)
{
    os << et.str();
    return os;
}



struct VRRStep
{
    RRStepType type;
    Quartet target;
    std::array<Quartet, 8> src;    
    XYZStep xyz;
    std::array<int, 4> ijkl;

    std::string str(void) const
    {
        std::stringstream ss;
        ss << target << " = ";
        for(const auto & it  : src)
        {
            if(it)
                ss << " + " << it;
        }

        return ss.str();
    }

    bool operator==(const VRRStep & rhs) const
    {
        // no need to compare ijkl - it's just to help
        return (target == rhs.target &&
                src == rhs.src &&
                xyz == rhs.xyz &&
                type == rhs.type);
    }

    bool operator<(const VRRStep & rhs) const
    {
        return this->target < rhs.target;
    }
};

inline std::ostream & operator<<(std::ostream & os, const VRRStep & vrr)
{
    os << vrr.str();
    return os;
}


#endif
