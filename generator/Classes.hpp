#ifndef GENERATOR_CLASSES_HPP
#define GENERATOR_CLASSES_HPP

#include <array>
#include <string>
#include <sstream>
#include <ostream>
#include <vector>
#include <set>
#include <map>
#include <utility>

static const char * amchar = "spdfghijklmnoqrtuvwxyzabceSPDFGHIJKLMNOQRTUVWXYZABCE0123456789";

typedef std::array<int, 2> DAM;
typedef std::array<int, 4> QAM;
typedef std::array<int, 3> ExpList;

// In Helpers.cpp
int GaussianOrder(const ExpList & ijk);

enum class DoubletType
{
    BRA,
    KET
};

enum class XYZStep
{
    STEP_X,
    STEP_Y,
    STEP_Z
};

inline std::string XYZStepToStr(XYZStep xyz)
{
    if(xyz == XYZStep::STEP_X)
        return "x";
    else if(xyz == XYZStep::STEP_Y)
        return "y";
    else
        return "z";
    
}

inline std::ostream & operator<<(std::ostream & os, const XYZStep xyz)
{
    os << XYZStepToStr(xyz);
    return os;
}

inline XYZStep IdxToXYZStep(int xyz)
{
    if(xyz == 0) 
        return XYZStep::STEP_X;
    else if(xyz == 1) 
        return XYZStep::STEP_Y;
    else
        return XYZStep::STEP_Z;
}

inline int XYZStepToIdx(XYZStep s)
{
    if(s == XYZStep::STEP_X)
        return 0;
    else if(s == XYZStep::STEP_Y)
        return 1;
    else //(s == XYZStep::STEP_Z)
        return 2;
}




struct Gaussian
{
    ExpList ijk;

    int am(void) const { return ijk[0] + ijk[1] + ijk[2]; }
    int idx(void) const { return GaussianOrder(ijk); } 
    int ncart(void) const { return ((am()+1)*(am()+2))/2; }

    std::string str(void) const
    {
        std::stringstream ss;
        
        if(*this)
            ss << amchar[am()] << "_" << ijk[0] << "_" << ijk[1] << "_" << ijk[2];
        else
            ss << "?_" << ijk[0] << ijk[1] << ijk[2];
        return ss.str(); 
    }


    // notice the reversing of the exponents
    bool operator<(const Gaussian & rhs) const
    {
        if(am() < rhs.am())
            return true;
        else if(am() == rhs.am())
        {
            if(ijk > rhs.ijk)
                return true;
        }

        return false;
    }

    bool operator==(const Gaussian & rhs) const
    {
        return (ijk == rhs.ijk);
    }

    operator bool(void) const
    {
        return (ijk[0] >= 0 && ijk[1] >= 0 && ijk[2] >= 0);
    }
    
    bool Iterate(void)
    {
        if(ijk[2] == am())  // at the end
            return false;

        if(ijk[2] < (am() - ijk[0]))
            ijk = {ijk[0],   ijk[1]-1,      ijk[2]+1 };
        else
            ijk = {ijk[0]-1, am()-ijk[0]+1, 0        };
        return *this;
    }


    Gaussian StepUp(XYZStep step, int n = 1) const
    {
        return StepUp(XYZStepToIdx(step), n);
    }

    Gaussian StepUp(int idx, int n = 1) const
    {
        Gaussian g(*this);
        g.ijk[idx] += n;
        return g;
    }

    Gaussian StepDown(XYZStep step, int n = 1) const
    {
        return StepDown(XYZStepToIdx(step), n);
    }

    Gaussian StepDown(int idx, int n = 1) const
    {
        Gaussian g(*this);
        g.ijk[idx] -= n;
        return g;
    }
};


inline std::ostream & operator<<(std::ostream & os, const Gaussian & g)
{
    os << g.str();
    return os;
}



// A single bra or ket, containing two gaussians
struct Doublet
{
    DoubletType type;
    Gaussian left;
    Gaussian right;

    int am(void) const { return left.am() + right.am(); }    
    int idx(void) const { return left.idx() * right.ncart() + right.idx(); }
    int ncart(void) const { return left.ncart() * right.ncart(); }
    DAM amlist(void) const { return {left.am(), right.am()}; }

    std::string str(void) const
    {
        std::stringstream ss;
        if(type == DoubletType::BRA)
          ss << "(" << left << "  " << right << "|";
        else
          ss << "|" << left << "  " << right << ")";

        return ss.str(); 
    }

    bool operator<(const Doublet & rhs) const
    {
        if(am() < rhs.am())
            return true;
        else if(am() == rhs.am())
        {
            if(left < rhs.left)
                return true;
            else if(left == rhs.left)
            {
                if(right < rhs.right)
                    return true;
            }
        }
        return false;
    }

    bool operator==(const Doublet & rhs) const
    {
        return (left == rhs.left &&
                right == rhs.right &&
                type == rhs.type);
    }

    operator bool(void) const
    {
        return (left && right);
    }
};

inline std::ostream & operator<<(std::ostream & os, const Doublet & d)
{
    os << d.str();
    return os;
}


struct Quartet
{
    Doublet bra;
    Doublet ket;
    int m;

    int am(void) const { return bra.am() + ket.am(); }
    int idx(void) const { return bra.idx() * ket.ncart() + ket.idx(); }
    int ncart(void) const { return bra.ncart() * ket.ncart(); }
    QAM amlist(void) const { return { bra.left.am(), bra.right.am(),
                                          ket.left.am(), ket.right.am() }; }

    Doublet get(DoubletType type) const
    {
        return (type == DoubletType::BRA ? bra : ket);
    }

    std::string str(void) const
    {
        std::stringstream ss;
        ss << "( " << bra.left << "  " << bra.right << " | "
                   << ket.left << "  " << ket.right << " )^"
                   << m;
        return ss.str();
    }


    bool operator<(const Quartet & rhs) const
    {
        if(am() < rhs.am())
            return true;
        else if(am() == rhs.am())
        {
            if(bra < rhs.bra)
                return true;
            else if(bra == rhs.bra)
            {
                if(ket < rhs.ket)
                    return true;
                else if(ket == rhs.ket)
                {
                    if(m < rhs.m)
                        return true;
                }
            }
        }

        return false;
    }

    bool operator==(const Quartet & rhs) const
    {
        return (bra == rhs.bra && ket == rhs.ket && m == rhs.m);
    }

    operator bool(void) const
    {
        return bra && ket && (m >= 0);
    }

};


inline std::ostream & operator<<(std::ostream & os, const Quartet & q)
{
    os << q.str();
    return os;
}


struct HRRDoubletStep
{
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
    char pq;
    DoubletType direction;

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
                direction == rhs.direction &&
                ik == rhs.ik &&
                pq == rhs.pq);
    }
};



inline std::ostream & operator<<(std::ostream & os, const ETStep & et)
{
    os << et.str();
    return os;
}



struct VRRStep
{
    DoubletType type;
    Quartet target;
    std::array<Quartet, 8> src;    
    XYZStep xyz;
    std::array<int, 4> ijkl;

    std::string str(void) const
    {
        return "TODO";
    }

    bool operator==(const VRRStep & rhs) const
    {
        // no need to compare ijkl - it's just to help
        return (target == rhs.target &&
                src == rhs.src &&
                xyz == rhs.xyz &&
                type == rhs.type);
    }
};


typedef std::set<std::string> StringSet;
typedef std::set<int> IntSet;

typedef std::set<QAM> QAMSet;
typedef std::set<DAM> DAMSet;

typedef std::vector<QAM> QAMList;
typedef std::vector<DAM> DAMList;

typedef std::set<Quartet> QuartetSet;
typedef std::set<Doublet> DoubletSet;


typedef std::set<Gaussian> GaussianSet;






#endif
