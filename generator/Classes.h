#ifndef GENERATOR_CLASSES_H
#define GENERATOR_CLASSES_H

#include <array>
#include <string>
#include <sstream>
#include <ostream>
#include <vector>
#include <set>

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


inline std::ostream & operator<<(std::ostream & os, const XYZStep xyz)
{
    if(xyz == XYZStep::STEP_X)
        os << "X";
    else if(xyz == XYZStep::STEP_Y)
        os << "Y";
    else if(xyz == XYZStep::STEP_Z)
        os << "Z";
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



struct Gaussian
{
    std::array<int, 3> ijk;

    int am(void) const { return ijk[0] + ijk[1] + ijk[2]; }

    std::string str(void) const
    {
       const char * amchar = "SPDFGHIJKLMNOQRTUVWXYZABCE";
       std::stringstream ss;
       ss << amchar[am()] << "_" << ijk[0] << ijk[1] << ijk[2];
       return ss.str(); 
    }


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

    void Iterate(void)
    {
        if(ijk[2] < (am() - ijk[0]))
            ijk = std::array<int, 3>{ijk[0],   ijk[1]-1,      ijk[2]+1 };
        else
            ijk = std::array<int, 3>{ijk[0]-1, am()-ijk[0]+1, 0        };
    }
};


inline std::ostream & operator<<(std::ostream & os, const Gaussian & g)
{
    os << g.str();
    return os;
}




struct Doublet
{
    DoubletType type;
    Gaussian left;
    Gaussian right;

    int am(void) const { return left.am() + right.am(); }    

    std::string str(void) const
    {
        std::stringstream ss;
        if(type == DoubletType::BRA)
          ss << "(" << left << " " << right << "|";
        else
          ss << "|" << left << " " << right << ")";
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

    std::string str(void) const
    {
        std::stringstream ss;
        ss << "( " << bra.left << " " << bra.right << " | "
                   << ket.left << " " << ket.right << " )";
        return ss.str();
    }


    std::string code_var(void) const
    {
        std::stringstream ss;
        ss << "Q_" 
           << bra.left.ijk[0] << bra.left.ijk[1] << bra.left.ijk[2]      << "_"
           << bra.right.ijk[0] << bra.right.ijk[1] << bra.right.ijk[2]   << "_"
           << ket.left.ijk[0] << ket.left.ijk[1] << ket.left.ijk[2]      << "_"
           << ket.right.ijk[0] << ket.right.ijk[1] << ket.right.ijk[2]   << "_"
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

};

inline std::ostream & operator<<(std::ostream & os, const Quartet & q)
{
    os << q.str();
    return os;
}


struct HRRStep
{
    Quartet target;
    Quartet src1;
    Quartet src2;
    DoubletType steptype;  // bra or ket being stepped
    XYZStep xyz;    


    std::string str(void) const
    {
        const char * xyztype = (steptype == DoubletType::BRA ? "ab" : "cd");
        std::stringstream ss;
        ss << target << " = " << src1 << " + " << xyz << "_" << xyztype << " * " << src2;
        return ss.str();
    }

    bool operator==(const HRRStep & rhs) const
    {
        return (target == rhs.target &&
                src1 == rhs.src1 &&
                src2 == rhs.src2 &&
                steptype == rhs.steptype &&
                xyz == rhs.xyz);
    }


    std::string code_line(void) const
    {
        const char * xyztype = (steptype == DoubletType::BRA ? "ab" : "cd");

        std::stringstream ss;
        ss << "const double " << target.code_var() << " = " << src1.code_var() 
                                                   << " + (" << xyz << "_" << xyztype 
                                                   << " * " << src2.code_var() << ");";
        return ss.str();
    }
};


inline std::ostream & operator<<(std::ostream & os, const HRRStep & hrr)
{
    os << hrr.str();
    return os;
}


typedef std::vector<HRRStep> HRRList;
typedef std::set<Quartet> QuartetSet;



#endif
