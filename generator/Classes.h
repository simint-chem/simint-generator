#ifndef GENERATOR_CLASSES_H
#define GENERATOR_CLASSES_H

#include <array>
#include <string>
#include <sstream>


enum class DoubletType
{
    BRA,
    KET
};


const char * amchar = "SPDFGHIJKLMNOQRTUVWXYZABCE";


struct Gaussian
{
    std::array<int, 3> ijk;

    int am(void) const { return ijk[0] + ijk[1] + ijk[2]; }

    std::string str(void) const
    {
       std::stringstream ss;
       ss << amchar[am()] << "_" << ijk[0] << ijk[1] << ijk[2];
       return ss.str(); 
    }


    bool operator<(const Gaussian & rhs)
    {
        if(am() < rhs.am())
            return true;
        else if(ijk < rhs.ijk)
            return true;
        else
            return false;
    }

    bool operator==(const Gaussian & rhs)
    {
        return (ijk == rhs.ijk);
    }
};



struct Doublet
{
    DoubletType type;
    Gaussian left;
    Gaussian right;

    int am(void) { return left.am() + right.am(); }    

    std::string str(void) const
    {
        std::stringstream ss;
        if(type == DoubletType::BRA)
          ss << "(" << left.str() << " " << right.str() << "|";
        else
          ss << "|" << left.str() << " " << right.str() << ")";
        return ss.str(); 
    }

    bool operator<(const Doublet & rhs)
    {
        if(left < rhs.left)
            return true;
        else if(left == rhs.left)
            if(right < rhs.right)
                return true;
        return false;
    }

    bool operator==(const Doublet & rhs)
    {
        return (left == rhs.left && right == rhs.right);
    }
};


#endif
