/*! \file
 *
 * \brief Common types used in the generator, plus some helper functions (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__TYPES_HPP_
#define SIMINT_GUARD_GENERATOR__TYPES_HPP_

#include <array>
#include <string>
#include <sstream>
#include <ostream>
#include <vector>
#include <set>
#include <map>
#include <utility>

static const char * amchar = "spdfghijklmnoqrtuvwxyzabceSPDFGHIJKLMNOQRTUVWXYZABCE0123456789";

//! Doublet AM (ie, AM pair for bra or ket)
typedef std::array<int, 2> DAM;

//! Quartet AM (ie, AM pair for bra or ket)
typedef std::array<int, 4> QAM;

//! \brief Exponents on the xyz prefactor of a gaussian
typedef std::array<int, 3> ExpList;

extern
std::map<int, std::vector<ExpList>> gorder_map; // in Ordering.cpp


/*! \brief Get the index of this gaussian in the standard order
 *
 * Ie, xxx = 0, xxy = 1, etc.
 *
 * \param [in] ijk Exponents
 */
int GaussianOrder(const ExpList & ijk);


//! Type of doublet being stored/used
enum class DoubletType
{
    BRA,
    KET
};


/*! \brief Direction of the recurrence relation
 *
 * IJKL = Four angular momenta of a quartet. The letter
 * represents the AM being increased (ie, K means form K+1).
 */
enum class RRStepType
{
    I = 0,
    J = 1,
    K = 2,
    L = 3
};


/*! \brief Cartesian direction of a recurrence step
 */
enum class XYZStep
{
    STEP_X = 0,
    STEP_Y = 1,
    STEP_Z = 2
};


/*! \brief Convert a step to a string
 *
 * \param [in] xyz Step to convert
 * \return Step as a string ("x", "y", "z")
 */
inline std::string XYZStepToStr(XYZStep xyz)
{
    if(xyz == XYZStep::STEP_X)
        return "x";
    else if(xyz == XYZStep::STEP_Y)
        return "y";
    else
        return "z";
}


/*! \brief Write a step to an ostream
 */
inline std::ostream & operator<<(std::ostream & os, const XYZStep xyz)
{
    os << XYZStepToStr(xyz);
    return os;
}


/*! \brief Convert an index (0,1,2) to an XYZ direction
 */
inline XYZStep IdxToXYZStep(int xyz)
{
    if(xyz == 0)
        return XYZStep::STEP_X;
    else if(xyz == 1)
        return XYZStep::STEP_Y;
    else
        return XYZStep::STEP_Z;
}


/*! \brief Convert an XYZ step direction to an index
 */
inline int XYZStepToIdx(XYZStep s)
{
    if(s == XYZStep::STEP_X)
        return 0;
    else if(s == XYZStep::STEP_Y)
        return 1;
    else //(s == XYZStep::STEP_Z)
        return 2;
}


/*! \brief A gaussian center of a particular AM
 *
 * \f[
 * \chi_A = (r_x - A_x)^{i} (r_y-A_y)^{j} (r_z - A_z)^{k} e^{-a r^2}
 * \f]
 *
 * Really we just need the ijk exponents here.
 */
struct Gaussian
{
    //! \brief Exponents on the XYZ prefactor of the gaussian
    ExpList ijk;

    //! Get the AM of the Gaussian center
    int am(void) const { return ijk[0] + ijk[1] + ijk[2]; }


    //! Get index of this gaussian in the standard Gaussian ordering
    int index(void) const { return GaussianOrder(ijk); }


    //! Get the number of cartesian functions in a Gaussian of this AM
    int ncart(void) const { return ((am()+1)*(am()+2))/2; }

    //! Get a string representing this Gaussian
    std::string str(void) const
    {
        std::stringstream ss;

        if(*this)
            ss << amchar[am()] << "_" << ijk[0] << "_" << ijk[1] << "_" << ijk[2];
        else
            ss << "?_" << ijk[0] << ijk[1] << ijk[2];
        return ss.str();
    }


    //! Compare if this Gaussian comes before another in the standard ordering
    bool operator<(const Gaussian & rhs) const
    {
        if(am() < rhs.am())
            return true;
        else if(am() == rhs.am())
            return index() < rhs.index();

        return false;
    }


    //! Compare if this Gaussian is the same as another
    bool operator==(const Gaussian & rhs) const
    {
        return (ijk == rhs.ijk);
    }


    /*! \brief Check if this Gaussian is valid
     *
     * That is, are all ijkl >= 0
     */
    operator bool(void) const
    {
        return (ijk[0] >= 0 && ijk[1] >= 0 && ijk[2] >= 0);
    }


    /*! \brief Make this the next gaussian in the standard ordering
     *
     * \return True if the resulting gaussian is valid
     */
    bool Iterate(void)
    {
        size_t idx = static_cast<size_t>(index()) + 1;
        const auto & v = gorder_map.at(am());

        if(idx >= v.size())
            return false;
        else
        {
            ijk = v.at(idx);
            return *this;
        }
    }


    /*! \brief Increase the exponent of a cartesian direction
     *
     * \param [in] step Which direction to increase
     * \param [in] n Amount to increase by
     * \return New Gaussian identical to this one except the exponent has been increased
     */
    Gaussian StepUp(XYZStep step, int n = 1) const
    {
        return StepUp(XYZStepToIdx(step), n);
    }


    /*! \brief Increase the exponent of a cartesian direction
     *
     * \param [in] idx Index of the cartesian direction to increase
     * \param [in] n Amount to increase by
     * \return New Gaussian identical to this one except the exponent has been increased
     */
    Gaussian StepUp(int idx, int n = 1) const
    {
        Gaussian g(*this);
        g.ijk[idx] += n;
        return g;
    }

    /*! \brief Decrease the exponent of a cartesian direction
     *
     * \param [in] step Which direction to decrease
     * \param [in] n Amount to increase by
     * \return New Gaussian identical to this one except the exponent has been decreased
     */
    Gaussian StepDown(XYZStep step, int n = 1) const
    {
        return StepDown(XYZStepToIdx(step), n);
    }


    /*! \brief Decrease the exponent of a cartesian direction
     *
     * \param [in] idx Index of the cartesian direction to decrease
     * \param [in] n Amount to increase by
     * \return New Gaussian identical to this one except the exponent has been decrease
     */
    Gaussian StepDown(int idx, int n = 1) const
    {
        Gaussian g(*this);
        g.ijk[idx] -= n;
        return g;
    }
};


/*! \brief Output a gaussian to an ostream
 */
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
    int index(void) const { return left.index() * right.ncart() + right.index(); }
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
    int index(void) const { return bra.index() * ket.ncart() + ket.index(); }
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



typedef std::set<std::string> StringSet;
typedef std::set<int> IntSet;

typedef std::set<QAM> QAMSet;
typedef std::set<DAM> DAMSet;

typedef std::vector<QAM> QAMList;
typedef std::vector<DAM> DAMList;

typedef std::set<Quartet> QuartetSet;
typedef std::set<Doublet> DoubletSet;


typedef std::set<Gaussian> GaussianSet;


// Other typedefs
typedef StringSet IncludeSet;
typedef std::map<std::string, std::string> ConstantMap;


// Some helper functions
QuartetSet GenerateInitialQuartetTargets(QAM amlst);

DoubletSet GenerateInitialDoubletTargets(DAM amlst, DoubletType type);

int GaussianOrder(const QAM & ijk);

GaussianSet AllGaussiansForAM(int am);

bool ValidQAM(QAM am);


#endif
