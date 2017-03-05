/*! \file
 *
 * \brief Calculation of number of cartesian integrals for a given angular momentum
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__NCART_HPP_
#define SIMINT_GUARD_GENERATOR__NCART_HPP_

#include "generator/Types.hpp"

#define NCART_(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)
#define NSPH_(am) (2*am+1)

inline int NCART(int am)
{
    return NCART_(am);
}

inline int NCART(int am1, int am2)
{
    return NCART_(am1) * NCART_(am2);
}

inline int NCART(int am1, int am2, int am3)
{
    return NCART_(am1) * NCART_(am2) * NCART_(am3);
}

inline int NCART(int am1, int am2, int am3, int am4)
{
    return NCART_(am1) * NCART_(am2) * NCART_(am3) * NCART_(am4);
}


inline int NCART(DAM am)
{
    int ncart = 1;

    for(auto & it : am)
        ncart *= NCART_(it);

    return ncart;
}

inline int NCART(QAM am)
{
    int ncart = 1;

    for(auto & it : am)
        ncart *= NCART_(it);

    return ncart;
}


inline int NSPH(int am)
{
    return NSPH_(am);
}

inline int NSPH(int am1, int am2)
{
    return NSPH_(am1) * NSPH_(am2);
}

inline int NSPH(int am1, int am2, int am3)
{
    return NSPH_(am1) * NSPH_(am2) * NSPH_(am3);
}

inline int NSPH(int am1, int am2, int am3, int am4)
{
    return NSPH_(am1) * NSPH_(am2) * NSPH_(am3) * NSPH_(am4);
}


inline int NSPH(DAM am)
{
    int ncart = 1;

    for(auto & it : am)
        ncart *= NSPH_(it);

    return ncart;
}

inline int NSPH(QAM am)
{
    int ncart = 1;

    for(auto & it : am)
        ncart *= NSPH_(it);

    return ncart;
}


#endif
