/*! \file
 *
 * \brief Calculation of number of cartesian integrals for a given angular momentum
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__NCART_HPP_
#define SIMINT_GUARD_GENERATOR__NCART_HPP_

#include "generator/Types.hpp"

#define NCART_(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)

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


#endif
