#ifndef NCART_HPP
#define NCART_HPP

#include <array>

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
