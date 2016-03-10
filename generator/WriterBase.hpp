/*! \file
 *
 * \brief Base class for output writers
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__WRITERBASE_HPP_
#define SIMINT_GUARD_GENERATOR__WRITERBASE_HPP_

#include <map>
#include <set>
#include <string>


class WriterBase
{
public:
    virtual bool IsInline(void) const = 0;
    virtual bool IsExternal(void) const = 0;
};


#endif
