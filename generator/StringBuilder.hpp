/*! \file
 *
 * \brief String creation
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__STRINGBUILDER_HPP_
#define SIMINT_GUARD_GENERATOR__STRINGBUILDER_HPP_

#include <string>
#include <sstream>

// by default, use std::to_string
template<typename T>
std::string StringBuilder(const T & arg)
{
    return std::to_string(arg);
}


// overload for a double
inline std::string StringBuilder(const double & arg)
{
    std::stringstream ss;
    ss.precision(18);
    ss << arg;
    return ss.str();
}


// overload for string, since std::to_string(std::string) doesn't exist
inline std::string StringBuilder(const std::string & arg)
{
    return arg;
}

// overload for const char *, since std::to_string(const char *) doesn't exist
inline std::string StringBuilder(const char * arg)
{
    return std::string(arg);
}

// overload for char, so it doesn't get converted to an integer type
inline std::string StringBuilder(const char & arg)
{
    return std::string(1, arg);
}


/*! \brief Concatenates all arguments to form a string
 *
 * All arguments are converted to a string type first (through std::to_string or something
 * similar
 */
template<typename T, typename ... Targs>
std::string StringBuilder(const T & arg, Targs && ... args)
{
    return StringBuilder(arg) + StringBuilder(std::forward<Targs>(args)...);
}


#endif
