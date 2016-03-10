/*! \file
 *
 * \brief Concatenation of different types to form a string
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__STRINGBUILDER_HPP_
#define SIMINT_GUARD_GENERATOR__STRINGBUILDER_HPP_

#include <string>
#include <sstream>

/*! \brief Concatenates all arguments to form a string
 *
 * Used for a single argument that can be converted via std::to_string
 */
template<typename T>
std::string StringBuilder(const T & arg)
{
    return std::to_string(arg);
}


/*! \brief Concatenates all arguments to form a string
 *
 * Used for a double type
 */
inline std::string StringBuilder(const double & arg)
{
    std::stringstream ss;
    ss.precision(18);
    ss << arg;
    return ss.str();
}


/*! \brief Concatenates all arguments to form a string
 *
 * Used for a std::string, since no std::to_string exists for that
 */
inline std::string StringBuilder(const std::string & arg)
{
    return arg;
}


/*! \brief Concatenates all arguments to form a string
 *
 * Used for a const char *, since no std::to_string exists for that
 */
inline std::string StringBuilder(const char * arg)
{
    return std::string(arg);
}

/*! \brief Concatenates all arguments to form a string
 *
 * Used for a single character, so it doesn't get converted to an integer type
 */
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
