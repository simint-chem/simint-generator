#ifndef STRINGBUILDER_HPP
#define STRINGBUILDER_HPP

#include <string>
#include <sstream>

template<typename T>
std::string StringBuilder(const T & arg)
{
    return std::to_string(arg);
}

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

// overload for char, so it doesn't get converted to an integer type
inline std::string StringBuilder(const char & arg)
{
    return std::string(1, arg);
}

// overload for string, since std::to_string(const char *) doesn't exist
inline std::string StringBuilder(const char * arg)
{
    return std::string(arg);
}


template<typename T, typename ... Targs>
std::string StringBuilder(const T & arg, Targs && ... args)
{
    return StringBuilder(arg) + StringBuilder(std::forward<Targs>(args)...);
}


#endif
