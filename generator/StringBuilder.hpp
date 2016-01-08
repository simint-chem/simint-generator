#ifndef STRINGBUILDER_HPP
#define STRINGBUILDER_HPP

#include <string>
#include <sstream>

template<typename T>
std::string StringBuilder(const T & arg)
{
    return std::to_string(arg);
}

std::string StringBuilder(const double & arg)
{
    std::stringstream ss;
    ss.precision(18);
    ss << arg;
    return ss.str();
}

// overload for string, since std::to_string(std::string) doesn't exist
std::string StringBuilder(const std::string & arg)
{
    return arg;
}


template<typename T, typename ... Targs>
std::string StringBuilder(const T & arg, Targs && ... args)
{
    return StringBuilder(arg) + StringBuilder(std::forward<Targs>(args)...);
}


#endif
