#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <map>

enum class Option
{
    StackMem,
    InlineVRR,
    InlineET,
    InlineHRR,
    Scalar,
    NoSingleET,
    NoET
};


enum class Compiler
{
    Intel,
    GCC
};



typedef std::map<Option, int> OptionMap;


#endif
