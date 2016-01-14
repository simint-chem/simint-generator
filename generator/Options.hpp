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


//! \brief free function to create a map with default options
inline OptionMap DefaultOptions(void)
{
    return OptionMap{
                      {Option::StackMem, 0},
                      {Option::InlineVRR, 1},
                      {Option::InlineET, 1},
                      {Option::InlineHRR, 1},
                      {Option::Scalar, 0},
                      {Option::NoSingleET, 0},
                      {Option::NoET, 1}
                    };
}



#endif
