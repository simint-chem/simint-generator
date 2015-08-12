#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <map>

typedef std::map<int, int> OptionsMap;


// options
#define OPTION_STACKMEM   1
#define OPTION_INLINEVRR  2
#define OPTION_INLINEET   3
#define OPTION_INLINEHRR  4
#define OPTION_INTRINSICS 5
#define OPTION_SCALAR     6

#define OPTION_NOSINGLEET 7
#endif
