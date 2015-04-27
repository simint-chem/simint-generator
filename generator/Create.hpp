#ifndef CREATE_HPP
#define CREATE_HPP

#include <memory>
#include <array>

class HRR_Algorithm_Base;

void Create_Unrolled(std::array<int, 4> amlist,
                     std::unique_ptr<HRR_Algorithm_Base> & hrralgo,
                     std::ostream & out);

void Create_Looped(std::array<int, 4> amlist,
                   std::unique_ptr<HRR_Algorithm_Base> & hrralgo,
                   std::ostream & out);

#endif
