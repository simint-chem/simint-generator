#ifndef CREATE_HPP
#define CREATE_HPP

#include <memory>
#include <array>

class HRR_Algorithm_Base;
typedef std::array<int, 4> QAMList;

HRRBraKetStepList Create_DoubletStepLists(QAMList amlist, std::unique_ptr<HRR_Algorithm_Base> & hrralgo);

HRRQuartetStepList Create_QuartetStepList(QAMList amlist, std::unique_ptr<HRR_Algorithm_Base> & hrralg);

#endif
