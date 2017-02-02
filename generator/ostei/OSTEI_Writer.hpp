#pragma once


#include <ostream>


class OSTEI_VRR_Writer;
class OSTEI_ET_Writer;
class OSTEI_HRR_Writer;
class OSTEI_GeneratorInfo;
class VectorInfo;


class OSTEI_Writer_Base
{
public:
    OSTEI_Writer_Base(std::ostream & os,
                      std::ostream & osh,
                      const OSTEI_GeneratorInfo & info,
                      const OSTEI_VRR_Writer & vrr_writer,
                      const OSTEI_HRR_Writer & hrr_writer);


    OSTEI_Writer_Base(const OSTEI_Writer_Base &) = default;
    OSTEI_Writer_Base(OSTEI_Writer_Base &&) = default;
    OSTEI_Writer_Base & operator=(const OSTEI_Writer_Base &) = default;
    OSTEI_Writer_Base & operator=(OSTEI_Writer_Base &&) = default;

    void DeclarePrimPointers(void) const;

    virtual void WriteFile(void) const = 0;


protected:
    std::ostream & os_;
    std::ostream & osh_;
    const OSTEI_GeneratorInfo & info_;
    const OSTEI_VRR_Writer & vrr_writer_;
    const OSTEI_HRR_Writer & hrr_writer_;
};




class OSTEI_Writer : public OSTEI_Writer_Base
{
public:
    using OSTEI_Writer_Base::OSTEI_Writer_Base;
    using OSTEI_Writer_Base::operator=;

    virtual void WriteFile(void) const;

private:
    void DeclareContwork(void) const;
    void WriteShellOffsets(void) const;
    void WriteAccumulation(void) const;

    std::string FunctionName_(QAM am) const;
    std::string FunctionPrototype_(QAM am) const;

    bool IsSpecialPermutation_(QAM am) const;
    void Write_Full_(void) const;
    void Write_Permutations_(void) const;
    void Write_Permute_(QAM am, bool swap12, bool swap34) const;
};



class OSTEIDeriv1_Writer : public OSTEI_Writer_Base
{
public:
    using OSTEI_Writer_Base::OSTEI_Writer_Base;
    using OSTEI_Writer_Base::operator=;

    virtual void WriteFile(void) const;

private:
    void DeclareContwork(void) const;
    void WriteShellOffsets(void) const;
    void WriteAccumulation(void) const;

    std::string FunctionName_(QAM am) const;
    std::string FunctionPrototype_(QAM am) const;

    bool IsSpecialPermutation_(QAM am) const;
    void Write_Full_(void) const;
    void Write_Permutations_(void) const;
    void Write_Permute_(QAM am, bool swap12, bool swap34) const;
};



