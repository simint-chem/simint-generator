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
    void ZeroContwork(void) const;
    void FreeContwork(void) const;
    void WriteShellOffsets(void) const;
    void WriteShellOffsets_Scalar(void) const;
    void WriteAccumulation(void) const;

    std::string FunctionName_(QAM am) const;
    std::string FunctionPrototype_(QAM am) const;

    void WriteFile_Full_(void) const;
    void WriteFile_SpecialPermute_(void) const;

    void WriteFile_Permutations_(void) const;
    void WriteFile_SinglePermutation_(void) const;
};



class OSTEIDeriv1_Writer : public OSTEI_Writer_Base
{
public:
    using OSTEI_Writer_Base::OSTEI_Writer_Base;
    using OSTEI_Writer_Base::operator=;

    virtual void WriteFile(void) const;

private:
    void DeclareContwork(void) const;
    void ZeroContwork(void) const;
    void FreeContwork(void) const;
    void WriteShellOffsets(void) const;
    void WriteShellOffsets_Scalar(void) const;
    void WriteAccumulation(void) const;

    void WriteFile_Permute_(void) const;
    void WriteFile_NoPermute_(void) const;
};



