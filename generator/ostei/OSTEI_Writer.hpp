#pragma once


#include <ostream>


class OSTEI_VRR_Writer;
class OSTEI_ET_Writer;
class OSTEI_HRR_Writer;
class OSTEI_GeneratorInfo;
class VectorInfo;


class OSTEI_Writer
{
public:
    OSTEI_Writer(std::ostream & os,
               std::ostream & osh,
               const OSTEI_GeneratorInfo & info,
               const OSTEI_VRR_Writer & vrr_writer,
               const OSTEI_ET_Writer & et_writer,
               const OSTEI_HRR_Writer & hrr_writer);


    OSTEI_Writer(const OSTEI_Writer &) = default;
    OSTEI_Writer(OSTEI_Writer &&) = default;
    OSTEI_Writer & operator=(const OSTEI_Writer &) = default;
    OSTEI_Writer & operator=(OSTEI_Writer &&) = default;



    virtual void WriteFile(void) const = 0;


protected:
    std::ostream & os_;
    std::ostream & osh_;
    const OSTEI_GeneratorInfo & info_;
    const VectorInfo & vinfo_;
    const OSTEI_VRR_Writer & vrr_writer_;
    const OSTEI_ET_Writer & et_writer_;
    const OSTEI_HRR_Writer & hrr_writer_;

};




class OSTEI_Writer_Basic : public OSTEI_Writer
{
public:
    using OSTEI_Writer::OSTEI_Writer;
    using OSTEI_Writer::operator=;

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



