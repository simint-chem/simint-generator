#ifndef SIMINT_GUARD_GENERATOR__ERI_WRITER_HPP_
#define SIMINT_GUARD_GENERATOR__ERI_WRITER_HPP_

#include <ostream>


class BoysGen;
class VRR_Writer;
class ET_Writer;
class HRR_Writer;
class ERIGeneratorInfo;
class VectorInfo;


class ERI_Writer
{
public:
    ERI_Writer(std::ostream & os,
               const ERIGeneratorInfo & info,
               const BoysGen & bg,
               const VRR_Writer & vrr_writer,
               const ET_Writer & et_writer,
               const HRR_Writer & hrr_writer);


    ERI_Writer(const ERI_Writer &) = default;
    ERI_Writer(ERI_Writer &&) = default;
    ERI_Writer & operator=(const ERI_Writer &) = default;
    ERI_Writer & operator=(ERI_Writer &&) = default;



    virtual void WriteFile(void) const = 0;


protected:
    std::ostream & os_;
    const ERIGeneratorInfo & info_;
    const VectorInfo & vinfo_;
    const BoysGen & bg_;
    const VRR_Writer & vrr_writer_;
    const ET_Writer & et_writer_;
    const HRR_Writer & hrr_writer_;

};




class ERI_Writer_Basic : public ERI_Writer
{
public:
    using ERI_Writer::ERI_Writer;
    using ERI_Writer::operator=;

    virtual void WriteFile(void) const;

private:
    void DeclareContwork(void) const;
    void ZeroContwork(void) const;
    void FreeContwork(void) const;
    void WriteShellOffsets(void) const;
    void WriteShellOffsets_Scalar(void) const;
    void WriteAccumulation(void) const;
    void WriteFile_Permute(void) const;
};


#endif
