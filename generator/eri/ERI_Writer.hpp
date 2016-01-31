#ifndef SIMINT_GUARD_ERI__ERI_WRITER_HPP_
#define SIMINT_GUARD_ERI__ERI_WRITER_HPP_

#include <ostream>


class BoysGenerator;
class ERI_VRR_Writer;
class ERI_ET_Writer;
class ERI_HRR_Writer;
class ERIGeneratorInfo;
class VectorInfo;


class ERI_Writer
{
public:
    ERI_Writer(std::ostream & os,
               std::ostream & osh,
               const ERIGeneratorInfo & info,
               const BoysGenerator & bg,
               const ERI_VRR_Writer & vrr_writer,
               const ERI_ET_Writer & et_writer,
               const ERI_HRR_Writer & hrr_writer);


    ERI_Writer(const ERI_Writer &) = default;
    ERI_Writer(ERI_Writer &&) = default;
    ERI_Writer & operator=(const ERI_Writer &) = default;
    ERI_Writer & operator=(ERI_Writer &&) = default;



    virtual void WriteFile(void) const = 0;


protected:
    std::ostream & os_;
    std::ostream & osh_;
    const ERIGeneratorInfo & info_;
    const VectorInfo & vinfo_;
    const BoysGenerator & bg_;
    const ERI_VRR_Writer & vrr_writer_;
    const ERI_ET_Writer & et_writer_;
    const ERI_HRR_Writer & hrr_writer_;

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

    void WriteFile_Permute_(void) const;
    void WriteFile_NoPermute_(void) const;
};


#endif
