#ifndef VRRWRITER_HPP
#define VRRWRITER_HPP

#include <iostream>

#include "generator/Classes.hpp"
#include "generator/VRR_Algorithm_Base.hpp"

// foward declare
class ERIGeneratorInfo;



class VRR_Writer
{   
    public:
        VRR_Writer(const VRR_Algorithm_Base & vrr_algo);

        bool HasVRR(void) const;
        bool HasBraVRR(void) const;
        bool HasKetVRR(void) const;

        void DeclarePrimArrays(std::ostream & os, const ERIGeneratorInfo & info) const;
        void DeclarePrimPointers(std::ostream & os, const ERIGeneratorInfo & info) const;


        virtual void AddConstants(ERIGeneratorInfo & info) const = 0;
        virtual void WriteVRR(std::ostream & os, const ERIGeneratorInfo & info) const = 0;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const = 0;

    protected:
        void WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepSet & vs, const std::string & num_n, const ERIGeneratorInfo & info) const;

        const VRR_Algorithm_Base & vrr_algo_;
};



class VRR_Writer_Inline : public VRR_Writer
{
    public:
        VRR_Writer_Inline(const VRR_Algorithm_Base & vrr_algo);

        virtual void WriteVRR(std::ostream & os, const ERIGeneratorInfo & info) const;
        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const;
};


class VRR_Writer_External : public VRR_Writer
{
    public:
        VRR_Writer_External(const VRR_Algorithm_Base & vrr_algo);

        virtual void WriteVRR(std::ostream & os, const ERIGeneratorInfo & info) const;
        virtual void AddConstants(ERIGeneratorInfo & info) const;
        virtual void WriteVRRFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const;
};



#endif
