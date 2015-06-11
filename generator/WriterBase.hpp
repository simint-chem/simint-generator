#ifndef WRITERBASE_H
#define WRITERBASE_H

#include "generator/Options.hpp"
#include "generator/Classes.hpp"

class WriterBase
{
    public:
        WriterBase(const OptionsMap & options, const std::string & prefix, const QAM & finalam);

        void ReadCPUFlags(const std::string & file);

        void SetContQ(const QAMSet & topquartets);

        void WriteIncludes(std::ostream & os) const;

        void DeclareContwork(std::ostream & os) const;
        void ZeroContWork(std::ostream & os, const std::string & nshell) const;
        void FreeContwork(std::ostream & os) const;

        int GetOption(int option) const;
        bool HasCPUFlag(const std::string & flag) const;
        bool Intrinsics(void) const;

        bool IsContArray(const QAM & am) const;
        
        QAM FinalAM(void) const;
        const std::string & Prefix(void) const;

        bool IsFinalAM(const QAM & am) const;
        //bool Permute(void) const;

        bool HasVRR(void) const;
        bool HasET(void) const;
        bool HasHRR(void) const;
        bool HasBraHRR(void) const;
        bool HasKetHRR(void) const;

        size_t MemoryReq(void) const;

        int L(void) const;

        int SimdLen(void) const;
        std::string DoubleType(void) const;
        std::string ConstDoubleType(void) const;
        std::string DoubleSet(const std::string & dbl) const;
        std::string DoubleLoad(const std::string & ptr, const std::string & idx) const;
        std::string DoubleStore(const std::string & var, const std::string & ptr, const std::string & idx) const;
        std::string NewDoubleSet(const std::string & var, const std::string & val) const;
        std::string NewDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx) const;
        std::string NewConstDoubleSet(const std::string & var, const std::string & val) const;
        std::string NewConstDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx) const;

        std::string Sqrt(const std::string & val) const;
        std::string RSqrt(const std::string & val) const;
        std::string Power(const std::string & base, const std::string & exp) const;
        std::string Exp(const std::string & exp) const;

        //void PermuteResult(std::ostream & os, const std::string & src) const;

        static std::string ArrVarName(const QAM & am, const std::string & prefix = "");
        static std::string ArrVarName(int am1, int am2, const std::string & ketstr, const std::string & prefix = "");
        static std::string ArrVarName(const std::string & brastr, int am3, int am4, const std::string & prefix = "");

        std::string HRRVarName(const QAM & am);
        std::string ArrVarName(int am1, int am2, const std::string & ketstr);
        std::string ArrVarName(const std::string & brastr, int am3, int am4);

        static std::string PrimVarName(const QAM & am);
        static std::string PrimPtrName(const QAM & am);

    private:
        QAMSet contq_;  // set of contracted integral AM
        std::string prefix_;
        OptionsMap options_;
        size_t memory_;
        QAM finalam_;

        int simdlen_;
        bool useintrinsics_;
        std::map<std::string, std::string> intrinsicmap_;

        std::set<std::string> cpuflags_;
};


#endif
