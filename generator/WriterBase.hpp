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

        void DeclareContwork(std::ostream & os) const;
        void ZeroContWork(std::ostream & os, const std::string & nshell) const;
        void FreeContwork(std::ostream & os) const;

        int GetOption(int option) const;
        bool HasCPUFlag(const std::string & flag) const;

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
        std::set<std::string> cpuflags_;
};


#endif
