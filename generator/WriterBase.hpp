#ifndef WRITERBASE_H
#define WRITERBASE_H

#include "generator/Options.hpp"
#include "generator/Classes.hpp"

class WriterBase
{
    public:
        WriterBase(const OptionsMap & options, const QAM & finalam);

        void SetContQ(const QuartetSet & topquartets);

        void DeclarePrimArrays(std::ostream & os) const;
        void DeclareContwork(std::ostream & os) const;
        void FreeContwork(std::ostream & os) const;

        int GetOption(int option) const;

        bool IsContArray(const QAM & am) const;
        
        QAM FinalAM(void) const;

        bool IsFinalAM(const QAM & am) const;

        size_t MemoryReq(void) const;

        int L(void) const;

        void PermuteResult(std::ostream & os, const std::string & src) const;

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
        OptionsMap options_;
        size_t memory_;
        QAM finalam_;
};


#endif
