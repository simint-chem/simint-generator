#ifndef WRITERINFO_H
#define WRITERINFO_H

#include "generator/Options.hpp"
#include "generator/Classes.hpp"

namespace WriterInfo {

        void Init(const OptionsMap & options, const QAM & finalam,
                  const std::string & cpuinfofile);

        void ReadCPUFlags(const std::string & file);

        void SetContQ(const QAMSet & topquartets);

        void WriteIncludes(std::ostream & os);
        void AddNamedConstant(const std::string & name, const std::string & val);
        void AddIntConstant(int i);
        void WriteConstants(std::ostream & os);

        void DeclareContwork(std::ostream & os);
        void ZeroContWork(std::ostream & os);
        void FreeContwork(std::ostream & os);

        int GetOption(int option);
        bool HasCPUFlag(const std::string & flag);
        bool HasFMA(void);
        bool Intrinsics(void);
        bool Scalar(void);

        bool IsContArray(const QAM & am);
        
        QAM FinalAM(void);

        bool IsFinalAM(const QAM & am);

        bool HasVRR(void);
        bool HasET(void);
        bool HasHRR(void);
        bool HasBraHRR(void);
        bool HasKetHRR(void);

        size_t MemoryReq(void);

        int L(void);

        int SimdLen(void);
        int ByteAlign(void);

        void WriteAccumulation(std::ostream & os, QAM qam, int ncart);

        std::string NamedConstant(const std::string & cname);
        std::string IntConstant(int i);

        std::string DoubleType(void);
        std::string ConstDoubleType(void);
        std::string DoubleSet(const std::vector<std::string> & dbls);
        std::string DoubleSet1(const std::string & dbl);
        std::string DoubleLoad(const std::string & ptr, const std::string & idx);
        std::string DoubleStore(const std::string & var, const std::string & ptr, const std::string & idx);
        std::string NewDoubleSet1(const std::string & var, const std::string & val);
        std::string NewDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx);
        std::string NewConstDoubleSet1(const std::string & var, const std::string & val);
        std::string NewConstDoubleLoad(const std::string & var, const std::string & ptr, const std::string & idx);

        std::string UnionType(void);
        std::string ConstUnionType(void);

        std::string FMAdd(const std::string & a, const std::string & b, const std::string & c);
        std::string FMSub(const std::string & a, const std::string & b, const std::string & c);
        std::string Sqrt(const std::string & val);
        std::string RSqrt(const std::string & val);
        std::string Power(const std::string & base, const std::string & exp);
        std::string Exp(const std::string & exp);

        std::string ArrVarName(const QAM & am, const std::string & prefix = "");
        std::string ArrVarName(int am1, int am2, const std::string & ketstr, const std::string & prefix = "");
        std::string ArrVarName(const std::string & brastr, int am3, int am4, const std::string & prefix = "");

        std::string HRRVarName(const QAM & am);
        std::string HRRVarName(int am1, int am2, const std::string & ketstr);
        std::string HRRVarName(const std::string & brastr, int am3, int am4);

        std::string PrimVarName(const QAM & am);
        std::string PrimPtrName(const QAM & am);
}


#endif
