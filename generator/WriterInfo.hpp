#ifndef WRITERINFO_H
#define WRITERINFO_H

#include "generator/Options.hpp"
#include "generator/Classes.hpp"

namespace WriterInfo {

        void Init(const OptionsMap & options, const QAM & finalam,
                  const std::string & cpuflags);

        void SetContQ(const QAMSet & topquartets);

        void WriteIncludes(std::ostream & os);
        void AddNamedConstant(const std::string & name, const std::string & val);
        void AddIntConstant(int i);
        void WriteConstants(std::ostream & os);

        void DeclareContwork(std::ostream & os);
        void ZeroContWork(std::ostream & os);
        void FreeContwork(std::ostream & os);

        bool IsContArray(const QAM & am);
        
        size_t MemoryReq(void);

        void WriteAccumulation(std::ostream & os);
        void WriteShellOffsets(std::ostream & os);
        void WriteShellOffsets_Scalar(std::ostream & os);

        std::string NamedConstant(const std::string & cname);
        std::string IntConstant(int i);


}


#endif
