#ifndef WRITERBASE_H
#define WRITERBASE_H

#include "generator/Options.hpp"
#include "generator/Classes.hpp"

class WriterBase
{
    public:
        WriterBase(const OptionsMap & options, const QAMList & finalam);

        void SetContQ(const QuartetSet & topquartets, const DoubletSetMap & topkets);

        void DeclareContwork(std::ostream & os) const;
        void FreeContwork(std::ostream & os) const;

        int GetOption(int option) const;

        bool IsContArray(const QAMList & am) const;
        
        QAMList FinalAM(void) const;

        size_t MemoryReq(void) const;

        int L(void) const;

        static std::string ArrVarName(const QAMList & am);

        static std::string AuxName(int i);



    private:
        QAMListSet contq_;  // set of contracted integral AM
        OptionsMap options_;
        size_t memory_;
        QAMList finalam_;

        std::string eri_prefix_;
        std::string vrr_prefix_;
        std::string et_prefix_;
        std::string hrr_prefix_;
};


#endif
