/*! \file
 *
 * \brief Holds information about the requested generation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__GENERATORINFOBASE_HPP_
#define SIMINT_GUARD_GENERATOR__GENERATORINFOBASE_HPP_

#include <memory>
#include "generator/Types.hpp"
#include "generator/VectorInfo.hpp"
#include "generator/Ncart.hpp"
#include "generator/Options.hpp"



/*! \brief Information about the current generation run
 *
 * This class holds information about the intrinsics, final requested
 * quartet, CPU capabilities, etc.
 */
class GeneratorInfoBase
{
public:    

    /*! \brief Constructor
     * 
     * Vector information is is determined based on the CPU flags passed in through \p cpuflagsstr
     * 
     * \param [in] finalam The final target AM quartet desired
     * \param [in] compiler The compiler to target for the generated code
     * \param [in] cpuflagsstr A comma-separated list of CPU flags to consider when generating the code
     * \param [in] options Options for code generation
     */
    GeneratorInfoBase(QAM finalam,
                      int deriv,
                      Compiler compiler,
                      const std::string & cpuflagsstr,
                      const OptionMap & options);


    virtual ~GeneratorInfoBase() = default;

    /*! \brief The final target AM for this generation run
     */   
    QAM FinalAM(void) const
    {
        return finalam_;
    }

    /*! \brief The derivative we are calculating
     */
    int Deriv(void) const
    {
        return deriv_;
    }

    /*! \brief Retrieve an option
     * 
     * Will throw an exception if the option has not been set yet
     *
     * \param [in] opt Option to retrieve
     * \return Value for option \p opt
     */   
    int GetOption(Option opt) const
    {
        return options_.at(opt);
    }

    /*! \brief Total angular momentum quantum number (sum of AM for all four centers)
     */   
    int L(void) const
    {
        return finalam_[0] + finalam_[1] + finalam_[2] + finalam_[3];
    }

    /*! \brief Compare an AM quartet against the final target AM
     * 
     * \param [in] am  The am quartet to check
     * \return True if \p am is the same as the target AM quartet
     */   
    bool IsFinalAM(const QAM & am) const
    {
        return am == finalam_;
    }

    /*! \brief Check for the existance of a CPU flags
     *
     * \param [in] flag Flag to check. Must be lowercase
     * \return True if that CPU flag has been set.
     */   
    bool HasCPUFlag(const std::string & flag) const
    {
        return cpuflags_.count(flag);
    }

    /*! \brief Check if FMA is available (based on CPU flags)
     * 
     * \return True if the fma CPU flags has been set
     */   
    bool HasFMA(void) const
    {
        return HasCPUFlag("fma");
    }

    /*! \brief Get the information regarding the vector types
     */   
    const VectorInfo & GetVectorInfo(void) const
    {
        return *(vector_);
    }

    /*! \brief Generate scalar code?
     */   
    bool Scalar(void) const
    {
        return scalar_;
    }

    /*! \brief Generate vectorized code?
     */   
    bool Vectorized(void) const
    {
        return !Scalar();
    }


private:
    //! The requested AM quartet
    QAM finalam_;

    //! The requested derivative
    int deriv_;

    //! The compiler the generated source is to be compiled with
    Compiler compiler_;

    //! CPU flags to consider for intrinsics and CPU capabilities
    StringSet cpuflags_;

    //! Code generation options
    OptionMap options_;

    //! Generate scalar code?
    bool scalar_;

    //! Get the vector intrinsics, etc, for this generation run
    std::unique_ptr<VectorInfo> vector_;

    /*! \brief Split a string of cpuflags into a set
     */
    static StringSet ConvertCPUFlags(std::string cpuflagsstr);
};



#endif
