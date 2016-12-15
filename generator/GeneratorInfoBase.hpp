/*! \file
 *
 * \brief Holds information about the requested generation (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__GENERATORINFOBASE_HPP_
#define SIMINT_GUARD_GENERATOR__GENERATORINFOBASE_HPP_

#include "generator/Types.hpp"
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
     * \param [in] deriv Derivative to generate code for
     * \param [in] options Options for code generation
     */
    GeneratorInfoBase(QAM finalam, int deriv, const OptionMap & options)
        : finalam_(finalam), deriv_(deriv), options_(options) { }


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

    /*! \brief Generate scalar code?
     */   
    bool Scalar(void) const
    {
        return GetOption(Option::Scalar);
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

    //! Code generation options
    OptionMap options_;

    //! Generate scalar code?
    bool scalar_;
};



#endif
