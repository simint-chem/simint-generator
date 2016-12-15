/*! \file
 *
 * \brief Options for code generation
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__OPTIONS_HPP_
#define SIMINT_GUARD_GENERATOR__OPTIONS_HPP_

#include <map>

/*! \brief Collection of options available for generation
 */
enum class Option
{
    StackMem,     //!< Limit of the amount of memory to consider when allocating on the stack
    ExternalVRR,  //!< Write external VRR at this L value and above
    GeneralVRR,   //!< Write general VRR at this L value and above
    ExternalHRR,  //!< Write external HRR at this L value and above
    GeneralHRR,   //!< Write general HRR at this L value and above
    Scalar,       //!< Generate Scalar code
};


/*! \brief A map of options to values (integers)
 */
typedef std::map<Option, int> OptionMap;


/*! \brief free function to create a map with default options
 */
inline OptionMap DefaultOptions(void)
{
    return OptionMap{
                      {Option::StackMem, 0},
                      {Option::ExternalVRR, 0},
                      {Option::GeneralVRR, 0},
                      {Option::ExternalHRR, 0},
                      {Option::GeneralHRR, 0},
                      {Option::Scalar, 0},
                    };
}



#endif
