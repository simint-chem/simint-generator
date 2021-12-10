/*! \file
 *
 * \brief Some helpers for parsing the command line (header)
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__COMMANDLINE_HPP_
#define SIMINT_GUARD_GENERATOR__COMMANDLINE_HPP_

#include <vector>
#include <string>
#include "generator/Options.hpp"


/*! \brief Get the next argument on the command line
 *
 * The parameter \p i will be incremented if possible
 *
 * \throw std::runtime_error if there is no more arguments
 * \param [inout] i Index to the argument on the command line
 * \param [in] argc Argument count from the command line
 * \param [in] argv Arguments from the command line
 * \return The next argument on the command line (as a string)
 */
std::string GetNextArg(int & i, int argc, char ** argv);


/*! \brief Get the next argument on the command line (and convert to an integer)
 *
 * The parameter \p i will be incremented if possible
 *
 * \throw std::runtime_error if there is no more arguments
 * \param [inout] i Index to the argument on the command line
 * \param [in] argc Argument count from the command line
 * \param [in] argv Arguments from the command line
 * \return The next argument on the command line (as an int)
 */
int GetIArg(int & i, int argc, char ** argv);


/*! \brief Get the next argument on the command line
 *
 * Version for command lines that have been converted to a vector of strings).
 *
 * The parameter \p i will be incremented if possible
 *
 * \throw std::runtime_error if there is no more arguments
 * \param [inout] i Index to the argument on the command line
 * \param [in] opt Options to parse 
 * \return The next argument on the command line (as an int)
 */
std::string GetNextArg(size_t & i, const std::vector<std::string> & opt);

/*! \brief Get the next argument on the command line (and convert to an integer)
 *
 * Version for command lines that have been converted to a vector of strings).
 *
 * The parameter \p i will be incremented if possible
 *
 * \throw std::runtime_error if there is no more arguments
 * \param [inout] i Index to the argument on the command line
 * \param [in] opt Vector of strings containing the command line arguments
 * \return The next argument on the command line (as an int)
 */
int GetIArg(size_t & i, const std::vector<std::string> & opt);


/*! \brief Parse common options on the command line
 *
 * Results will be stored in the \p options parameter.
 *
 * \param [inout] options Options object to set the options in
 * \param [in] argc Argument count from the command line
 * \param [in] argv Arguments from the command line
 * \return Any leftover options that haven't been parsed
 */
std::vector<std::string> ParseCommonOptions(OptionMap & options, int argc, char ** argv);


//! An assert to check if something has been set
#define CMDLINE_ASSERT(val, desc) if( !(val) ) { std::cout << "\n" << (desc) << "\n\n"; return 1; }


#endif
