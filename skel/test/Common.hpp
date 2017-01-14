#pragma once

#include <array>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <atomic>

#include "test/Timer.h" 
#include "simint/shell/shell.h"

#define RND_ZERO 1.0e-15

#define NCART(am) ((am>=0)?((((am)+2)*((am)+1))>>1):0)


//! Contributions to timings
struct TimeContrib
{
    std::atomic<TimerType> ticks_shell_pair;
    std::atomic<TimerType> ticks_integrals;
    std::atomic<TimerType> fullticks_integrals;
    std::atomic<TimerType> time_shell_pair;
    std::atomic<TimerType> time_integrals;
    std::atomic<TimerType> fulltime_integrals;

    TimeContrib(const TimeContrib & rhs) noexcept
      : ticks_shell_pair(rhs.ticks_shell_pair.load()),
        ticks_integrals(rhs.ticks_integrals.load()),
        fullticks_integrals(rhs.fullticks_integrals.load()),
        time_shell_pair(rhs.time_shell_pair.load()),
        time_integrals(rhs.time_integrals.load()),
        fulltime_integrals(rhs.fulltime_integrals.load())
    { }

    TimeContrib(TimeContrib && rhs) noexcept
        : TimeContrib(static_cast<const TimeContrib &>(rhs))
    { }

    TimeContrib & operator=(const TimeContrib &) noexcept = delete;
    TimeContrib & operator=(TimeContrib &&) noexcept = delete;

    TimeContrib(void) noexcept
    {
        ticks_shell_pair = ticks_integrals = fullticks_integrals = 0;
        time_shell_pair = time_integrals = fulltime_integrals = 0;
    }

    TimeContrib & operator+=(const TimeContrib & rhs) noexcept
    {
        ticks_shell_pair += rhs.ticks_shell_pair;
        ticks_integrals += rhs.ticks_integrals;
        fullticks_integrals += rhs.fullticks_integrals;
        time_shell_pair += rhs.time_shell_pair;
        time_integrals += rhs.time_integrals;
        fulltime_integrals += rhs.fulltime_integrals;
        return *this;
    }
};


/*! \brief Checks to see if a quartet is a unique permutation
 *
 * ie, checks that i >= j, k >= l, and (i+j) >= (k+l)
 */
bool UniqueQuartet(int i, int j, int k, int l);




//! Just a vector of gaussian shells
typedef std::vector<simint_shell> GaussianVec;

//! Maps a vector of gaussians to their am
typedef std::map<int, GaussianVec> ShellMap;


/*! \brief Check if the integer exponents of a gaussian are valid
 *
 * ie, ijk must all be >= 0
 */
bool ValidGaussian(const std::array<int, 3> & g);

/*! \brief Iterate over cartesian components of a gaussian
 *
 * The elements of the array represent the exponents on x, y, and z.
 * This increments \p g to be the next gaussian. If the next
 * gaussian is valid, the function returns true. Otherwise,
 * it returns false.
 */
bool IterateGaussian(std::array<int, 3> & g);


/*! \brief Deep copies a vector of gaussians
 */
GaussianVec CopyGaussianVec(const GaussianVec & v);

/*! \brief Frees memory associated with a vector of gaussians
 */
void FreeGaussianVec(GaussianVec & agv);

/*! \brief Deep copies a map of gaussian vectors
 */
ShellMap CopyShellMap(const ShellMap & m);

/*! \brief Frees memory associated with a map of gaussian vectors
 */
void FreeShellMap(ShellMap & m);


/*! \brief Reads a basis set from a file into a ShellMap
 */
ShellMap ReadBasis(const std::string & file);


/*! \brief Find some maximum parameters from a basis set (shell map)
 *
 * Returns an array containing
 *   0.) maximum angular momentum
 *   1.) maximum storage needed for cartesians of a shell of a given angular momentum
 */
std::pair<int, int> FindMaxParams(const ShellMap & m);




/*! \brief Chop very small values to zero
 *
 * \param [in] calc Array of doubles to chop
 * \param [in] ncalc Length of the array
 */
void Chop(double * const restrict calc, int ncalc);

/*! \brief Calculate the error relative to a reference
 *
 * \param [in] calc Calculated values
 * \param [in] ref Reference values
 * \param [in] ncalc Number of values in \p calc and \p ref
 * \return A pair representing the maximum absolute and relative error found
 */
std::pair<double, double> CalcError(double const * const restrict calc,
                                    double const * const restrict ref,
                                    int ncalc);


/*! \brief Print the header for timings */
void PrintTimingHeader(void);

/*! \brief Print a line with the timings for a particular AM quartet */
void PrintAMTimingInfo(int i, int j, int k, int l, size_t nshell1234, size_t nprim1234, const TimeContrib & info);

