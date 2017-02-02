/*! \file
 *
 * \brief Standard, consistent naming of arrays of integrals
 * \author Benjamin Pritchard (ben@bennyp.org)
 */


#ifndef SIMINT_GUARD_GENERATOR__NAMING_HPP_
#define SIMINT_GUARD_GENERATOR__NAMING_HPP_

#include <string>
#include <sstream>
#include "generator/StringBuilder.hpp"

/*! \brief Generate a general name of an array of integrals
 *
 * \param [in] tam Tagged AM of the integral to be stored
 * \param [in] prefix Optional prefix to use
 * \return A standard name of an array
 */
inline std::string ArrVarName(const QAM & tam, const std::string & prefix = "")
{
    auto am = tam.qam;

    std::string tag = (tam.tag.size() ? tam.tag + "_" : "");

    std::string str = StringBuilder("INT__" , tag, amchar[am[0]], "_", 
                                    amchar[am[1]], "_", amchar[am[2]], "_", amchar[am[3]]);
    
    if(prefix.size())
        return StringBuilder(prefix, "_", str);
    else
        return str;
}


/*! \brief Generate a general name of an array of integrals (split with custom ket string)
 *
 * \param [in] am1 AM of the first center of the bra
 * \param [in] am2 AM of the second center of the bra
 * \param [in] ketstr Custom string for the ket part
 * \param [in] prefix Optional prefix to use
 * \return A standard name of an array
 */
inline std::string ArrVarName(int am1, int am2, const std::string & ketstr, const std::string & prefix = "")
{
    std::string str = StringBuilder("INT__" , amchar[am1], "_", amchar[am2], "_", ketstr);
    
    if(prefix.size())
        return StringBuilder(prefix, "_", str);
    else
        return str;
}

/*! \brief Generate a general name of an array of integrals (doublet with custom ket string)
 *
 * \param [in] am1 AM of the first center of the bra
 * \param [in] am2 AM of the second center of the bra
 * \param [in] ketstr Custom string for the ket part
 * \param [in] prefix Optional prefix to use
 * \return A standard name of an array
 */
inline std::string ArrVarName(Doublet dbra, const std::string & ketstr, const std::string & prefix = "")
{
    std::string tag = (dbra.tag.size() ? dbra.tag + "_" : "");
    std::string str = StringBuilder("INT__" , tag, amchar[dbra.left.am()],
                                    "_", amchar[dbra.right.am()], "_", ketstr);
    
    if(prefix.size())
        return StringBuilder(prefix, "_", str);
    else
        return str;
}


/*! \brief Generate a general name of an array of integrals (split with custom bra string)
 *
 * \param [in] brastr Custom string for the bra part
 * \param [in] am3 AM of the first center of the bra
 * \param [in] am4 AM of the second center of the bra
 * \param [in] prefix Optional prefix to use
 * \return A standard name of an array
 */
inline std::string ArrVarName(const std::string & brastr, int am3, int am4, const std::string & prefix = "")
{
    std::string str = StringBuilder("INT__" , brastr, "_", amchar[am3], "_", amchar[am4]);
    
    if(prefix.size())
        return StringBuilder(prefix, "_", str);
    else
        return str;
}


/*! \brief Generate a general name of an array of integrals (doublet with custom bra string)
 *
 * \param [in] brastr Custom string for the bra part
 * \param [in] am3 AM of the first center of the bra
 * \param [in] am4 AM of the second center of the bra
 * \param [in] prefix Optional prefix to use
 * \return A standard name of an array
 */
inline std::string ArrVarName(const std::string & brastr, Doublet dket, const std::string & prefix = "")
{
    std::string tag = (dket.tag.size() ? dket.tag + "_" : "");
    std::string str = StringBuilder("INT__" , tag, brastr, "_", amchar[dket.left.am()],
                                    "_", amchar[dket.right.am()]);
    
    if(prefix.size())
        return StringBuilder(prefix, "_", str);
    else
        return str;
}



/*! \brief Generate a general name of an array of integrals used in HRR
 *
 * \param [in] am AM of the integral to be stored
 * \return A standard name of an array used in HRR
 */
inline std::string HRRVarName(const QAM & am)
{
    return ArrVarName(am, "HRR");
}



/*! \brief Generate a general name of an array of integrals used in HRR (split with custom ket string)
 *
 * \param [in] am1 AM of the first center of the bra
 * \param [in] am2 AM of the second center of the bra
 * \param [in] ketstr Custom string for the ket part
 * \return A standard name of an array used in HRR
 */
inline std::string HRRVarName(int am1, int am2, const std::string & ketstr)
{
    return ArrVarName(am1, am2, ketstr, "HRR");
}


/*! \brief Generate a general name of an array of integrals used in HRR (split with custom bra string)
 *
 * \param [in] brastr Custom string for the bra part
 * \param [in] am3 AM of the first center of the bra
 * \param [in] am4 AM of the second center of the bra
 * \return A standard name of an array used in HRR
 */
inline std::string HRRVarName(const std::string & brastr, int am3, int am4)
{
    return ArrVarName(brastr, am3, am4, "HRR");
}


/*! \brief Generate a name for an array of integrals used within the primitive loop
 *
 * \param [in] am AM of the integral to be stored
 * \return A standard name of an array used in the primitive loop
 */
inline std::string PrimVarName(const QAM & am)
{
    return ArrVarName(am, "PRIM");
}


/*! \brief Generate a name for a pointer to an array of integrals used within the primitive loop
 *
 * \param [in] am AM of the integral to be stored
 * \return A standard name of a pointer to an array used in the primitive loop
 */
inline std::string PrimPtrName(const QAM & am)
{
    return ArrVarName(am, "PRIM_PTR");
}


#endif
