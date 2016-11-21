#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Initializes the OSTEI functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_ostei_init(void);


/*! \brief Initializes the OSTEI 1st derivative functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_ostei_deriv1_init(void);


/*! \brief Finalizes the OSTEI functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_ostei_finalize(void);


/*! \brief Finalizes the OSTEI 1st derivative functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_ostei_deriv1_finalize(void);


#ifdef __cplusplus
}
#endif

