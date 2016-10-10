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


/*! \brief Finalizes the OSTEI functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_ostei_finalize(void);


#ifdef __cplusplus
}
#endif

