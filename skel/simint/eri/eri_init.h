#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Initializes the ERI functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_eri_init(void);


/*! \brief Finalizes the ERI functionality
 *
 * \warning This is not expected to be called directly from
 *          outside the library
 */
void simint_eri_finalize(void);


#ifdef __cplusplus
}
#endif

