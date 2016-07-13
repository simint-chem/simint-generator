#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Initializes the simint library
 *
 * This function performs any required initialization of variables, data
 * structures, etc, within libint.
 *
 * The library must be initialized for a given process (not necessarily a given
 * thread).
 *
 * It is safe to initialize the library multiple times. After the first call,
 * this function has no effect.
 * 
 * \warning Is not thread safe (ie, should not be called at the same time2
 *          from multiple threads).
 */
void simint_init(void);


/*! \brief Cleans up any resources used by the simint library
 *
 * This function cleans up any data structures within the simint library.
 * This must be run when you are no longer using any simint functionality.
 * Afterwards, simint may crash if you attempt to use some functions within
 * simint.
 *
 * It is safe to finalize the library multiple times. After the first call, this
 * function has no effect.
 *
 * \warning Is not thread safe (ie, should not be called at the same time2
 *          from multiple threads).
 */
void simint_finalize(void);


#ifdef __cplusplus
}
#endif
