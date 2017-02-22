#pragma once

#include <stddef.h>
#include "simint/ostei/ostei_config.h"


#define SIMINT_SCREEN_NONE         0
#define SIMINT_SCREEN_SCHWARZ      1
#define SIMINT_SCREEN_FASTSCHWARZ  2


#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Information about a gaussian shell */
struct simint_shell
{
    int am;          //!< Angular momentum (0 = s, etc)
    int nprim;       //!< Number of primitives in this shell

    double x;        //!< X coordinate (in bohr)
    double y;        //!< Y coordinate (in bohr)
    double z;        //!< Z coordinate (in bohr)

    double * alpha;  //!< Exponents of the gaussian functions (for each primitive)
    double * coef;   //!< Contraction coefficients (for each primitive)

    size_t memsize;     //!< Total memory for storing various data in this structure (in bytes)
    void * ptr;      //!< Pointer to all the allocated memory within this structure (in bytes)
};


/*! \brief A structure holding information about multiple shell pair
 *
 * A shell pair is a combination of two gaussian shells occuring in the bra
 * or ket of an integral. Several factors used in calculation of integrals
 * can be computed solely from a pair of shells (that is, it doesn't require
 * information from all four centers). This include some factors from the
 * gaussian product theorem.
 *
 * The only real requirement is that all shells in the first position have the
 * same angular momentum, and all shells in the second position have the same
 * angular momentum.
 */
struct simint_multi_shellpair
{
    int am1;            //!< Angular momentum of the first position
    int am2;            //!< Angular momentum of the second position
    int nprim;          //!< Total number of primitive combinations stored (not including padding)

    int nshell12;       //!< Total number of shell pair (nshell1 * nshell2)
    int nshell12_clip;  //!< Total number of shell pair to actual calculate (should be <= nshell12)
    int * nprim12;      //!< Number of primitive combinations for each shell pair, not including padding (length nshell12)

    double * AB_x;      //!< X distance between the centers for a shell (Ax - Bx) (length nshell12).
    double * AB_y;      //!< Y distance between the centers for a shell (Ay - By) (length nshell12).
    double * AB_z;      //!< Z distance between the centers for a shell (Az - Bz) (length nshell12).

    double * x;         //!< X coordinate of each primitive pair (from GPT)
    double * y;         //!< y coordinate of each primitive pair (from GPT)
    double * z;         //!< z coordinate of each primitive pair (from GPT)
    double * PA_x;      //!< Px - Ax
    double * PA_y;      //!< Py - Ay
    double * PA_z;      //!< Pz - Az
    double * PB_x;      //!< Px - Bx
    double * PB_y;      //!< Py - By
    double * PB_z;      //!< Pz - Bz

    double * alpha;     //!< New coefficients (from GPT)

    #if SIMINT_OSTEI_MAXDER > 0
    double * alpha2;    //!< 2*exponent on the first center
    double * beta2;     //!< 2*exponent on the second center
    #endif

    double * prefac;    //!< Prefactors for each primitive pair, including coefficients and other factors
    double * screen;    //!< Screening information (value of g_{abab} for all primitive shell pair)
    double screen_max;  //!< Maximum value in the screen array


    size_t memsize;     //!< Total memory for storing various data in this structure (in bytes)
    void * ptr;         //!< Pointer to all the allocated memory within this structure (length memsize)
};


/*! \brief See if two shells are equivalent */
static inline
int compare_shell(struct simint_shell const * A,
                  struct simint_shell const * B)
{
    return A->nprim == B->nprim && A->ptr == B->ptr;
}


/*! \brief Initialize a shell structure
 *
 * This sets certain values to zero so that they can be used with
 * simint_allocate_shell, etc.
 *
 * \param [inout] G The structure to initialize
 */
void simint_initialize_shell(struct simint_shell * G);


/*! \brief Allocate memory in a shell
 *
 * Allocate enough memory for the exponents and coefficients, and sets up
 * the \p alpha and \coef pointers. This also sets the \p ptr and \p memsize
 * members of \p G
 */
void simint_allocate_shell(int nprim, struct simint_shell * G);


/*! \brief Frees memory associated with a shell structure
 *
 * The \p ptr member will be set to NULL and the memsize will be set to zero.
 */
void simint_free_shell(struct simint_shell * G);


/*! \brief Copies a shell structure
 *
 * Memory will be allocated for a new shell and all the data copied from \p src
 * to \p dest
 */
void simint_copy_shell(struct simint_shell const * src,
                       struct simint_shell * dest);


/*! \brief Normalize the coefficients of shells
 *
 * The normalization in this function is what is expected in the simint library
 *
 * \param [in] n The number of shells to normalize
 * \param [inout] G Pointer to the shells to normalize
 */
void simint_normalize_shells(int n, struct simint_shell * G);


/*! \brief Create a simint shell
 *
 * Just for convenience. You may also allocate and then
 * manually fill in the data members yourself.
 *
 * Data will be copied from the \p alpha and \p coef
 * pointers.
 *
 * The shell must be freed later (via simint_free_shell)
 *
 * \param [in] nprim Number of primitives in the shell
 * \param [in] am Angular momentum of the shell
 * \param [in] x The x-coordinate of the shell
 * \param [in] y The y-coordinate of the shell
 * \param [in] z The z-coordinate of the shell
 * \param [in] alpha The alpha parameters of the primitives
 * \param [in] coef The coefficients of the primitives
 * \param [inout] G The gaussian shell to use
 */
void simint_create_shell(int nprim, int am, double x, double y, double z,
                         double const * alpha,
                         double const * coef,
                         struct simint_shell * G);


/*! \brief Initialize a shell pair structure
 *
 * This sets certain values to zero so that they can be used with
 * simint_allocate_multi_shellpair, etc.
 *
 * \param [inout] P The structure to initialize
 */
void simint_initialize_multi_shellpair(struct simint_multi_shellpair * P);


/*! \brief Allocate space in a multi shellpair structure
 *
 * Only sets up the pointers in \P and fills in the \p memsize and \p ptr members.
 *
 * \param [in] na Number of shells in the first position
 * \param [in] A Shells that will be stored in the first position of this multi_shellpair
 * \param [in] nb Number of shells in the second position
 * \param [in] B Shells that will be stored in the second position of this multi_shellpair
 * \param [inout] P The structure in which to allocate the memory
 */
void simint_allocate_multi_shellpair(int na, struct simint_shell const * A,
                                     int nb, struct simint_shell const * B,
                                     struct simint_multi_shellpair * P,
                                     int screen_method);


/*! \brief Allocate space in a multi shellpair structure
 *
 * Only sets up the pointers in \P and fills in the \p memsize and \p ptr members.
 *
 * \param [in] npair Number of shell pairs in the array
 * \param [in] AB Pairs of shells to place in the shell pair
 * \param [inout] P The structure in which to allocate the memory
 */
void simint_allocate_multi_shellpair2(int npair, struct simint_shell const * AB,
                                      struct simint_multi_shellpair * P,
                                      int screen_method);


/*! \brief Frees memory associated with a multi shellpair structure
 *
 * The \p ptr member will be set to NULL and the memsize will be set to zero.
 */
void simint_free_multi_shellpair(struct simint_multi_shellpair * P);


/*! \brief Computes and fills in values for a multi shellpair
 *
 * This calculates the values for all the members of a simint_multishellpair structure by
 * looping over all combinations of A and B and forming the shell pair information.
 *
 * \param [in] na Number of shells in the first position
 * \param [in] A Shells in the first position of this multi_shellpair
 * \param [in] nb Number of shells in the second position
 * \param [in] B Shells in the second position of this multi_shellpair
 * \param [inout] P The structure that will hold the shell pair data
 *
 * \warning \p P must already be allocated (via simint_allocate_multi_shellpair)
 */
void simint_fill_multi_shellpair(int na, struct simint_shell const * A,
                                 int nb, struct simint_shell const * B,
                                 struct simint_multi_shellpair * P,
                                 int screen_method);


/*! \brief Computes and fills in values for a multi shellpair
 *
 * This calculates the values for all the members of a simint_multishellpair structure by
 * looping through the array of pairs.
 *
 * The array is expected to be 2*npair in length, with the pairs being adjacent. Ie,
 *
 * AB = [A1 B1 A1 B2 A2 B1 A2 B2]  with npair = 4
 *
 * \param [in] npair Number of shell pairs in the array
 * \param [in] AB Pairs of shells to place in the shell pair
 * \param [inout] P The structure that will hold the shell pair data
 * \param [in] screen_method Screening method for primitives
 *
 * \warning \p P must already be allocated (via simint_allocate_multi_shellpair)
 */
void simint_fill_multi_shellpair2(int npair, struct simint_shell const * AB,
                                  struct simint_multi_shellpair * P,
                                  int screen_method);


/*! \brief Allocates and fills a multi shellpair structure
 *
 * For convenience. Creates a new simint_multi_shellpair structure,
 * allocates it, and then calculates all the data.
 *
 * \param [in] na Number of shells in the first position
 * \param [in] A Shells in the first position of this multi_shellpair
 * \param [in] nb Number of shells in the second position
 * \param [in] B Shells in the second position of this multi_shellpair
 * \param [inout] P The structure that will hold the shell pair data
 * \param [in] screen_method Screening method for primitives
 */
void
simint_create_multi_shellpair(int na, struct simint_shell const * A,
                              int nb, struct simint_shell const * B,
                              struct simint_multi_shellpair * P,
                              int screen_method);


/*! \brief Allocates and fills a multi shellpair structure
 *
 * For convenience. Creates a new simint_multi_shellpair structure,
 * allocates it, and then calculates all the data.
 *
 * \param [in] npair Number of shell pairs in the array
 * \param [in] AB Pairs of shells to place in the shell pair
 * \param [inout] P The structure that will hold the shell pair data
 * \param [in] screen_method Screening method for primitives
 */
void
simint_create_multi_shellpair2(int npair,
                               struct simint_shell const * AB,
                               struct simint_multi_shellpair * P,
                               int screen_method);


/*! \brief Combine existing multi shellpair structures into a new one
 *
 * Existing information in \Pout will be erased
 *
 * \param [in] nmpair Number of multi shellpair pointers in \p Pin
 * \param [in] Pin Array of pointers to multi shellpair to combine
 * \param [inout] Pout The structure that will hold the new shell pair data
 * \param [in] screen_method Screening method for primitives
 */
void simint_cat_multi_shellpair(int nmpair,
                                struct simint_multi_shellpair const ** Pin,
                                struct simint_multi_shellpair * Pout,
                                int screen_method);


/*! \brief Remove all insignificant primitive pairs
 *
 *
 * \param [in] npair Number of shell pairs in the array
 * \param [in] AB Pairs of shells to place in the shell pair
 * \param [inout] P The structure that will hold the shell pair data
 */
/*
void
simint_prune_multi_shellpair(struct simint_multi_shellpair const * P,
                             struct simint_multi_shellpair * out,
                             double screen_max, double screen_tol);
*/

#ifdef __cplusplus
}
#endif

