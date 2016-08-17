#pragma once

#include "simint/simint_config.h" // for USE_ET define

#include <stddef.h>
#include <stdint.h>

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
    int nprim;          //!< Total number of primitive combinations stored
    int nprim_length;   //!< Acutal length of primitive arrays (alpha, etc), including padding

    int nshell12;       //!< Total number of shell pair (nshell1 * nshell2)
    int * nprim12;      //!< Number of primitive combinations for each shell pair (length nshell12) 

    int nbatch;         //!< Number of batches of shell pairs
    int * nbatchprim;   //!< Number of primitives in each batch, including padding (length nbatch)

    double * AB_x;      //!< X distance between the centers for a shell (Ax - Bx) (length nshell12).
    double * AB_y;      //!< Y distance between the centers for a shell (Ay - By) (length nshell12).
    double * AB_z;      //!< Z distance between the centers for a shell (Az - Bz) (length nshell12).

    double * x;         //!< X coordinate of each primitive pair (from GPT) (length nprim)
    double * y;         //!< y coordinate of each primitive pair (from GPT) (length nprim)
    double * z;         //!< z coordinate of each primitive pair (from GPT) (length nprim)
    double * PA_x;      //!< Px - Ax (length nprim)
    double * PA_y;      //!< Py - Ay (length nprim)
    double * PA_z;      //!< Pz - Az (length nprim)
    double * PB_x;      //!< Px - Bx (length nprim)
    double * PB_y;      //!< Py - By (length nprim)
    double * PB_z;      //!< Pz - Bz (length nprim)

    #ifdef SIMINT_ERI_USE_ET
    double * bAB_x;     //!< Balpha * (Ax - Bx) (length nprim)
    double * bAB_y;     //!< Balpha * (Ay - By) (length nprim)
    double * bAB_z;     //!< Balpha * (Az - Bz) (length nprim)
    #endif

    double * alpha;     //!< New coefficients (from GPT) (length nprim)
    double * prefac;    //!< Prefactors for each primitive pair, including coefficients and other factors


    size_t memsize;     //!< Total memory for storing various data in this structure (in bytes)
    void * ptr;         //!< Pointer to all the allocated memory within this structure (length memsize)
};


/*! \brief Allocate memory in a shell
 *
 * Allocate enough memory for the exponents and coefficients, and sets up
 * the \p alpha and \coef pointers. This also sets the \p ptr and \p memsize
 * members of \p G
 */
void simint_allocate_shell(int nprim, struct simint_shell * const restrict G);


/*! \brief Frees memory associated with a shell structure
 *
 * The \p ptr member will be set to NULL and the memsize will be set to zero.
 */
void simint_free_shell(struct simint_shell G);


/*! \brief Copies a shell structure
 *
 * Memory will be allocated for a new shell and all the data copied from \p G
 */
struct simint_shell simint_copy_shell(const struct simint_shell G);


/*! \brief Normalize the coefficients of shells
 *
 * The normalization in this function is what is expected in the simint library
 *
 * \param [in] n The number of shells to normalize
 * \param [inout] G Pointer to the shells to normalize
 */
void simint_normalize_shells(int n, struct simint_shell * const restrict G);


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
void simint_allocate_multi_shellpair(int na, struct simint_shell const * const restrict A,
                                     int nb, struct simint_shell const * const restrict B,
                                     struct simint_multi_shellpair * const restrict P);


/*! \brief Allocate space in a multi shellpair structure
 *
 * Only sets up the pointers in \P and fills in the \p memsize and \p ptr members.
 *
 * \param [in] npair Number of shell pairs in the array
 * \param [in] AB Pairs of shells to place in the shell pair
 * \param [inout] P The structure in which to allocate the memory
 */
void simint_allocate_multi_shellpair2(int npair, struct simint_shell const * const restrict AB,
                                      struct simint_multi_shellpair * const restrict P);


/*! \brief Frees memory associated with a multi shellpair structure
 *
 * The \p ptr member will be set to NULL and the memsize will be set to zero.
 */
void simint_free_multi_shellpair(struct simint_multi_shellpair P);


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
void simint_fill_multi_shellpair(int na, struct simint_shell const * const restrict A,
                                 int nb, struct simint_shell const * const restrict B,
                                 struct simint_multi_shellpair * const restrict P);



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
 *
 * \warning \p P must already be allocated (via simint_allocate_multi_shellpair)
 */
void simint_fill_multi_shellpair2(int npair, struct simint_shell const * const restrict AB,
                                  struct simint_multi_shellpair * const restrict P);




/*! \brief Allocates and fills a multi shellpair structure
 *
 * For convenience. Creates a new simint_multi_shellpair structure,
 * allocates it, and then calculates all the data.
 *
 * \param [in] na Number of shells in the first position
 * \param [in] A Shells in the first position of this multi_shellpair
 * \param [in] nb Number of shells in the second position
 * \param [in] B Shells in the second position of this multi_shellpair
 */
struct simint_multi_shellpair
simint_create_multi_shellpair(int na, struct simint_shell const * const restrict A,
                              int nb, struct simint_shell const * const restrict B);

/*! \brief Allocates and fills a multi shellpair structure
 *
 * For convenience. Creates a new simint_multi_shellpair structure,
 * allocates it, and then calculates all the data.
 *
 * \param [in] npair Number of shell pairs in the array
 * \param [in] AB Pairs of shells to place in the shell pair
 */
struct simint_multi_shellpair
simint_create_multi_shellpair2(int npair, struct simint_shell const * const restrict AB);


#ifdef __cplusplus
}
#endif

