#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : main.c                                         */
/*                                                                           */
/* Created:       2011/11/29 (JLe)                                           */
/* Last modified: 2012/12/28 (JLe)                                           */
/* Version:       2.1.12                                                     */
/*                                                                           */
/* Description: Finalizes MPI before exiting program                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "FinalizeMPI:"

/*****************************************************************************/

void FinalizeMPI()
{
#ifdef MPI

  /* Synchronise */

  MPI_Barrier(MPI_COMM_WORLD);

  /* Finalize */

  if (MPI_Finalize() != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");

#endif
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
