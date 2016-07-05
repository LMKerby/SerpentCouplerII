/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collectparalleldata.c                          */
/*                                                                           */
/* Created:       2010/11/23 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Collects results from parallel MPI tasks                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectParallelData:"

/*****************************************************************************/

void CollectParallelData()
{

#ifdef MPI

  double *buff;
  long sz;

  /* Check if access to private arrays is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "Private arrays not ready for access");

  /* Check number of tasks */

  if (mpitasks == 1)
    return;

  fprintf(out, "Waiting for results from other MPI tasks...\n");

  /* Start timers */

  StartTimer(TIMER_MPI_OVERHEAD);
  StartTimer(TIMER_MPI_OVERHEAD_TOTAL);
  
  /***************************************************************************/
  
  /***** Data in RES1 array **************************************************/
  
  /* Check MPI reproducibility option */

  if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
    {
      /* Get size of results data block */

      sz = (long)RDB[DATA_ALLOC_RES1_SIZE];
      
      /* Allocate memory for results */
      
      if (mpiid == 0)
	buff = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
      else
	buff = NULL;
      
      /* Synchronise */
      
      MPI_Barrier(MPI_COMM_WORLD);
      
      /* Reduce data */
      
      MPITransfer(RES1, buff, sz, 0, MPI_METH_RED);
      
      /* Move data to original block */
      
      if (mpiid == 0)
	{  
	  /* Copy data */
	  
	  memcpy(RES1, buff, sz*sizeof(double));
	  
	  /* Free buffer */
	  
	  Mem(MEM_FREE, buff);     
	}
      
      /* Synchronise */
      
      MPI_Barrier(MPI_COMM_WORLD);
      
      /* Broadcast data to other tasks */
      
      MPITransfer(RES1, NULL, sz, 0, MPI_METH_BC);
    }

  /***************************************************************************/

  /***** Data in RES2 array **************************************************/

  /* Reduce private results */

  ReducePrivateRes();

  /* Get size of results data block */

  sz = (long)RDB[DATA_ALLOC_RES2_SIZE];

  /* Allocate memory for results */

  if (mpiid == 0)
    buff = (double *)Mem(MEM_ALLOC, sz, sizeof(double));
  else
    buff = NULL;

  /* Synchronise */

  MPI_Barrier(MPI_COMM_WORLD);

  /* Reduce data */

  MPITransfer(RES2, buff, sz, 0, MPI_METH_RED);

  /* Move data to original block */

  if (mpiid == 0)
    {  
      /* Copy data */
      
      memcpy(RES2, buff, sz*sizeof(double));
  
      /* Free buffer */

      Mem(MEM_FREE, buff);     
    }

  /* Synchronise */

  MPI_Barrier(MPI_COMM_WORLD);

  /* Broadcast data to other tasks */

  MPITransfer(RES2, NULL, sz, 0, MPI_METH_BC);

  /* Synchronise */

  MPI_Barrier(MPI_COMM_WORLD);

  /***************************************************************************/

  fprintf(out, "OK.\n\n");

  /* Stop timers */

  StopTimer(TIMER_MPI_OVERHEAD);
  StopTimer(TIMER_MPI_OVERHEAD_TOTAL);

#endif
}

/*****************************************************************************/
