/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : rendezvous.c                                   */
/*                                                                           */
/* Created:       2011/12/22 (JLe)                                           */
/* Last modified: 2014/11/27 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Collect source size, weight etc. data from MPI parallel      */
/*              tasks in criticality source mode                             */
/*                                                                           */
/* Comments: - Needed for MPI reproducibility                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Rendezvous:"

/*****************************************************************************/

void Rendezvous(long *nsrc0, double *wgt0)
{
#ifdef OLD_HIST
#ifdef MPI

  long n, ntot, *ntsk, *buf, ptr, i, i0, sz1, sz2, nhist, pth, id;
  double totwgt, *src1, *src2;

  /* Check mode */

  if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
    return;

  Die(FUNCTION_NAME, "Tää ei toimi uuden historiarakenteen kanssa");

  /* Check number of mpi tasks */

  if (mpitasks == 1)
    return;

  /* Start timers */

  StartTimer(TIMER_MPI_OVERHEAD);
  StartTimer(TIMER_MPI_OVERHEAD_TOTAL);

  /* Allocate memory for task-wise sizes */

  ntsk = (long *)Mem(MEM_ALLOC, mpitasks, sizeof(long));
  buf = (long *)Mem(MEM_ALLOC, mpitasks, sizeof(long));

  /* Put task-wise value */

  buf[mpiid] = *nsrc0;

  /* Reduce data */

  MPI_Barrier(MPI_COMM_WORLD);
  if (MPI_Reduce(buf, ntsk, mpitasks, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD)
      != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");

  /* Broadcast data */

  MPI_Barrier(MPI_COMM_WORLD);
  if (MPI_Bcast(ntsk, mpitasks, MPI_LONG, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Die(FUNCTION_NAME, "MPI Error");    

  /* Calculate total size */
      
  ntot = 0;
  for (n = 0; n < mpitasks; n++)
    ntot = ntot + ntsk[n];
    
  /* Free buffer */

  Mem(MEM_FREE, buf);

  /* Get size of particle and history data blocks */

  sz1 = PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE;
  sz2 = HIST_BLOCK_SIZE - LIST_DATA_SIZE;

  /* Get size of history array */

  if ((nhist = (long)RDB[DATA_HIST_LIST_SIZE]) < 0)
    nhist = 0;

  /* Allocate memory for histories */

  src1 = Mem(MEM_ALLOC, (sz1 + nhist*sz2)*ntot, sizeof(double));
  
  /* Calculate starting point */

  i0 = 0;
  for (i = 0; i < mpiid; i++)
    i0 = i0 + ntsk[i]*(sz1 + nhist*sz2);

  /* Reset OpenMP id */

  id = 0;

  /* Loop over distribution and read data into block */
  
  while (1 != 2)
    {
      /* Pointer to source distribution */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Pointer to first after dummy */
      
      if ((ptr = NextItem(ptr)) < VALID_PTR)
	break;

      /* Remove particle from source */

      RemoveItem(ptr);

      /* Copy particle data */

      memcpy(&src1[i0], &RDB[ptr + LIST_DATA_SIZE], sz1*sizeof(double));

      /* Update pointer */

      i0 = i0 + sz1;

      /* Pointer to history data */

      pth = (long)RDB[ptr + PARTICLE_PTR_HIST];

      /* Loop over histories */

      for (i = 0; i < nhist; i++)
	{
	  /* Check pointer */

	  CheckPointer(FUNCTION_NAME, "(pth1)", DATA_ARRAY, pth);

	  /* Copy history data */

	  memcpy(&src1[i0], &RDB[pth + LIST_DATA_SIZE], sz2*sizeof(double));

	  /* Update pointer */

	  i0 = i0 + sz2;

	  /* Next */

	  pth = NextItem(pth);
	}

      /* Put particle in stack */

      ToStack(ptr, id++);

      /* Check OpenMP id */

      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
	id = 0;
    }

  /* Allocate memory for temporary data */

  if (mpiid == 0)
    src2 = Mem(MEM_ALLOC, (sz1 + nhist*sz2)*ntot, sizeof(double));
  else
    src2 = NULL;

  /* Reduce data */

  MPI_Barrier(MPI_COMM_WORLD);
  MPITransfer(src1, src2, (sz1 + nhist*sz2)*ntot, 0, MPI_METH_RED);

  /* Move data back to original */

  if (mpiid == 0)
    memcpy(src1, src2, (sz1 + nhist*sz2)*ntot*sizeof(double));

  /* Free temporary array */

  if (mpiid == 0)
    Mem(MEM_FREE, src2);

  /* Broadcast data */

  MPI_Barrier(MPI_COMM_WORLD);
  MPITransfer(src1, NULL, (sz1 + nhist*sz2)*ntot, 0, MPI_METH_BC);

  /* Reset pointer and  OpenMP id */

  i0 = 0;
  id = 0;

  /* Read data back to source */

  for (n = 0; n < ntot; n++)
    {
      /* Get new particle from stack */

      ptr = FromStack(PARTICLE_TYPE_NEUTRON, id++);

      /* Check OpenMP id */

      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
	id = 0;

      /* Get pointer to history data (done here to avoid overwrite) */

      pth = (long)RDB[ptr + PARTICLE_PTR_HIST];      

      /* Copy particle data */

      memcpy(&WDB[ptr + LIST_DATA_SIZE], &src1[i0], sz1*sizeof(double));

      /* Put pointer to history data */

      WDB[ptr + PARTICLE_PTR_HIST] = (double)pth;

      /* Update pointer */

      i0 = i0 + sz1;

      /* Loop over histories */

      for (i = 0; i < nhist; i++)
	{
	  /* Check pointer */

	  CheckPointer(FUNCTION_NAME, "(pth2)", DATA_ARRAY, pth);

	  /* Copy history data */

	  memcpy(&WDB[pth + LIST_DATA_SIZE], &src1[i0], sz2*sizeof(double));

	  /* Update pointer */

	  i0 = i0 + sz2;

	  /* Next */

	  pth = NextItem(pth);
	}      

      /* Put particle back to source */
      
      AddItem(DATA_PART_PTR_SOURCE, ptr);
    }

  /* Free memory */

  Mem(MEM_FREE, src1);
  Mem(MEM_FREE, ntsk);

  /* Calculate number of particles and total weight */
  
  totwgt = 0.0;
  n = 0;

  /* Pointer to first after dummy */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  ptr = NextItem(ptr);

  /* Loop over remaining */

  while (ptr > VALID_PTR)
    {
      /* Add to counter and weight */

      n++;
      totwgt = totwgt + RDB[ptr + PARTICLE_WGT];

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Check */

  if (n != ntot)
    Die(FUNCTION_NAME, "Error in count");

  /* Put values */

  *nsrc0 = ntot;
  *wgt0 = totwgt;

  MPI_Barrier(MPI_COMM_WORLD);

  /* Stop timers */

  StopTimer(TIMER_MPI_OVERHEAD);
  StopTimer(TIMER_MPI_OVERHEAD_TOTAL);

#endif
#endif
}

/*****************************************************************************/
