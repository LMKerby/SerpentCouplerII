/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iteratecc.c                                    */
/*                                                                           */
/* Created:       2014/07/07 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Iterates external/internal solvers for coupled calculation   */
/*                                                                           */
/* Comments: -Handles solution relaxation and related calculations           */
/*           -Also transfers updated interfaces between MPI-tasks            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IterateCC:"

/*****************************************************************************/

void IterateCC()
{
  long id, mat;

#ifdef MPI
  long loc0, sz;
#endif

  /* Return if not running coupled calculation */

  if(RDB[DATA_RUN_CC] == NO)
    return;

  /* If on non-first predictor and running SIE, return */
  if((((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) && 
     !((long)RDB[DATA_BURN_STEP] == 0)) && (RDB[DATA_BURN_SIE] == (double)YES))
    return;

  /* Calculate relaxation factor alpha */

  CalculateRelAlpha();

  /* Calculate population size for next iteration */

  CalculateRelPopSize();

  /* Relax the tallied power distribution */

  RelaxInterfacePower();

  /* Intra-step relaxation of tallied transmutation cross sections */

  /* Get OpenMP id */

  id = OMP_THREAD_NUM;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
	{

	  /* Calculate momentary transmutation XS*/

	  CalculateTransmuXS(mat, id);

	  /* Relax transmutation XS inside burnup step */

	  RelaxTransmuXS(mat, id);

	}

      /* Next material */

      mat = NextItem(mat);

    }

  /* Only call the coupled solvers from MPI task 0 */

  if(mpiid == 0)
    {

      IterateFinix();

      /* IterateCOSY();*/

      IterateExternal();

      IterateInternal();
    }

  /* Broadcast updated interfaces from 0 to other tasks */
#ifdef MPI

  /* Loop over interfaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];

  while(loc0 > VALID_PTR)
    {
      sz = (long)RDB[loc0 + IFC_MEM_SIZE];

      /* Synchronise */

      MPI_Barrier(MPI_COMM_WORLD);

      /* Broadcast data from task 0 to other tasks */

      MPITransfer(&WDB[loc0], NULL, sz, 0, MPI_METH_BC);

      /* Synchronise */

      MPI_Barrier(MPI_COMM_WORLD);

      /* Next interface */

      loc0 = NextItem(loc0);

    }

#endif

  /* Increase iteration number */

  WDB[DATA_SOL_REL_ITER] = RDB[DATA_SOL_REL_ITER] + 1.0;

  /* Check stopping criterion */

  StopCCIter();

  return;
}
