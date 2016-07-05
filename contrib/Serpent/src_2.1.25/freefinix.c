/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : freefinix.c                                    */
/*                                                                           */
/* Created:       2014/11/05 (VVa)                                           */
/* Last modified: 2016/02/17 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Frees all finix arrays                                       */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "FreeFinix:"

/*****************************************************************************/

void FreeFinix()
{
  long fib, fpe, ifc, tbi;
  Rod *rod;
  Boundary_conditions *bc;
  Scenario *scenario;
  Results *results;
  Options *options;

  /* Check that some finix pins are defined */

  if ((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  /* MPI-id 0 is the task that has allocated these */
  /* let it free them */

  if (mpiid > 0)
    return;
  
  /* Loop over finix blocks to free arrays from them */

  while(fib > VALID_PTR)
    {
      /* Get Finix pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
      bc = (Boundary_conditions *)((long)RDB[fib + FINIX_PTR_BC]);
      scenario = (Scenario *)((long)RDB[fib + FINIX_PTR_SCENARIO]);

      /* Free arrays */

      finix_data_structures_destruct(rod, bc, scenario, results, options);

      /* Next finix pin*/

      fib = NextItem(fib);
    }

  if ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
    {

      /* Get pointer to first Finix block */

      fib = (long)RDB[DATA_PTR_FIN0];

      /* Get pointer to interface */

      ifc = (long)RDB[fib + FINIX_PTR_IFC];

      /* Pointer to FUEP BLOCK */

      fpe = (long)RDB[ifc + IFC_PTR_FUEP];

      /* Loop over pins to free arrays from interface */

      while (fpe > VALID_PTR)
	{

	  /* Loop over time steps */

	  tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

	  while(tbi > VALID_PTR)
	    {

	      /* Get Finix pointers for the time step*/
	  
	      rod = (Rod *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_ROD]);
	      scenario = (Scenario *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_SCENARIO]);
	      results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	      options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
	      bc = (Boundary_conditions *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BC]);

	      /* Free arrays */

	      finix_data_structures_destruct(rod, bc, scenario, results, options);

	      /* Next time bin */

	      tbi = NextItem(tbi);
	    }

	  /* Next rod */

	  fpe = NextItem(fpe);
	}
    }

}	  

#endif

/*****************************************************************************/
