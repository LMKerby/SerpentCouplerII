/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : runfinix.c                                     */
/*                                                                           */
/* Created:       2013/03/27 (VVa)                                           */
/* Last modified: 2016/02/17 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Runs Finix for steady state or transient                     */
/*                                                                           */
/* Comments:  - Transientti pitää vielä kirjoittaa                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"
#include "./FINIX/initial.h"
#include "./FINIX/transient.h"
#include "./FINIX/aux_functions.h"

#define FUNCTION_NAME "RunFinix:"

/*****************************************************************************/

void RunFinix(long fib, long fpe)
{
  long tbi, prev, tme, nt, maxt, nz, nr, i, j;
  double tming, tmaxg, tmin, tmax;
  Rod *rod;
  Boundary_conditions *bc;
  Results *results, *results_prev;
  Options *options;
  char **err=NULL;

  /* Check steady state */

  if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {

      /* Criticality source simulation */

      /* Get Finix pointers for steady state */

      /* Get steady state Finix pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
      bc = (Boundary_conditions *)((long)RDB[fib + FINIX_PTR_BC]);

      /* Solve steady state solution */

      fprintf(out, "Solving FINIX steady state for rod %s\n", 
	      GetText(fib + FINIX_PTR_UNI_NAME));

      err = finix_solve_initial_steady_state(rod, bc, results, options);

      /* Handle errors */

      if (err != NULL)
	{
	  finix_printf_err(err);
	  Die(FUNCTION_NAME, "Error in executing FINIX");
	}

      finix_free_err(&err);

    }
  else
    {

      /* Time dependent simulation */

      /* Get transport time bin limits */

      tming = RDB[DATA_TIME_CUT_TMIN];

      /* Get current time bin index */

      nt = (long)RDB[DATA_DYN_TB];

      /* Get number of time bins */

      maxt = (long)RDB[DATA_DYN_NB];

      /* Check if this is the last transport time interval */

      if (nt == maxt - 1)
	{
	  /* If this is the last transport time interval, solve only until */
	  /* end of this interval */
      
	  tmaxg = RDB[DATA_TIME_CUT_TMAX];

	}
      else
	{
	  /* This is not the last transport time interval, solve also for */
	  /* next time interval to have a better initial guess for that   */

	  /* Get pointer to time bin structure */

	  tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
	  CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

	  /* Get pointer to bins */

	  tme = (long)RDB[tme + TME_PTR_BINS];
	  CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);    

	  /* Get time cut-off from next bin  */

	  tmaxg = RDB[tme + nt + 2];

	}
     
      

      /* Get pointer to first ifc time bin */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

      /* Loop over time bins */

      while (tbi > VALID_PTR)
	{

	  /* If this IFC time bin ends before current transport time interval */
	  /* get next IFC time interval */

	  if (RDB[tbi + IFC_FUEP_T_TMAX] < tming + 1E-15)
	    {
	      /* Next time bin */
	      
	      tbi = NextItem(tbi);

	      /* Cycle loop */

	      continue;
	    }

	  /* If this IFC time bin begins after the requested end time          */
	  /* We have looped over all relevant ifc bins, break out of loop      */

	  if (RDB[tbi + IFC_FUEP_T_TMIN] > tmaxg - 1E-15)
	    break;
 
	  /* Get Finix pointers for the time step*/
	  
	  rod = (Rod *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_ROD]);
	  results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	  options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
	  bc = (Boundary_conditions *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BC]);

	  /* Get temperature at the beginning of time step if available */

	  if ((prev = PrevItem(tbi)) > VALID_PTR)
	    {
	      results_prev = 
		(Results *)((long)RDB[prev + IFC_FUEP_FINIX_PTR_RESULTS]);
	    }
	  else
	    {
	      /* This is the first step, get initial steady state temperature */
	      /* via steady state calculation*/
	  
	      results_prev = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);

	    }

	  /* Set temperature distribution at beginning of step to equal that */
	  /* of the end of previous step */

	  nz = options->axial_nodes;
	  nr = options->pellet_radial_nodes + options->clad_radial_nodes;

	  for (i = 0 ; i < nz ; i++)
	    for (j = 0 ; j < nr ; j++)
	      results->temperature[i][j] = results_prev->temperature[i][j];

	  /* Get IFC time bin limits*/

	  tmin = RDB[tbi + IFC_FUEP_T_TMIN];
	  tmax = RDB[tbi + IFC_FUEP_T_TMAX];

	  /* Set up boundary conditions */

	  /* Set beginning time */

	  bc->time = tmin;

	  /* Set timestep length */

	  bc->dt = tmax - tmin;

	  /* Boundary conditions should be updated before runfinix.c if needed */

	  /* Run Finix for the current step */

	  fprintf(out,"Solving FINIX transient from %E to %E (%E s)\n", tmin, tmax, tmax-tmin);

	  /* Solve transient */

	  finix_solve_transient(bc->dt, rod, bc, results, options);

	  /* Handle errors */

	  if (err != NULL)
	    {
	      finix_printf_err(err);
	      Die(FUNCTION_NAME, "Error in transient FINIX solution");
	    }

	  finix_free_err(&err);

	  /* Next time bin */

	  tbi = NextItem(tbi);
	}

    }
}

#endif

/*****************************************************************************/
