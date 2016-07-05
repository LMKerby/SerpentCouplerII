/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : distributefinix.c                              */
/*                                                                           */
/* Created:       2013/10/11 (VVa)                                           */
/* Last modified: 2016/02/17 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Sets up FINIX time bin data to interface structures          */
/*              this can not be done in processfinix as the interface        */
/*              structures are not present then                              */
/*                                                                           */
/* Comments:   -Needs both, the FINIX blocks and the interface blocks to be  */
/*              created prior to this routine.                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"
#include "./FINIX/finix_initialization.h"
#include "./FINIX/aux_functions.h"
#include "./FINIX/initial.h"
#include "./FINIX/finix_output.h"

#define FUNCTION_NAME "DistributeFinix:"

#define FILE_TYPE_ROD      0
#define FILE_TYPE_OPTIONS  1
#define FILE_TYPE_SCENARIO 2

/*****************************************************************************/

void DistributeFinix()
{
  long i, j, fib, ifc;
  long fpe, found, nu, ptr, tbi, count;
  double val, frac, last;
  Rod *rod;
  Boundary_conditions *bc;
  Scenario *scenario;
  Results *results;
  Options *options;
  char **err=NULL;

  if ((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  fprintf(out,"Distributing FINIX data to timesteps\n");

  /* Loop over FINIX definitions */

  while (fib > VALID_PTR)
    {
      /* Get Finix pointers */

      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);
      scenario = (Scenario *)((long)RDB[fib + FINIX_PTR_SCENARIO]);
      bc = (Boundary_conditions *)((long)RDB[fib + FINIX_PTR_BC]);

      /* Get pointer to corresponding interface */

      ifc = (long)RDB[fib + FINIX_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

      /* Find correct pin from interface */
      /* Pointer to ifc-rod block */

      fpe = (long)RDB[ifc + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(fpe)", DATA_ARRAY, fpe);

      /* Loop over interface pins to find match */

      while (fpe > VALID_PTR)
	{
	  /* Get Finix pin uni */

	  found = 0;

	  /* Get number of universes in this pin */

	  nu = RDB[fpe + IFC_FUEP_N_UNI];

	  /* Get pointer to rod segments */

	  ptr = (long)RDB[fpe + IFC_FUEP_PTR_UNI_LIST];

	  /* Compare to the universes of the interface FPE block */

	  for(i=0; i < nu; i++)
	    {
	      if(CompareStr(ptr + i, fib + FINIX_PTR_UNI_NAME))
		found=1;
	    }

	  /* Check if found */

	  if (found == 1)
	    break;
	  else
	    fpe = NextItem(fpe);	 

	}

      /* Check if found */
      if (fpe < VALID_PTR)
	Die(FUNCTION_NAME, "Could not find universe");
      
      last =  0.0;
      count = 0; 

      /* Get time bin pointer */
	  
      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

      /* Loop over time bins to initialize a separate FINIX for each time step */

      while (tbi > VALID_PTR)
	{

	  /* Calculate fraction and print progress */
	  
	  frac = (double)(count++)/((double)ListSize(tbi));
	  
	  if (frac - last > 0.10)
	    {
	      fprintf(out, " %3.0f%% complete\n", 100.0*frac);
	      last = frac;
	    }


	  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	    {
	      /* In criticality source mode, just copy steady state pointers here */

	      WDB[tbi + IFC_FUEP_FINIX_PTR_ROD] = RDB[fib + FINIX_PTR_ROD];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_SCENARIO] = RDB[fib + FINIX_PTR_SCENARIO];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS] = RDB[fib + FINIX_PTR_RESULTS];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS] = RDB[fib + FINIX_PTR_OPTIONS];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_BC] = RDB[fib + FINIX_PTR_BC];

	      /* Next time bin (there should not be more than one) */

	      tbi = NextItem(tbi);

	      /* Cycle loop */

	      continue;
	    }

	  /* TODO: This initialization is mostly the same as in processfinix.c */
	  /* but has to be done for each FINIX entity. Currently there are     */
	  /* separate FINIXes for each time step */

	  /* Construct arrays */

	  rod = finix_rod_construct();
	  bc = finix_bc_construct();
	  scenario = finix_scenario_construct();
	  results = finix_results_construct();
	  options = finix_options_construct();

	  /* Store pointers */

	  WDB[tbi + IFC_FUEP_FINIX_PTR_ROD] = (double)((long)rod);
	  WDB[tbi + IFC_FUEP_FINIX_PTR_SCENARIO] = (double)((long)scenario);
	  WDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS] = (double)((long)results);
	  WDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS] = (double)((long)options);
	  WDB[tbi + IFC_FUEP_FINIX_PTR_BC] = (double)((long)bc);

	  /* Create finix_rod.inp for this rod */

	  WriteFinixInputFile(fib, FILE_TYPE_ROD);

	  /* Create finix_options.inp for this rod */

	  WriteFinixInputFile(fib, FILE_TYPE_OPTIONS);

	  /* Create finix_scenario.inp for this rod */

	  WriteFinixInputFile(fib, FILE_TYPE_SCENARIO);

	  /* Initialize data structures */

	  err = finix_initialize_data_structures(rod, bc, scenario, results, options);

	  /* Handle errors */

	  if (err != NULL)
	    {
	      finix_printf_err(err);
	      Die(FUNCTION_NAME, "Error in initialization of FINIX data structures");
	    }

	  finix_free_err(&err);

	  /* Put zero linear power to all axial segments */

	  for (i = 0 ; i < options->axial_nodes ; i++)
	    bc->linear_power[i] = 0.0;

	  /* Set radial nodes to equal-distance, rather than equal area */

	  for (i = 0 ; i < options->axial_nodes ; i++)
	    {

	      /* Pellet */
	      for (j = 0; j < options->pellet_radial_nodes ; j++)
		{
		  val = rod->pellet_inner_radius 
		    + (double)j/(double)(options->pellet_radial_nodes - 1)*
		    (rod->pellet_outer_radius - rod->pellet_inner_radius);

		  results->radial_node_position[i][j] = val;
		  results->radial_node_position_cold[i][j] = val;

		}

	      /* Clad */

	      for (; j < options->pellet_radial_nodes 
		     + options->clad_radial_nodes ; j++)
		{
		  val = rod->clad_inner_radius 
		    + (double)(j - options->pellet_radial_nodes)/
		    (double)(options->clad_radial_nodes - 1)*
		    (rod->clad_outer_radius - rod->clad_inner_radius);

		  results->radial_node_position[i][j] = val;
		  results->radial_node_position_cold[i][j] = val;

		}
	    }

	  /* Solve initial steady state (HZP) */

	  err = finix_solve_initial_steady_state(rod, bc, results, options);

	  /* Handle errors */

	  if (err != NULL)
	    {
	      finix_printf_err(err);
	      Die(FUNCTION_NAME, "Error in initialization of FINIX data structures");
	    }

	  finix_free_err(&err);

	  /* Next time bin */
	  
	  tbi = NextItem(tbi);

	}

      /* Next fuel rod */

      fib = NextItem(fib);
    }

  fprintf(out, " %3.0f%% complete\n\n", 100.0);

}

#endif

/*****************************************************************************/
