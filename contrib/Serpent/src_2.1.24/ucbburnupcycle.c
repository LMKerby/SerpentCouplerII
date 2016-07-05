/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ucbburnupcycle.c                               */
/*                                                                           */
/* Created:       2012/07/05 (JLe)                                           */
/* Last modified: 2012/08/30 (JLe)                                           */
/* Version:       2.1.8                                                      */
/*                                                                           */
/* Description: UC Berkeley version of burnupcycle.c                         */
/*                                                                           */
/* Comments: Stripped-down version of BurnupCycle() for testing and          */
/*           development purposes, from 7.5.2012 (2.1.8) version of the main */
/*           routine.                                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UCBBurnupCycle:"

/*****************************************************************************/

void UCBBurnupCycle()
{
#ifdef EI_TOIMI_ENAA

  long dep, type, step, steps, mat, uni1, uni2;
  char tmpstr[MAX_STR];

  /* Remove binary work file */

  sprintf(tmpstr, "%s.wrk", GetText(DATA_PTR_INPUT_FNAME));
  remove(tmpstr);

  /* Reset predictor-corrector mode (the pc-loop was removed and */
  /* this is just to make sure that the code doesn't think its   */
  /* running a corrector calculation) */

  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_NONE;

  /* First loop is over intervals, each interval corresponds to a      */
  /* "dep"-card in the input. Step type (day/burnup) and normalization */ 
  /* may change between the steps */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Copy normalization from dep-structure to a variable used */
      /* in the transport and burnup routines */

      WDB[DATA_PTR_NORM] = RDB[dep + DEP_HIS_PTR_NORM];

      /* Get step type and number of steps */

      type = (long)RDB[dep + DEP_HIS_STEP_TYPE];
      steps = (long)RDB[dep + DEP_HIS_N_STEPS];

      /* Put type */

      WDB[DATA_BURN_STEP_TYPE] = (double)type;

      /* If this is the last interval (pointer to next interval is NULL), */
      /* add one more step to run transport calculation for the final     */
      /* composition */

      if (NextItem(dep) < VALID_PTR)
	steps++;

      /* Second loop is over steps listed in the dep-card */

      for (step = 0; step < steps; step++)
	{
	  /*******************************************************************/
	  
	  /***** Experimentation *********************************************/

	  /* Material compositions are written into a binary work file */
	  /* using StoreComposition(), that takes material pointer and */
	  /* a numerical index as arguments. Burnup step counter is    */
	  /* used here as the index, but it can be any number that     */
	  /* identifies the stored composition (iteration index, etc.) */
	  /* Function RetrieveComposition() retrieves the composition  */
	  /* of a previously stored material. Nothing is removed from  */
	  /* the file, so the index-pointer pair is unique, and if two */
	  /* compositions of the same material are stored with the     */
	  /* same index, only the first one can be retrieved.          */

	  /* Loop over materials and write compositions of burnable */
	  /* materials into binary work file */

	  mat = (long)RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Check burn-flag and write composition */

	      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
		StoreComposition(mat, (long)RDB[DATA_BURN_STEP]);

	      /* Next material */

	      mat = NextItem(mat);
	    }

	  /* Retrieve initial composition at 5th step (0 = initial) */

	  if ((long)RDB[DATA_BURN_STEP] == 5)
	    {
	      /* Loop over materials */

	      mat = (long)RDB[DATA_PTR_M0];
	      while (mat > VALID_PTR)
		{
		  /* Material name can be accessed with GetText() */

		  fprintf(out, "Retrieving material %s\n", 
			  GetText(mat + MATERIAL_PTR_NAME));

		  /* Check burn-flag and write composition */
		  
		  if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
		    RetrieveComposition(mat, 0);
		  
		  /* Next material */
		  
		  mat = NextItem(mat);
		}
	    }

	  /* Find pointers to universes named "1" and "2" (see the */
	  /* example input file) */

	  uni1 = (long)RDB[DATA_PTR_U0];
	  uni1 = SeekListStr(uni1, UNIVERSE_PTR_NAME, "1");

	  uni2 = (long)RDB[DATA_PTR_U0];
	  uni2 = SeekListStr(uni2, UNIVERSE_PTR_NAME, "2");

	  /* Swap the two universes */

	  SwapUniverses(uni1, uni2);

	  /* These routines must be called to re-create the geometry */
      
	  CellCount(-1, -1, 0, 1);
	  MaterialVolumes();
	  SumDivCompositions();

	  /* This checks material volumes by MC (it's probably a good idea */
	  /* to call this and confirm that the calculation sees correct    */
	  /* volumes and masses, since incorrect volume and mass is the    */
	  /* most probable cause for errors in burnup calculation.)        */

	  WDB[DATA_VOLUME_MC_NMAX] = 100000000;
	  WDB[DATA_VOLUME_CALCULATION_MODE] = (double)YES;
	  VolumesMC();
	  WDB[DATA_VOLUME_CALCULATION_MODE] = (double)NO;

	  /***** End experimentation *****************************************/

	  /*******************************************************************/

	  /* Step type is predictor, since predictor-corrector calculation */
	  /* was disabled */
	  
	  WDB[DATA_BURN_STEP_PC] = PREDICTOR_STEP;
              
	  /* This routine pre-calculates cross sections, clears results  */
	  /* from previous steps, etc. and makes everything ready to for */
	  /* the transport calculation */

	  PrepareTransportCycle();

	  /* Print material compositions (<input>.bumat<n> -files) */
	  
	  PrintCompositions((long)RDB[DATA_BURN_STEP]);
              
	  /* Run transport calculation cycle if step type is not decay, */
	  /* otherwise just print the output */

	  if ((type != DEP_STEP_DEC_STEP) && (type != DEP_STEP_DEC_TOT))
	    TransportCycle();
	  else
	    MatlabOutput();

	  /* K-eff from the transport calculation can be accessed like this */

	  fprintf(out, "k-eff = %1.5f +/- %1.5f\n", 
		  RDB[DATA_BURN_PREV_KEFF], RDB[DATA_BURN_PREV_DKEFF]);
	  
	  /* Start burnup timers - these are used for timing the burnup */
          /* calculation loop (BURNUP_CYCLE_TIME in the _res.m output)  */
	  
	  ResetTimer(TIMER_BURNUP);
	  StartTimer(TIMER_BURNUP);
	  StartTimer(TIMER_BURNUP_TOTAL);
	  
	  /* Print depletion output - The depletion output is done in two  */
	  /* stages: 1) write a binary file with all nuclide densities and */
	  /* 2) read the file and print out the nuclide compositions given */
	  /* in the inventory list. The file is used for output only, and  */
	  /* the idea is that any composition can be requested afterwards  */
	  /* without re-running the entire simulation (the -rdep command   */
	  /* line option */
	  
	  fprintf(out, "Writing depletion output...\n");
          
	  /* 1) Write binary depletion file */
	  
	  WriteDepFile();       
	  
	  /* 2) Read the compositions and print depletion output */
	      
	  PrintDepOutput();     
	  
	  fprintf(out, "OK.\n\n");
	  
	  /* Break here if final step (the remainder of the loop is for */
	  /* preparing and running the next depletion calculation) */

	  if (step == RDB[dep + DEP_HIS_N_STEPS])
	    break;
              
	  /* Calculate coefficients for the fit to xs/flux/power */
	  /* (this is predictor-corrector-stuff) */

	  DepletionPolyFit(dep, step);
	  
	  /* Set depletion step size (converts step given in units of */
	  /* burnup or days into seconds) */

	  SetDepStepSize(dep, step);

	  /* The step length can be overriden using these two variables. */
	  /* Time step is the one actually used in the calculation,      */
	  /* burnup step is for printing only. */

	  /*
	  WDB[DATA_BURN_TIME_INTERVAL] = time_step_in_seconds;
	  WDB[DATA_BURN_BURNUP_INTERVAL] = burnup_step_in_MWDkgU;
	  */

	  /* Burnup calculation - this is the main depletion solver    */
          /* routine that loops over burnable materials, depletes them */
	  /* and updates the composition */
	  
	  BurnMaterials(dep, step);
          
	  /* Collect material compositions from MPI parallel tasks */
	  /* (in MPI calculation the materials are divided between */
	  /* the parallel tasks, so all results must be gathered   */
	  /* and re-distributed before running the next transport  */
	  /* simulation) */

	  CollectBurnData();

	  /* Stop burnup timers (the two timers started above) */
	  
	  StopTimer(TIMER_BURNUP);
	  StopTimer(TIMER_BURNUP_TOTAL);
          
	  /* Add to number of predictor and corrector cycles */
	  
	  WDB[DATA_BURN_PRED_STEP] = RDB[DATA_BURN_PRED_STEP] + 1.0;
	  
	  /* Update cumulative burnup and time */
	  
	  WDB[DATA_BURN_CUM_BURNTIME] = RDB[DATA_BURN_CUM_BURNTIME] 
	    + WDB[DATA_BURN_TIME_INTERVAL];
	  WDB[DATA_BURN_CUM_BURNUP] = RDB[DATA_BURN_CUM_BURNUP]
	    +  WDB[DATA_BURN_BURNUP_INTERVAL];
	  
	  /* Update burnup step (total, not in this interval) */
	  
	  WDB[DATA_BURN_STEP] = RDB[DATA_BURN_STEP] + 1.0;
	}

      /* Next interval */

      dep = NextItem(dep);
    }

  /* Check total time */

  if (RDB[DATA_BURN_CUM_BURNTIME] == 0.0)
    Die(FUNCTION_NAME, "No burnup calculation performed");  

#endif
}

/*****************************************************************************/
