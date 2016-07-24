#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printcycleoutput.c                             */
/*                                                                           */
/* Created:       2011/04/03 (JLe)                                           */
/* Last modified: 2016/04/13 (TVi)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Prints cycle-wise data & info to standard output             */
/*                                                                           */
/* Comments:  - Kytketyssä laskennassa tätä pitää päivittää                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintCycleOutput:"

/*****************************************************************************/

void PrintCycleOutput()
{
  long i, cycles, skip, skip1, pop, ptr, tb, tbmax, tme;
  double estimt, tott, tming, tmaxg;
  char tmpstr[MAX_STR];

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Get cycle index, number of cycles and population */
  
  i = (long)RDB[DATA_CYCLE_IDX] + 1;

  /* Set number of inactive batches */

  if(RDB[DATA_USE_FSP] == (double)NO)
    {
      /* No fission source passing*/

      skip = (long)RDB[DATA_CRIT_SKIP];
      skip1 = skip;
    }
  else if ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
         && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    {
      /* First transportcycle with fission source passing */

      skip = (long)RDB[DATA_CRIT_SKIP];
      skip1 = skip;
    }
  else
    {
      /* Subsequent transportcycle with fission source passing */

      skip = (long)RDB[DATA_CRIT_SKIP];
      skip1 = (long)RDB[DATA_FSP_CRIT_SKIP];
    }

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {

#ifdef MPI_MODE2

      /* Check number of tasks */
  
      if (mpitasks > 1)
	{
	  /* Calculate number of particles per task */

	  cycles = (long)(RDB[DATA_CRIT_CYCLES]/((double)mpitasks));

	}
      else
	cycles = (long)RDB[DATA_CRIT_CYCLES];

#else
      
      cycles = (long)RDB[DATA_CRIT_CYCLES];

#endif

      pop = (long)RDB[DATA_CRIT_POP];
    }
  else
    {
      cycles = (long)RDB[DATA_SRC_BATCHES];
      pop = (long)RDB[DATA_SRC_POP];
    }

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      /* Get time bin index */

      tb = (long)RDB[DATA_DYN_TB];
      tbmax = (long)RDB[DATA_DYN_NB];

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);
      
      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme ptr)", DATA_ARRAY, tme);    

      /* Get transport time interval */

      tming = RDB[tme + tb];
      tmaxg = RDB[tme + tb + 1];
    }
  else
    {
      /* Criticality source mode  */
      /* Set interval to infinity */

      tb = 0;
      tbmax = 1;

      tming = -INFTY;
      tmaxg = INFTY;
    }

  /* Get estimater running time */

  EstimateRuntime();
  
  estimt = RDB[DATA_ESTIM_CYCLE_TIME];
  tott = RDB[DATA_ESTIM_TOT_TIME];

  /***************************************************************************/
	 
  /***** Print inactive cycle output *****************************************/

  if (i < skip + 1)
    {
      if (RDB[DATA_USE_FSP] == (double)NO)
	fprintf(out, "Inactive cycle %3ld / %3ld: ", i, skip);
      else
	fprintf(out, "Inactive cycle %3ld / %3ld: ", i + skip1 - skip, skip1);

      if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
	fprintf(out, "k-eff = %1.5f\n", RDB[DATA_WIELANDT_KP]);
      else
	fprintf(out, "k-eff = %1.5f\n", RDB[DATA_CYCLE_KEFF]);
    }

  /***************************************************************************/

  /***** Print active cycle output *******************************************/

  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT) &&
      (i == skip + 1))
    fprintf(out, "\n----- Begin active cycles -----\n\n");

  if (i > skip)
    {    
      fprintf(out, "------------------------------------------------------------\n");

      fprintf(out, "\nSerpent %s", CODE_VERSION);

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	fprintf(out, " -- Static criticality source simulation\n");
      else if (((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES) &&
	       ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES))
	fprintf(out, " -- Combined neutron / photon transport simulation\n");
      else if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
	{
	  if ((long)RDB[DATA_DYN_NB] == 1)
	    fprintf(out, " -- Neutron external source simulation\n");
	  else
	    fprintf(out, " -- Dynamic neutron external source simulation\n");
	}
      else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
	fprintf(out, " -- Photon transport simulation\n");

      if ((long)RDB[DATA_PTR_TITLE] > VALID_PTR)
	fprintf(out, "\nTitle: \"%s\"\n", GetText(DATA_PTR_TITLE));
      else
	fprintf(out, "\nInput file: \"%s\"\n", GetText(DATA_PTR_INPUT_FNAME));

      if ((long)RDB[DATA_COEF_CALC_IDX] > 0)
	fprintf(out, "\nCoefficient calculation: restart = %ld / %ld\n", 
		(long)RDB[DATA_COEF_CALC_RUN_IDX], 
		(long)RDB[DATA_COEF_CALC_TOT_RUNS]);
      else if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
	{
	  fprintf(out, "\nTransport calculation: step = %ld / %ld ", 
		 (long)RDB[DATA_BURN_STEP] + 1, 
		 (long)RDB[DATA_BURN_TOT_STEPS] + 1);
	  
	  if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE)
	    fprintf(out, "\n");
	  else if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
	    fprintf(out, "(predictor)\n");
	  
	  else if ((long)RDB[DATA_BURN_CI_MAXI] > 1)

	    fprintf(out, "(corrector %ld/%ld)\n",
                    (long)RDB[DATA_BURN_CI_I]+1, (long)RDB[DATA_BURN_CI_MAXI]);

	  else
	    fprintf(out, "(corrector)\n");
	  
	  if (RDB[DATA_INI_FMASS] > 0.0)
	    fprintf(out, "                       BU   = %1.2f MWd/kgU\n", 
		   RDB[DATA_BURN_CUM_BURNUP]);
	  
	  if (RDB[DATA_BURN_CUM_BURNTIME] == 0.0)    
	    fprintf(out, "                       time = 0.00 days\n");
	  else
	    fprintf(out, "                       time = %s\n", 
		   TimeIntervalStr(RDB[DATA_BURN_CUM_BURNTIME]));
	}

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	{
	  fprintf(out, "\nActive cycle %4ld / %ld  Source neutrons : %5ld\n\n", 
		   i - skip, cycles, (long)RDB[DATA_CYCLE_BATCH_SIZE]);

	  if(RDB[DATA_RUN_CC] == (double)YES)
	    {
	      fprintf(out, "Coupled calculation iteration:");
	      fprintf(out, "%3ld / %3ld \n\n",
		      (long)RDB[DATA_SOL_REL_ITER] + 1,
		      (long)RDB[DATA_SOL_REL_MAX_ITER]);
	    }
	}
      else
	{
	  fprintf(out, "\nSource batch %2ld / %ld (%ld histories per batch)\n", 
		  i, cycles, pop);
	 
	  
	  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
	    {
	      /* Print current time interval */

	      fprintf(out, "Time interval %ld / %ld from %6E s to %6E s\n", tb + 1, tbmax, tming, tmaxg);
	    }
	  
	}

      /* Print running time */
      
      fprintf(out, "Running time :                %9s\n", 
	      TimeStr((long)TimerVal(TIMER_RUNTIME)));

      if (tott > 0.0)
	{
	  if (TimerVal(TIMER_TRANSPORT_ACTIVE) > 2.0)
	    {
	      strcpy(tmpstr, TimeStr((long)tott));
	      fprintf(out, "Estimated running time :      %9s %9s\n", 
		     TimeStr((long)estimt), tmpstr);

	      strcpy(tmpstr, TimeStr((long)(tott- TimerVal(TIMER_RUNTIME))));
	      fprintf(out, "Estimated running time left : %9s %9s\n", 
		     TimeStr((long)(estimt - TimerVal(TIMER_RUNTIME))),
		     tmpstr);
	    }
	  else
	    {
	      fprintf(out, "Estimated running time :        -:--:--   -:--:--\n");
	      fprintf(out, "Estimated running time left :   -:--:--   -:--:--\n");
	    }
	}
      else
	{
	  if (TimerVal(TIMER_TRANSPORT_ACTIVE) > 2.0)
	    {
	      fprintf(out, "Estimated running time :      %9s\n", 
		     TimeStr((long)estimt));
	      fprintf(out, "Estimated running time left : %9s\n", 
		     TimeStr((long)(estimt - TimerVal(TIMER_RUNTIME))));
	    }
	  else
	    {
	      fprintf(out, "Estimated running time :        -:--:--\n");
	      fprintf(out, "Estimated running time left :   -:--:--\n");
	    }
	}

      fprintf(out, "\nEstimated relative CPU usage : %7.1f%%\n", 
	     100.0*TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
	     TimerVal(TIMER_TRANSPORT_CYCLE));

      /* K-eff estimates */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	{
	  ptr = (long)RDB[RES_ANA_KEFF];
	  fprintf(out, "\nk-eff (analog)    = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
		 Mean(ptr, 0), StdDev(ptr, 0),
		 Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
		 Mean(ptr, 0) + 1.96*StdDev(ptr, 0));
	  
	  if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == YES)
	    {
	      ptr = (long)RDB[RES_IMP_KEFF];
	      fprintf(out, "k-eff (implicit)  = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
		     Mean(ptr, 0), StdDev(ptr, 0),
		     Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
		     Mean(ptr, 0) + 1.96*StdDev(ptr, 0));
	    }
	  else
	    {
	      ptr = (long)RDB[RES_COL_KEFF];
	      fprintf(out, "k-eff (collision) = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
		     Mean(ptr, 0), StdDev(ptr, 0),
		     Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
		     Mean(ptr, 0) + 1.96*StdDev(ptr, 0));
	    }	    
	}
      else if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
	{
	  ptr = (long)RDB[RES_SRC_MULT];
	  fprintf(out, "\nMultiplication   = %1.4E (%1.5f)\n",
		  Mean(ptr, 0), RelErr(ptr, 0));

	  ptr = (long)RDB[RES_ANA_KEFF];
	  fprintf(out, "k-eff (analog)   = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
		  Mean(ptr, 0), StdDev(ptr, 0),
		  Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
		  Mean(ptr, 0) + 1.96*StdDev(ptr, 0));
	  
	  ptr = (long)RDB[RES_EXT_K];
	  fprintf(out, "k0    (source)   = %1.5f +/- %1.5f  [%1.5f  %1.5f]\n",
		  Mean(ptr, 0), 
		  StdDev(ptr, 0),
		  Mean(ptr, 0) - 1.96*StdDev(ptr, 0),
		  Mean(ptr, 0) + 1.96*StdDev(ptr, 0));

	  if (RDB[DATA_DYN_TMAX] != INFTY)
	    {
	      if ((long)RDB[DATA_DYN_NB] > 1)
		fprintf(out, "\nTime cut-off at %1.2E seconds, %ld intervals\n",
			RDB[DATA_DYN_TMAX], (long)RDB[DATA_DYN_NB]);
	      else
		fprintf(out, "\nTime cut-off at %1.2E seconds\n",
			RDB[DATA_DYN_TMAX]);
	    }

	}

      fprintf(out, "\n");

      /* Options */

      fprintf(out, "Options : ");

#ifdef DEBUG

      fprintf(out, "(DBG) ");

#endif

      if ((long)RDB[DATA_OPTI_REPLAY] == YES)
	fprintf(out, "(R) ");

      fprintf(out, "(O%ld) ", (long)RDB[DATA_OPTI_MODE]);

      if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
	fprintf(out, "(W) ");

      if ((long)RDB[DATA_ITER_MODE] == ITER_MODE_ALBEDO)
	fprintf(out, "(IA) ");

      if ((long)RDB[DATA_B1_CALC] == YES)
	fprintf(out, "(B1) ");

      /* Cells with union operator (this is temporary) */

      if ((long)RDB[DATA_N_UNION_CELLS] > 0)
	fprintf(out, "(U) ");

      if ((long)RDB[DATA_OPT_IMPL_CAPT] == YES)
	fprintf(out, "(IC) ");
      /*
      if ((long)RDB[DATA_OPT_IMPL_NXN] == YES)
	fprintf(out, "(IX) ");

      if ((long)RDB[DATA_OPT_IMPL_FISS] == YES)
	fprintf(out, "(IF) ");
      */
      if ((long)RDB[DATA_USE_URES] == YES)
	fprintf(out, "(UNR) ");

      if ((long)RDB[DATA_USE_DBRC] == YES)
	fprintf(out, "(DBRC) ");

      if ((long)RDB[DATA_TMS_MODE] != TMS_MODE_NONE)
	fprintf(out, "(TMS) ");
      
      if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] > -1)
	fprintf(out, "(XE) ");

      if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] > -1)
	fprintf(out, "(SM) ");

#ifdef MPI
      
      fprintf(out, "(MPI=%d) ", mpitasks);
      
#endif
      
#ifdef OPEN_MP

      fprintf(out, "(OMP=%ld) ", (long)RDB[DATA_OMP_MAX_THREADS]);

#endif

      if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)     
	{
	  if ((long)RDB[DATA_BURN_PRED_TYPE] == PRED_TYPE_CONSTANT)
	    fprintf(out, "(CE");
	  else
	    fprintf(out, "(LE");

          if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_CONSTANT)
	    fprintf(out, "/CE");
	  else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_LINEAR)
	    fprintf(out, "/LI");
	  else if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_QUADRATIC)
	    fprintf(out, "/QI");

	  if (((long)RDB[DATA_BURN_PRED_NSS] > 1) && 
	      ((long)RDB[DATA_BURN_CORR_NSS] < 1))
	    fprintf(out, " + %ld SS) ", (long)RDB[DATA_BURN_PRED_NSS]);
	  else if (((long)RDB[DATA_BURN_PRED_NSS] == 1) && 
		   ((long)RDB[DATA_BURN_CORR_NSS] > 1))
	    fprintf(out, " + %ld SS) ", (long)RDB[DATA_BURN_CORR_NSS]);
	  else if (((long)RDB[DATA_BURN_PRED_NSS] > 1) && 
		   ((long)RDB[DATA_BURN_CORR_NSS] > 1))
	    fprintf(out, " + %ld/%ld SS) ", (long)RDB[DATA_BURN_PRED_NSS],
		   (long)RDB[DATA_BURN_CORR_NSS]);
	  else
	    fprintf(out, ") ");
	}

      if ((long)RDB[DATA_RUN_CC] == YES)
	fprintf(out, "(CC) ");

      fprintf(out, "\n");

      fprintf(out, "------------------------------------------------------------\n");
    }

  /* Check if completed */

  if ((i == cycles + skip) && 
      ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_DYN))
    {
      fprintf(out, "\nTransport cycle completed in %s.\n\n", 
	      TimeIntervalStr(TimerVal(TIMER_TRANSPORT)));
      
      PrintTMSDiagnostics();
    }
  
  /* For EDo (18.1.2014) */
  
  DiffCoefED(5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  
  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
