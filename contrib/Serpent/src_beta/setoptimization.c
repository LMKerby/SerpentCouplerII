/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setoptimization.c                              */
/*                                                                           */
/* Created:       2011/07/17 (JLe)                                           */
/* Last modified: 2016/02/16 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: - Set various options based on optimization                  */
/*                                                                           */
/* Comments: - Tää aliohjelma tekee oikeastaan kaikkea muutakin kuin         */
/*             asettaa optimoinnin, mm. korjaa käyttäjän tekemiä valintoja   */
/*             jne... Nimen voisi muuttaa.                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetOptimization:"

/*****************************************************************************/

void SetOptimization()
{
  char tmpstr[MAX_STR];

  /***************************************************************************/

  /***** Set options etc. ****************************************************/

  /* Track plotter */

  if ((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_TRACKS)
    {
      /* Check number of tracks */

      if (RDB[DATA_TRACK_PLOTTER_HIS] < 1.0)
	Die(FUNCTION_NAME, "Invalid number of tracks");

      /* Reset burnup mode */

      WDB[DATA_BURNUP_CALCULATION_MODE] = NO;

      /* Set number of histories */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	{
	  WDB[DATA_CRIT_POP] = RDB[DATA_TRACK_PLOTTER_HIS];
	  WDB[DATA_CRIT_CYCLES] = 0.0;
	  WDB[DATA_CRIT_SKIP] = 1.0;
	}
      else
	{
	  WDB[DATA_SRC_POP] = RDB[DATA_TRACK_PLOTTER_HIS];
	  WDB[DATA_SRC_BATCHES] = 1.0;
	  WDB[DATA_CRIT_SKIP] = 0.0;
	}

      /* Set length of history */

      WDB[DATA_HIST_LIST_SIZE] = RDB[DATA_TRACK_PLOT_NHIS];

      /* Set flag */

      SetOption(DATA_EVENT_RECORD_FLAGS, RECORD_EVENT_PLOTTER);

      /* Reset minimum collision frequency */

      WDB[DATA_CFE_N_MIN_L] = INFTY;
      WDB[DATA_CFE_G_MIN_L] = INFTY;
    }

  /* RIA Simulation */

  if ((long)RDB[DATA_PTR_RIA0] > VALID_PTR)
    {
      /* Set simulation mode to criticality */

      WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_CRIT;

      /* Check parameters */

      if ((long)RDB[DATA_SRC_POP] < 1)
	Error(0, "Source batch size must be defined");
      if ((long)RDB[DATA_SRC_BATCHES] < 1)
	Error(0, "Number of source batches cycles must be defined");

      /* Write source to file */

      sprintf(tmpstr, "%s.src", GetText(DATA_PTR_INPUT_FNAME));
      WDB[DATA_PTR_CRIT_SRC_DET] = (double)PutText(tmpstr);
    }

  /* Add k-eff iteration cycles to skip cycles */

  if ((long)RDB[DATA_ITER_MODE] != ITER_MODE_NONE)
    {
      if ((long)RDB[DATA_ITER_FIX] == YES)
	WDB[DATA_CRIT_SKIP] = 2.0*RDB[DATA_CRIT_SKIP] + RDB[DATA_ITER_NCYC];
      else
	WDB[DATA_CRIT_SKIP] = RDB[DATA_CRIT_SKIP] + RDB[DATA_ITER_NCYC];
    }

  /* Check that simulation mode is set */

  if ((long)RDB[DATA_SIMULATION_MODE] < 0)
    Error(0, "Simulation mode must be set using \"pop\", \"nps\" or \"dyn\"");
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Check parameters */

      if ((long)RDB[DATA_CRIT_POP] < 1)
	Error(0, "Population size must be defined");
      if ((long)RDB[DATA_CRIT_CYCLES] + (long)RDB[DATA_CRIT_SKIP] < 1)
	Error(0, "Number of criticality cycles must be defined");
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
    {
      /* Check parameters */

      if ((long)RDB[DATA_SRC_POP] < 1)
	Error(0, "Source batch size must be defined");
      if ((long)RDB[DATA_SRC_BATCHES] < 1)
	Error(0, "Number of source batches cycles must be defined");

      /* Reset number of skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* External source mode + precursor tracking = DELDYN mode */
      /* External source mode + coupled calculation = DYN mode */
      
      if ((long)RDB[DATA_PTR_PREC_DET] > VALID_PTR) 
	WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_DELDYN;
      else if (RDB[DATA_RUN_CC] == (double)YES) 
	WDB[DATA_SIMULATION_MODE] = (double)SIMULATION_MODE_DYN;
    }
  else
    Die(FUNCTION_NAME, "Invalid simulation mode");

  /* Include total list in mode 0 */

  WDB[DATA_OPTI_MODE0_INCLUDE_TOTAL] = 1.0;

  /* Set unionized energy grid thinning tolerance */

  if (RDB[DATA_ERG_TOL] < 0.0)
    {
      if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
	WDB[DATA_ERG_TOL] = 0.0;
      else
	WDB[DATA_ERG_TOL] = 5E-5;
    }

  /* Check optimization mode */

  switch ((long)RDB[DATA_OPTI_MODE])
    {
    case 1:
      {
	WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)NO;
	WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)NO;
	WDB[DATA_OPTI_IMPLICIT_RR] = (double)NO;
	WDB[DATA_OPTI_GC_CALC] = (double)NO;

	if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
	  WDB[DATA_OPTI_MG_MODE] = (double)NO;	

	if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
	  Note(0, "Option 'set xscalc 2' ignored in optimization mode 1");

	WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)NO;

	break;
      }
    case 2:
      {
	WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)YES;
	WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)NO;
	WDB[DATA_OPTI_IMPLICIT_RR] = (double)NO;
	WDB[DATA_OPTI_GC_CALC] = (double)NO;

	if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
	  WDB[DATA_OPTI_MG_MODE] = (double)YES;	

	if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
	  Note(0, "Option 'set xscalc 2' ignored in optimization mode 2");

	WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)NO;
	
	break;
      }
    case 3:
      {
	WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)NO;
	WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)YES;
	WDB[DATA_OPTI_IMPLICIT_RR] = (double)YES;
	WDB[DATA_OPTI_GC_CALC] = (double)YES;

	if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
	  WDB[DATA_OPTI_MG_MODE] = (double)NO;	

	if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] < 0)
	  WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)YES;
      
	break;
      }
    case 4:
      {
	WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = (double)YES;
	WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = (double)YES;
	WDB[DATA_OPTI_IMPLICIT_RR] = (double)YES;
	WDB[DATA_OPTI_GC_CALC] = (double)YES;

	if ((long)RDB[DATA_OPTI_MG_MODE] < 0)
	  WDB[DATA_OPTI_MG_MODE] = (double)NO;	

	if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] < 0)
	  WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)YES;
      
	break;
      }
    default:
      Error(0, "Invalid optimization mode %ld", (long)RDB[DATA_OPTI_MODE]);
    }
  
  /* Switch to multi-group TMS mode if mg mode */

  if (((long)RDB[DATA_OPTI_MG_MODE] == YES) && 
      ((long)RDB[DATA_TMS_MODE] == TMS_MODE_CE))
    WDB[DATA_TMS_MODE] = (double)TMS_MODE_MG;

  /* Enforce delta-tracking if void cells are ignored */

  if ((long)RDB[DATA_IGNORE_VOID_CELLS] == YES)
    {
      WDB[DATA_DT_NTHRESH] = 0.0;
      WDB[DATA_DT_PTHRESH] = 0.0;
      WDB[DATA_OPT_USE_DT] = (double)YES;
    }

  /* Set unionization flag */

  if (((long)RDB[DATA_OPT_USE_DT] == YES) || 
      ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == YES) ||
      ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES))
    WDB[DATA_OPTI_UNIONIZE_GRID] = (double)YES;
  else
    WDB[DATA_OPTI_UNIONIZE_GRID] = (double)NO;

  /* Set group constant calculation on if universe is given and off if */
  /* set to null */

  if ((long)RDB[DATA_PTR_GCU0] > 0)
    WDB[DATA_OPTI_GC_CALC] = (double)YES;
  else if ((long)RDB[DATA_PTR_MORA0] > 0)
    WDB[DATA_OPTI_GC_CALC] = (double)YES;
  else if ((long)RDB[DATA_PTR_GCU0] < 0)
    WDB[DATA_OPTI_GC_CALC] = (double)NO;

  /* Set group constant calculation off if not neutron transport mode */

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == NO)
    WDB[DATA_OPTI_GC_CALC] = (double)NO;

  /* Use implicit reaction rates if group constants are generated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == YES)
    WDB[DATA_OPTI_IMPLICIT_RR] = (double)YES;

  /* Check that group constant calculation is used with fum */
  
  if (((long)RDB[DATA_OPTI_GC_CALC] == NO) && ((long)RDB[DATA_B1_CALC] == YES))
    Error(0, "Group constant calculation needed for B1 mode");

  /* Check on-the-fly mode */

  if ((long)RDB[DATA_TMS_MODE] != TMS_MODE_NONE)
    {
      /* Spectrum-collapse off */
      
      WDB[DATA_BU_SPECTRUM_COLLAPSE] = (double)NO;

      /* Grid thinning off */

      WDB[DATA_ERG_TOL] = 0.0;
    }

  /* Set delayed nubar flag if not set in input */

  if ((long)RDB[DATA_USE_DELNU] == -1)
    {
      /* Check simulation mode */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	WDB[DATA_USE_DELNU] = (double)YES;
      else
	WDB[DATA_USE_DELNU] = (double)NO;
    }

  /* Reset replay option if no reproducibility */

  if ((mpitasks > 1) && ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO))
    WDB[DATA_OPTI_REPLAY] = (double)NO;
  
  if (((long)RDB[DATA_OMP_MAX_THREADS] > 1) &&
      ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == NO))
    WDB[DATA_OPTI_REPLAY] = (double)NO;

  /* Disable grid thinning if problems */

  if (((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == NO) ||
      (((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == NO) &&
       ((long)WDB[DATA_OPTI_MG_MODE] == NO)))
    WDB[DATA_ERG_TOL] = 0.0;

  /* Disable shared results array if source biasing is in use (ei   */
  /* ole aavistustakaan miksi tämä tehdään, mutta se kaataa RES2-   */
  /* muistinvarauksen jos moodi on initdata.c:ssä YES, ja inputissa */
  /* on mesh plotteja. */

  /*
  if ((long)RDB[DATA_UFS_MODE] != UFS_MODE_NONE)
    WDB[DATA_OPTI_SHARED_RES2] = (double)NO;
  */
  /* Set explicit fission if criticality source simulation */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    WDB[DATA_OPT_IMPL_FISS] = (double)NO;

  /* Set neutron buffer factor */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    {
      /* Check value */

      if ((long)RDB[DATA_PART_NBUF_FACTOR] < 100)
	WDB[DATA_PART_NBUF_FACTOR] = 100.0;  
    }
  else if ((long)RDB[DATA_PART_NBUF_FACTOR] < 0.0)
    {
      /* Set default */

      WDB[DATA_PART_NBUF_FACTOR] = 5.0;

      /* Check if too low */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	if (RDB[DATA_CRIT_POP]/RDB[DATA_OMP_MAX_THREADS] < 1000.0)
	  WDB[DATA_PART_NBUF_FACTOR] = RDB[DATA_PART_NBUF_FACTOR]*
	    1000.0/RDB[DATA_CRIT_POP]*RDB[DATA_OMP_MAX_THREADS];
    }

  /* Set gamma buffer factor */

  if (RDB[DATA_PART_GBUF_FACTOR] < 0.0)
    {
      if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
	WDB[DATA_PART_GBUF_FACTOR] = 100.0;  
      else
	WDB[DATA_PART_GBUF_FACTOR] = 1.2;
    }

  /* Override micro-group batching for group constant generation */

  if ((long)RDB[DATA_MICRO_CALC_BATCH_SIZE] == -1)
    {
      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	WDB[DATA_MICRO_CALC_BATCH_SIZE] = RDB[DATA_CRIT_CYCLES];
      else
	WDB[DATA_MICRO_CALC_BATCH_SIZE] = RDB[DATA_SRC_BATCHES];
    }

  if ((long)RDB[DATA_BATCH_INTERVAL] > 1)
    WDB[DATA_MICRO_CALC_BATCH_SIZE] = RDB[DATA_BATCH_INTERVAL];

  /* Check implicit capture and multi-group mode */

  if (((long)RDB[DATA_OPTI_MG_MODE] == YES) && 
      ((long)RDB[DATA_OPT_IMPL_CAPT] == YES))
    Error(0, "Implicit capture does not work in this optimization mode");

  /* Wieland shift */
  
  if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
    {
      /* Sanity check of input */
      
      if ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
	Die(FUNCTION_NAME, 
	    "Wielandt shift can only be used with criticality source mode");

      /* Compute initial guess of shifted eigenvalue */

      if ((long)RDB[DATA_WIELANDT_MODE] == WIELANDT_MODE_FIX_K)
	{
	  WDB[DATA_CYCLE_KEFF] = 
	    (RDB[DATA_CYCLE_KEFF]*RDB[DATA_WIELANDT_KEFF])/
	    (RDB[DATA_WIELANDT_KEFF] - RDB[DATA_CYCLE_KEFF]);
	  
	  /* Check */

	  if (RDB[DATA_CYCLE_KEFF] < 0.3)
	    WDB[DATA_CYCLE_KEFF] = 0.3;
	}
    }

  /* Check analog (n,xn) with group constant generation */
  
  if (((long)RDB[DATA_OPTI_GC_CALC] == YES) &&
      ((long)RDB[DATA_OPT_IMPL_NXN] == NO))
    Error(0, "Analog (n,xn) does not work with group constant generation");

  /* Set weight cut-off if implicit capture is on */

  if ((long)RDB[DATA_OPT_IMPL_CAPT] == YES)
    if (RDB[DATA_OPT_ROULETTE_W0] < 0.0)
      WDB[DATA_OPT_ROULETTE_W0] = 0.001;

  /***************************************************************************/
}

/*****************************************************************************/
