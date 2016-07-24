#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : initdata.c                                     */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2016/04/04 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Inits values in main data block                              */
/*                                                                           */
/* Comments: - Tää vaatii vielä paljon työtä                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitData:"

/*****************************************************************************/

void InitData()
{
  long ptr;
  double val;
  char *path, *seed;
  time_t t0;
  clock_t cput0;

  /* Set error pointer */

  err = stdout;

  /* Check that system is 64 bit */

  if (sizeof(long) != 8)
    {
      /* Print error */

      fprintf(err, "\nSerpent 2 must be run in a 64-bit system. ");
      
      if (sizeof(long) == 4)
	fprintf(err, "Your system appears to be 32-bit.");

      fprintf(err, "\n\n");

      /* Abort */

      exit(-1);
    }

  /* Reset pointers */

  RDB = NULL;
  WDB = NULL;
  ACE = NULL;
  ASCII = NULL;
  PRIVA = NULL;
  BUF = NULL;
  RES1 = NULL;
  RES2 = NULL;
  SEED = NULL;
  SEED0 = NULL;

  /* Allocate memory for random number seed vectors */

  SEED = (unsigned long *)Mem(MEM_ALLOC, MAX_OMP_THREADS*RNG_SZ, 
			      sizeof(unsigned long));
  SEED0 = (unsigned long *)Mem(MEM_ALLOC, MAX_OMP_THREADS*RNG_SZ, 
			      sizeof(unsigned long));

  /* Set pointer to standard output */

  if (mpiid == 0)
    out = stdout;
  else
    {
      out = fopen("/dev/null", "w");
      return;
    }

  /* Set line-buffering for stdout */

  setlinebuf(out);

  /* Allocate memory for main data block */

  ReallocMem(DATA_ARRAY, DATA_FIXED_BLOCK_SIZE);

  /* Allow memory allocation */

  Mem(MEM_ALLOW);

  /* Allocate data at the beginning of RES1 array to simplify error testing */
  
  ReallocMem(RES1_ARRAY, DATA_FIXED_BLOCK_SIZE);

  /* Allocate memory for ASCII data block to simplify error testing */

  ASCII = (char *)Mem(MEM_ALLOC, (long)(VALID_PTR + 1), sizeof(char)); 
  WDB[DATA_ASCII_DATA_SIZE] = (double)(VALID_PTR + 1);

  /* Reset status of private memory blocks */

  WDB[DATA_PRIVA_MEM_READY] = (double)NO;

  /* Init list pointers */

  WDB[DATA_PTR_C0] = NULLPTR;
  WDB[DATA_PTR_S0] = NULLPTR;
  WDB[DATA_PTR_M0] = NULLPTR;
  WDB[DATA_PTR_L0] = NULLPTR;
  WDB[DATA_PTR_T0] = NULLPTR;
  WDB[DATA_PTR_NST0] = NULLPTR;
  WDB[DATA_PTR_GPL0] = NULLPTR;
  WDB[DATA_PTR_TR0] = NULLPTR;
  WDB[DATA_PTR_PB0] = NULLPTR;
  WDB[DATA_PTR_IFC0] = NULLPTR;
  WDB[DATA_PTR_DIV0] = NULLPTR;

  /* Init file names */

  WDB[DATA_PTR_ACEDATA_FNAME_LIST] = NULLPTR;
  WDB[DATA_PTR_DECDATA_FNAME_LIST] = NULLPTR;
  WDB[DATA_PTR_NFYDATA_FNAME_LIST] = NULLPTR;
  WDB[DATA_PTR_SFYDATA_FNAME_LIST] = NULLPTR;

 /* Maximum fraction of system memory to use */
  /* Try to get fraction from environmental variable */

  if ((path = getenv("SERPENT_MEM_FRAC")) != NULL)
    if (sscanf(path, "%lf", &val) == 1)
      if ((val > 0.0) && (val < 1.0))
	WDB[DATA_CPU_MEM_FRAC] = val;

  /* Set to default value if not set in environmental variable */

  if (RDB[DATA_CPU_MEM_FRAC] == 0.0)
    WDB[DATA_CPU_MEM_FRAC] = 0.8;

  /* Few-group constant generation */

  WDB[DATA_ERG_FG_PTR_PREDEF] = (double)PutText("default2");

  /* Thinning tolerance and energy boundaries (tarkista noi fotonirajat) */

  WDB[DATA_ERG_TOL] = -1.0;

  WDB[DATA_NEUTRON_EMIN] = 1E-11;
  WDB[DATA_NEUTRON_EMAX] = 20.0;

  WDB[DATA_PHOTON_EMIN] = 1E-3;
  WDB[DATA_PHOTON_EMAX] = 100.0;

  /* Set boundary conditions and albedos (NOTE: Toi STOP_AT_BOUNDARY */
  /* käännettiin pois juuri ennen 2.1.25 jakelua) */

  WDB[DATA_STOP_AT_BOUNDARY] = (double)NO;
  WDB[DATA_GEOM_BC0] = (double)BC_BLACK;
  WDB[DATA_GEOM_BC1] = (double)BC_BLACK;
  WDB[DATA_GEOM_BC2] = (double)BC_BLACK;
  WDB[DATA_GEOM_BC3] = (double)BC_BLACK;
  WDB[DATA_GEOM_ALBEDO1] = 1.0;
  WDB[DATA_GEOM_ALBEDO2] = 1.0;
  WDB[DATA_GEOM_ALBEDO3] = 1.0;

  /* K-eff iteration */

  WDB[DATA_ITER_MODE] = (double)ITER_MODE_NONE;
  WDB[DATA_ITER_KEFF] = 1.0;
  WDB[DATA_ITER_VAL] = -1.0;
  WDB[DATA_ITER_NCYC] = 50.0;
  WDB[DATA_ITER_FIX] = (double)NO;

  WDB[DATA_ITER_ALB_F1] = 1.0;
  WDB[DATA_ITER_ALB_F2] = 1.0;
  WDB[DATA_ITER_ALB_F3] = 1.0;

  /* Set initial date */

  WDB[DATA_PTR_DATE] = (double)PutText(TimeStamp());

  /* Get initial cpu time */

  cput0 = clock();
  WDB[DATA_CPU_T0] = (double)cput0;

  /* Reset calculation modes */

  WDB[DATA_NEUTRON_TRANSPORT_MODE] = (double)NO;
  WDB[DATA_PHOTON_TRANSPORT_MODE] = (double)NO;
  WDB[DATA_BURNUP_CALCULATION_MODE] = (double)NO;
  WDB[DATA_VOLUME_CALCULATION_MODE] = (double)NO;
  WDB[DATA_PARTICLE_DISPERSER_MODE] = (double)NO;
  WDB[DATA_PARTICLE_REDEPLETE_MODE] = (double)NO;
  WDB[DATA_QUICK_PLOT_MODE] = (double)NO;
  WDB[DATA_PHOTON_PRODUCTION] = (double)NO;

  /* Normalisation */

  WDB[DATA_PTR_NORM] = NULLPTR;
  WDB[DATA_NORM_U235_FISSE] = U235_FISSE;

  /* Replay mode */

  WDB[DATA_OPTI_REPLAY] = (double)NO;

  /* Set default random number seed */

  time(&t0);
  parent_seed = (unsigned long)t0;

  /* Override with environment variable */

  if ((seed = getenv("SERPENT_RNG_SEED")) != NULL)
    {
      parent_seed = (unsigned long)atoi(seed);
      WDB[DATA_OPTI_REPLAY] = (double)YES;
    }

  /* XS data plotter parameters */

  WDB[DATA_XSPLOT_NE] = -1.0;
  WDB[DATA_XSPLOT_EMIN] = -1.0;
  WDB[DATA_XSPLOT_EMAX] = -1.0;

  /* Print compositions */

  WDB[DATA_BURN_PRINT_COMP] = (double)NO;
  WDB[DATA_BURN_PRINT_COMP_LIM] = 1.0;
  WDB[DATA_BURN_MAT_OUTPUT] = (double)BURN_OUT_MAT_PARENT;
  WDB[DATA_BURN_PRINT_INTERMEDIATE] = (double)YES;

  /* Decay only mode */

  WDB[DATA_BURN_DECAY_CALC] = (double)NO;

  /* Print histories */
  
  WDB[DATA_OPTI_PRINT_HIS] = (double)NO;

  /* Depletion cut-offs */

  WDB[DATA_DEP_HALF_LIFE_CUTOFF] = INFTY;
  WDB[DATA_DEP_TTA_CUTOFF]       = 1E-28;

  /* Time cut-off for transport simulation */

  WDB[DATA_TIME_CUT_TMIN] = 0.0;
  WDB[DATA_TIME_CUT_TMAX] = INFTY;

  /* Reset generation cut-off */

  WDB[DATA_GEN_CUT] = (double)MAX_GENERATIONS;
  WDB[DATA_MAX_PROMPT_CHAIN_LENGTH] = -1.0;
  
  /* Fission source entropy */

  WDB[DATA_OPTI_ENTROPY_CALC] = (double)NO;

  WDB[DATA_ENTROPY_NX] = 5.0;
  WDB[DATA_ENTROPY_NY] = 5.0;
  WDB[DATA_ENTROPY_NZ] = 5.0;

  WDB[DATA_ENTROPY_XMIN] = -INFTY;
  WDB[DATA_ENTROPY_XMAX] =  INFTY;

  WDB[DATA_ENTROPY_YMIN] = -INFTY;
  WDB[DATA_ENTROPY_YMAX] =  INFTY;

  WDB[DATA_ENTROPY_ZMIN] = -INFTY;
  WDB[DATA_ENTROPY_ZMAX] =  INFTY;

  /* Source point animation */

  WDB[DATA_SOURCE_PT_ANIM] = (double)NO;
  WDB[DATA_SOURCE_PT_ANIM_F] = 1.0;
  WDB[DATA_SOURCE_PT_ANIM_PALETTE] = (double)PALETTE_HOT; 

  /* Microscopic partial total cross section limit for including nuclide */
  /* in reaction lists. */

  WDB[DATA_MIN_TOTXS] = 1E-5;

  /* Use unresolved resonance probability tables (turned off by default) */

  WDB[DATA_USE_URES] = (double)NO;
  WDB[DATA_URES_PTR_USE_LIST] = NULLPTR;

  /* Plotter */

  WDB[DATA_STOP_AFTER_PLOT] = (double)STOP_AFTER_PLOT_NONE;
  WDB[DATA_PLOTTER_MODE] = (double)NO;

  /* Optimization and memory/data options */

  WDB[DATA_OPTI_MODE] = 4.0;
  WDB[DATA_OPTI_UNIONIZE_GRID] = -1.0;
  WDB[DATA_OPTI_RECONSTRUCT_MICROXS] = -1.0;
  WDB[DATA_OPTI_RECONSTRUCT_MACROXS] = -1.0;
  WDB[DATA_OPTI_INCLUDE_SPECIALS] = (double)NO;

  /* Delta-tracking */

  WDB[DATA_OPT_USE_DT] = (double)YES;
  WDB[DATA_DT_NTHRESH] = 0.1;
  WDB[DATA_DT_PTHRESH] = 0.1;

  /* Implicit Monte Carlo */

  WDB[DATA_OPT_IMPL_FISS] = (double)NO;
  WDB[DATA_OPT_IMPL_CAPT] = (double)NO;
  WDB[DATA_OPT_IMPL_NXN] = (double)YES;

  /* Set roulette off by default */

  WDB[DATA_OPT_ROULETTE_W0] = -1.0;
  WDB[DATA_OPT_ROULETTE_P0] = 0.5;

  /* Run-time variables */

  WDB[DATA_CYCLE_KEFF] = 1.0;

  /* Simulation mode and batching interval */

  WDB[DATA_SIMULATION_MODE] = -1.0;
  WDB[DATA_BATCH_INTERVAL] = 1.0;

  /* OpenMP stuff */

  WDB[DATA_OMP_MAX_THREADS] = 1.0;

  /* Allocate empty memory ACE array to deal with zero pointers */

  WDB[DATA_PTR_ACE0] = NULLPTR;
  ReallocMem(ACE_ARRAY, DATA_FIXED_BLOCK_SIZE);

  WDB[DATA_PTR_ACE_NFY_DATA] = NULLPTR;
  WDB[DATA_PTR_ACE_SFY_DATA] = NULLPTR;

  /* Root universe */

  WDB[DATA_PTR_ROOT_UNIVERSE] = (double)PutText("0");

  /* History list size */

  WDB[DATA_HIST_LIST_SIZE] = -50.0;

  /* Event bank size and maximum number of generations */

  WDB[DATA_EVENT_BANK_SZ] = 100.0;
  WDB[DATA_EVENT_MAX_GEN] = 20.0;
  WDB[DATA_EVENT_RECORD_FLAGS] = 0.0;

  /* Monte Carlo volume calculation */

  WDB[DATA_VOLUME_MC_NMAX] = -1.0;
  WDB[DATA_VOLUME_MC_TMAX] = -1.0;
  WDB[DATA_VOLUME_MC_EMAX] = -1.0;

  /* Default confidence levels for TMP majorant generation */
  
#ifndef TRADMAJO  

  /* Revisited majorant */

  WDB[DATA_QPARAM_TMS] = 2.0E-5; 
  WDB[DATA_QPARAM_DBRC] = 2.0E-5; 

#else
  
  /* "f-parameter" of traditional majorant */

  WDB[DATA_QPARAM_TMS] = 3.0; 
  WDB[DATA_QPARAM_DBRC] = 3.0; 

#endif 

  /* DBRC */

  WDB[DATA_USE_DBRC] = (double)NO;
  WDB[DATA_PTR_DBRC] = NULLPTR;
  WDB[DATA_DBRC_EMIN] = 0.4E-6;
  WDB[DATA_DBRC_EMAX] = 210E-6;

  /* Energy cut-off for transmutation reactions */

  WDB[DATA_BURN_ENECUT] = 10.0;

  /* Normalization */

  WDB[DATA_NORM_BURN] = BURN_NORM_ALL;

  /* Burnup mode */

  WDB[DATA_BURN_BUMODE] = (double)BUMODE_CRAM;
  WDB[DATA_BU_SPECTRUM_COLLAPSE] = -1.0;
  
  /* CRAM:n asteluku */

  WDB[DATA_BURN_CRAM_K] = 14.0;

  /* Predictor-corrector calculation */

  WDB[DATA_BURN_PRED_TYPE] = (double)PRED_TYPE_CONSTANT;
  WDB[DATA_BURN_PRED_NSS] = 1.0;
  WDB[DATA_BURN_CORR_TYPE] = (double)CORR_TYPE_LINEAR;
  WDB[DATA_BURN_CORR_NSS] = 1.0;
  WDB[DATA_BURN_SIE] = (double)NO;

  /* Corrector iteration */
  
  WDB[DATA_BURN_CI_TYPE] = (double)CI_TYPE_OUTER;
  WDB[DATA_BURN_CI_MAXI] = 1.0;
  WDB[DATA_BURN_CI_LAST] = (double)NO;
  
  WDB[DATA_BURN_CI_NBATCH] = -1.0;
  WDB[DATA_BURN_CI_CYCLES] = -1.0;
  WDB[DATA_BURN_CI_SKIP]   = -1.0;

  /* Core power distribution */

  WDB[DATA_CORE_PDE_DEPTH] = -1.0;

  /* Uniform fission source method */

  WDB[DATA_UFS_PTR_SRC_MESH] = NULLPTR;
  WDB[DATA_UFS_ORDER] = 1.0;
  WDB[DATA_UFS_MIN] = 0.001;
  WDB[DATA_UFS_MAX] = 10.0;

  /* Energy grid for coarse multi-group xs */

  WDB[DATA_COARSE_MG_NE] = 4000.0;
  WDB[DATA_OPTI_MG_MODE] = -1.0;

  /* Reset minimum and maximum energies of cross section data */

  WDB[DATA_NEUTRON_XS_EMIN] = INFTY;
  WDB[DATA_NEUTRON_XS_EMAX] = -INFTY;

  WDB[DATA_PHOTON_XS_EMIN] = INFTY;
  WDB[DATA_PHOTON_XS_EMAX] = -INFTY;

  /* Delayed nubar flag */

  WDB[DATA_USE_DELNU] = -1.0;

  /* Doppler broadening mode */

  WDB[DATA_TMS_MODE] = TMS_MODE_NONE;
  WDB[DATA_USE_DENSITY_FACTOR] = (double)NO;
  WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)NO;

  /* Ures energy boundaries */

  WDB[DATA_URES_EMIN] = INFTY;
  WDB[DATA_URES_EMAX] = -INFTY;

  /* Fundamental mode calculation */

  WDB[DATA_B1_CALC] = (double)NO;
  WDB[DATA_B1_BURNUP_CORR] = (double)NO;
  WDB[DATA_B1_ERR_LIMIT] = 1E-5;

  /* Default micro-group structure */

  WDB[DATA_MICRO_PTR_EGRID] = (double)PutText("defaultmg");

  /* Batching for micro-group calculation (-1 = set equal to */
  /* number of cycles) */

  WDB[DATA_MICRO_CALC_BATCH_SIZE] = 20;

  /* Shared scoring buffer */

  WDB[DATA_OPTI_SHARED_BUF] = (double)NO;

  /* Shared RES2 array */

  WDB[DATA_OPTI_SHARED_RES2] = (double)YES;

  /* Reproducibility in OpenMP and MPI modes */

  WDB[DATA_OPTI_OMP_REPRODUCIBILITY] = (double)YES;
  WDB[DATA_OPTI_MPI_REPRODUCIBILITY] = (double)NO;

  /* Include scattering production in removal xs */

  WDB[DATA_GC_REMXS_MULT] = (double)YES;

  /* Solution relaxation factor for coupled calculation */

  WDB[DATA_SOL_REL_FACT] = 1.0;

  /* Solution relaxation maximum population for coupled calculation */

  WDB[DATA_SOL_REL_MAX_POP] = INFTY/1e6;

  /* Analog reaction rate calculation */

  WDB[DATA_ANA_RR_NCALC] = (double)ARR_MODE_NONE;
  WDB[DATA_ANA_RR_PCALC] = (double)ARR_MODE_NONE;

  /* Terminate on die */

  WDB[DATA_TERMINATE_ON_DIE] = (double)YES;

  /* Recorded tracks for plots */

  WDB[DATA_TRACK_PLOTTER_HIS] = -1.0;

  /* Track plots */

  WDB[DATA_TRACK_PLOT_TMIN] = 0.0;
  WDB[DATA_TRACK_PLOT_TMAX] = INFTY;
  WDB[DATA_TRACK_PLOT_FRAMES] = 1.0;
  WDB[DATA_TRACK_PLOT_NHIS] = 100.0;

  /* Number of parallel eigenvalue calculations */

  WDB[DATA_N_POP_EIG] = 1.0;

  /* UFS mode */

  WDB[DATA_UFS_MODE] = (double)UFS_MODE_NONE;

  /* Default cross section data library */

  if ((path = getenv("SERPENT_ACELIB")) != NULL)
    {
      ptr = ReallocMem(DATA_ARRAY, 2);
      WDB[DATA_PTR_ACEDATA_FNAME_LIST] = (double)ptr;
      WDB[ptr++] = (double)PutText(path);
      WDB[ptr] = NULLPTR;
    }

  /* Reset neutron, gamma and precursor buffer factors */

  WDB[DATA_PART_NBUF_FACTOR] = -1.0;  
  WDB[DATA_PART_GBUF_FACTOR] = -1.0;
  WDB[DATA_PART_PBUF_FACTOR] = -1.0;

  /* Reaction sampling */

  WDB[DATA_NPHYS_SAMPLE_FISS] = (double)YES;
  WDB[DATA_NPHYS_SAMPLE_CAPT] = (double)YES;
  WDB[DATA_NPHYS_SAMPLE_SCATT] = (double)YES;

  /* Number of time bins */

  WDB[DATA_DYN_NB] = 1.0;

  /* Alpha eigenvalue */

  WDB[DATA_ALPHA_EIG] = 0.0; /* 1.431E+04 = EPR-laskun kriittinen alpha */

  /* Minimum xs for CFE */

  WDB[DATA_CFE_N_MIN_L] = 100.0;
  WDB[DATA_CFE_N_MIN_T] = -1.0;
  WDB[DATA_CFE_G_MIN_L] = 100.0;
  WDB[DATA_CFE_G_MIN_T] = -1.0;

  /* Energy boundary for cache-optimized xs block */
  
  WDB[DATA_CACHE_OPTI_EMAX] = 1E-5;

  /* Ures infinite dilution cut-off */

  WDB[DATA_URES_DILU_CUT] = 1E-9;

  /* Poison calculation */

  WDB[DATA_OPTI_POISON_CALC] = (double)NO;

  /* Equilibrium poison calculation */

  WDB[DATA_XENON_EQUILIBRIUM_MODE] = -1.0;
  WDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] = -1.0;

  /* MPI batch size */

  WDB[DATA_OPTI_MPI_BATCH_SIZE] = 10000.0;

  /* Actinide limits for burnup calculation */

  WDB[DATA_BU_ACT_MIN_Z] = 90.0;
  WDB[DATA_BU_ACT_MAX_Z] = 96.0;

  /* Number of progenies for beta-eff and prompt lifetime calculation */

  WDB[DATA_IFP_OPT_PRINT_ALL] = (double)NO;
  WDB[DATA_IFP_CHAIN_LENGTH] = 15.0;
  WDB[DATA_PERT_VAR_A] = 1E-8;
  WDB[DATA_PERT_VAR_C] = 1E-6;
  WDB[DATA_PERT_N_BATCH] = 9.0;

  /* ICM stuff */

  WDB[DATA_ICM_CALC] = (double)NO;
  WDB[DATA_ICM_NSEG] = -1.0;
  WDB[DATA_ICM_NSUB] = 1.0;
  WDB[DATA_ICM_NMU0] = 1.0;
  WDB[DATA_ICM_NMU1] = 1.0;
  WDB[DATA_ICM_NMU2] = 1.0;

  /* Group constant generation at multiple levels */

  WDB[DATA_MULTI_LEVEL_GCU] = (double)NO;

  /* Void cells */

  WDB[DATA_IGNORE_VOID_CELLS] = (double)NO;

  /* Print interval */

  WDB[DATA_PRINT_INTERVAL] = 50.0;

  /* Write/write restart file */

  WDB[DATA_WRITE_RESTART_FILE] = (double)NO;
  WDB[DATA_READ_RESTART_FILE] = (double)NO;

  /* Statistical tests on group constants */

  WDB[DATA_GC_STAT_TESTS] = (double)NO;
  WDB[DATA_RUN_STAT_TESTS] = (double)NO;

  /* STL geometry stuff */
  
  WDB[DATA_STL_TEMP_ARRAY_SIZE] = 100.0;
  WDB[DATA_STL_TEST_N_PTS] = 0.0;
  WDB[DATA_STL_TEST_N_DIR] = 0.0;
  WDB[DATA_STL_GEOM_TEST_MODE] = (double)NO;

  /* NOTE: Jos tää on suurempi niin se voi aiheuttaa infinite loopin */
  /* ST-moodissa (esim pupu testikeissi). Pienempi arvo puolestaan   */
  /* voi aiheuttaa geometriaerroreita kun säteet osuu liian lähelle  */
  /* kolmion reunaa. */

  WDB[DATA_STL_FACET_EXD] = 1E-4;

  /* NOTE: Tää estää ST:n käytön pisteistä jotka lähtee facetien */
  /* bounding boxien sisältä. */

  WDB[DATA_STL_ENFORCE_DT] = (double)YES;

  /* Reset parameters for Wielandt method */

  WDB[DATA_WIELANDT_MODE] = (double)WIELANDT_MODE_NONE;
  WDB[DATA_WIELANDT_KEFF] = INFTY;
  WDB[DATA_WIELANDT_KP] = 0.0;
  WDB[DATA_WIELANDT_P] = 0.0;

  /* Reset coefficient calculation index and error mode */

  WDB[DATA_COEF_CALC_IDX] = -1.0;
  WDB[DATA_COEF_CALC_INCLUDE_ERRORS] = (double)NO;

  /* Reset flag for more coefficient calculations */

  WDB[DATA_MORE_COEF_CALC] = (double)NO;
  WDB[DATA_COEF_CALC_SPECIAL_MODE] = 0.0;

  /* Homogeneous flux solution */

  WDB[DATA_HOMOFLUX_SOLVER] = 1.0;
  WDB[DATA_ADF_TRAPZ_PT] = 100.0;

  /* Klein-Nishina threshold energy */

  WDB[DATA_PHOTON_EKN] = INFTY;

  /* Use Doppler broadening of Compton photons */

  WDB[DATA_PHOTON_USE_DOPPLER] = (double)YES;

  /* TODO: Statistics: Total number of bremsstrahlung photons created */
  /* WDB[DATA_PHOTON_BREM_TOT] = 0; */

  /* Use TTB-approximation */

  WDB[DATA_PHOTON_USE_TTB] = (double)YES;
  
  /* Separate TTB data for positrons */

  WDB[DATA_PHOTON_TTBPM] = (double)YES;

  /* Local energy conservation of TTB */

  WDB[DATA_PHOTON_TTBEC] = (double)NO;

  /* Detailed angular distribution for Compton electrons */

  WDB[DATA_PHOTON_COMP_EANG] = (double)NO;

  /* Photon data file names */

  WDB[DATA_PHOTON_COH_FNAME] = PutText("cohff.dat");
  WDB[DATA_PHOTON_INCOH_FNAME] = PutText("incohsf.dat");
  WDB[DATA_PHOTON_RELAX_FNAME] = PutText("relax.dat");
  WDB[DATA_PHOTON_PESS_FNAME] = PutText("xspess.dat");
  WDB[DATA_PHOTON_PETOT_FNAME] = PutText("xspetot.dat");
  WDB[DATA_PHOTON_CP_FNAME] = PutText("ComptonProfiles.dat");
  WDB[DATA_PHOTON_ELSP_FNAME] = PutText("el_stopping_power.dat");
  WDB[DATA_PHOTON_ELBR_FNAME] = PutText("pdebr.dat");

  /* Decay source flag */

  WDB[DATA_USE_DECAY_SRC] = (double)NO;

  /* Maximum tracking loop */

  WDB[DATA_NEUTRON_MAX_TRACK_LOOP] = 1000000.0;
  WDB[DATA_PHOTON_MAX_TRACK_LOOP] = 10000.0;

  /* Error flags */

  WDB[DATA_NEUTRON_MAX_TRACK_LOOP_ERR] = (double)YES;
  WDB[DATA_PHOTON_MAX_TRACK_LOOP_ERR] = (double)NO;

  /* Weight windows */
  
  WDB[DATA_USE_WEIGHT_WINDOWS] = (double)NO;

  WDB[DATA_WWD_LOWER_BOUND] = 0.5;
  WDB[DATA_WWD_UPPER_BOUND] = 2.0;

  /* Heat and gamma production cross sections */

  WDB[DATA_INCLUDE_HEAT_PROD_XS] = (double)NO;
  WDB[DATA_INCLUDE_PHOT_PROD_XS] = (double)NO;

  /* Factor for reading f*RDB[SRC_POP] precursors points for initial source */
  /* in transient simulations */

  WDB[DATA_PREC_SRC_FACT] = 10.0;

  /* Threshold for storing precursors generated by implicit estimator     */
  /* during transient simulation, negative value is based on prec. weight */
  /* Positive based on prec. emission during typical time interval */

  WDB[DATA_PREC_STORE_TRESH] = 1.0;

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
