#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : transportcycle.c                               */
/*                                                                           */
/* Created:       2011/05/23 (JLe)                                           */
/* Last modified: 2016/04/05 (VVa)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Prepares and runs the main transport cycle                   */
/*                                                                           */
/* Comments: - Tää on ihan kauhea sekasotku!!!                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TransportCycle:"

/*****************************************************************************/

void TransportCycle()
{
  long nb, nb0, maxb, skip, nn, nt, maxt, id, idx, ptr, tme;
  long nsrc, tosimulate;
  double t0, c0;

  /***************************************************************************/

  /***** Main transport cycle ************************************************/

  /* Reset total and active timer */

  ResetTimer(TIMER_TRANSPORT);
  ResetTimer(TIMER_TRANSPORT_ACTIVE);

  /* Start transport timer */
      
  StartTimer(TIMER_TRANSPORT);
  StartTimer(TIMER_TRANSPORT_TOTAL);

  /* Reset completed flag and sort counter */

  WDB[DATA_SIMULATION_COMPLETED] = (double)NO;
  WDB[DATA_SORT_COUNT] = 1.0;

  /* Check mode */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
    {
      /***********************************************************************/

      /***** External source simulation **************************************/

      /* Reset skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* Set cycle-wise batch size */

      WDB[DATA_CYCLE_BATCH_SIZE] = RDB[DATA_SRC_POP];
      
      /* Set number of source points in batch */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_SRC_POP];

      /* Get number of bins */

      maxt = (long)RDB[DATA_DYN_NB];
      CheckValue(FUNCTION_NAME, "maxt", "", maxt, 1, 1000000000);

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);
      
      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);    
     
      /* Start active transport timer */

      StartTimer(TIMER_TRANSPORT_ACTIVE);

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      /* Reset batch counters */
	      
      WDB[DATA_BATCH_COUNT] = 0.0;
      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;

      /* Check batching interval */

      if ((long)RDB[DATA_SRC_BATCHES] % (long)RDB[DATA_BATCH_INTERVAL])
	Error(0, 
	    "Total number of batches %ld is not a multiple of interval %ld",
	    (long)RDB[DATA_SRC_BATCHES], (long)RDB[DATA_BATCH_INTERVAL]);

      /* Start simulation */

      fprintf(out, "Starting external source simulation...\n\n");

      /* Loop over batches */
      
      for (nb = 0; nb < (long)RDB[DATA_SRC_BATCHES]; nb++)
	{
	  /*******************************************************************/
	  
	  /***** Loop over source batches ************************************/

	  /* Reset time bin index */

	  WDB[DATA_DYN_TB] = 0.0;

	  /* Set time cut-off for first cycle */
	  
	  WDB[DATA_TIME_CUT_TMIN] = RDB[tme];
	  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + 1];
	  
	  /* Stop tracks at outer boundary if time cut-off is used */

	  if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
	    WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

	  /* Reset cycle k-eff and put starting weight */

	  WDB[DATA_CYCLE_KEFF] = 1.0;
	  WDB[DATA_DYN_WGT0] = RDB[DATA_SRC_POP];

	  /* Put cycle index */

	  WDB[DATA_CYCLE_IDX] = (double)nb;

	  /* Get beginning time */

	  t0 = TimerVal(TIMER_TRANSPORT);
	  c0 = TimerCPUVal(TIMER_TRANSPORT);

	  /* Start parallel timer */

	  StartTimer(TIMER_OMP_PARA);

	  /* Parallel loop over histories */

#ifdef OPEN_MP
#pragma omp parallel private(id, idx, nn) 
#endif
	  {
	    /* Get Open MP thread id */
		
	    id = OMP_THREAD_NUM;
		
#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif	  
	    /* Loop over source neutrons */
	  
	    for (nn = 0; nn < (long)RDB[DATA_SRC_POP]; nn++)
	      {
		/* Calculate particle index */

		idx = (long)RDB[DATA_NHIST_TOT];
		idx = idx + (long)(nb*RDB[DATA_SRC_POP]) + nn;

		/* Sample source point */

		SampleSrcPoint(id, nn, idx);

		/* Track */

		Tracking(id);
	      }
	  }

	  /* Move collected pulse data to statistics */

	  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
	    PulseDet(-1, -1, -1.0, 0.0, 0.0, 0.0, -1.0, -1);
	
	  /* Stop parallel timer */

	  StopTimer(TIMER_OMP_PARA);

	  /* Collect data from interval */

	  CollectDynData();
	  
	  /* Reset normalization coefficients (pitää kutsua tässä) */
	  
	  WDB[DATA_NORM_COEF_N] = -1.0;
	  WDB[DATA_NORM_COEF_G] = -1.0;

	  /* Update batch count (JLe: siirretty tuolta alempaa tänne */
	  /* 1.9.2015 / 2.1.24) */
	      
	  WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

	  /* Check for remaining time intervals */

	  if (maxt > 1)
	    {
	      /* Fix normalization to first step */
	      
	      if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
		NormCoef(PARTICLE_TYPE_NEUTRON);
	      if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
		NormCoef(PARTICLE_TYPE_GAMMA);
	      
	      /* Re-open buffer for writing */

	      WDB[DATA_BUF_REDUCED] = (double)NO;
	      
	      /* Loop over remaining time intervals */
	      
	      for (nt = 1; nt < maxt; nt++)
		{
		  
		  /* Re-open buffer for writing */
		  
		  WDB[DATA_BUF_REDUCED] = (double)NO;

		  /* Set time bin index */

		  WDB[DATA_DYN_TB] = (double)nt;

		  /* Set new time cut-off */

		  WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
		  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];
		  
		  /* Normalize source */
		  
		  if (NormalizeDynSrc() < 0)
		    break;
		  
		  /* Start parallel timer */
		  
		  StartTimer(TIMER_OMP_PARA);
		  
		  /* Loop until source is empty */
		  
#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
		  {
		    /* Get Open MP thread id */
		    
		    id = OMP_THREAD_NUM;
		    
		    /* Loop over source */
		    
		    while(FromSrc(id) > VALID_PTR)
		      Tracking(id);      
		  }

		  /* Collect data from interval */

		  CollectDynData();
		  
		  /* Stop parallel timer */
		  
		  StopTimer(TIMER_OMP_PARA);
		}
	    }
	      
	  /* Get end time */

	  t0 = TimerVal(TIMER_TRANSPORT) - t0;
	  c0 = TimerCPUVal(TIMER_TRANSPORT) - c0;

	  /* Score time */

	  ptr = (long)RDB[RES_CYCLE_RUNTIME];
	  AddStat(t0, ptr, 0); 
	  AddStat(c0, ptr, 1); 
  
	  /* Add to micro batch counter */

	  WDB[DATA_MICRO_CALC_BATCH_COUNT] = 
	    RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;
	  
	  /* Check batch interval */

	  if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
	    {
	      /* Collect and clear buffered results */

	      CollectResults();
	      CollectDet();
	      PoisonEq();
	      CalcMicroGroupXS();
	      ClearBuf();

	      /* Reset batch counter */
	      
	      WDB[DATA_BATCH_COUNT] = 0.0;
	    }

	  /* Flush bank */

	  FlushBank();

	  /* Stop cycle-wise transport timer */

	  StopTimer(TIMER_TRANSPORT_CYCLE);
	  
	  if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
	    t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
	      TimerVal(TIMER_TRANSPORT_CYCLE);
	  else
	    t0 = 0.0;

	  /* CPU usage */

	  ptr = (long)RDB[RES_CPU_USAGE];
	  AddStat(t0, ptr, 0); 

	  /* Print cycle-wise output */

	  PrintCycleOutput();

	  /* Reset and restart cycle-wise transport timer */

	  ResetTimer(TIMER_TRANSPORT_CYCLE);
	  StartTimer(TIMER_TRANSPORT_CYCLE);

	  /* Sort lists */

	  SortAll();

	  /* Print results */
	  
	  if (!((nb + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
	    {
	      MatlabOutput();
	      DetectorOutput();
	      MeshPlotter();
	      PrintCoreDistr();
	      PrintHistoryOutput();
	      PrintPBData();
	      PrintInterfaceOutput();
	      PrintFinix();
	      RROutput();
	      /*
	      GeometryPlotter(NO);
	      */
	      FissMtxOutput();
	      MORAOutput();
	      WriteICMData();
	    }

	  /******************************************************************/
	}

      /* Update number of neutron histories run */

      WDB[DATA_NHIST_TOT] = WDB[DATA_NHIST_TOT] +
	RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN)
    {
      /***********************************************************************/

      /***** Time dependent simulation with delayed neutron ******************/
      /* NB: This will be combined with SIMULATION_MODE_DYN                  */
      /* This is only set on, if a precursor detector is explicitly set up in*/
      /* input and the mode would otherwise be SIMULATION_MODE_SRC           */

      /* Reset skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* Set cycle-wise batch size */

      WDB[DATA_CYCLE_BATCH_SIZE] = RDB[DATA_SRC_POP];
      
      /* Set number of source points in batch */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_SRC_POP];

      /* Get number of bins */

      maxt = (long)RDB[DATA_DYN_NB];
      CheckValue(FUNCTION_NAME, "maxt", "", maxt, 1, 1000000000);

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);
      
      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);    
     
      /* Start active transport timer */

      StartTimer(TIMER_TRANSPORT_ACTIVE);

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      /* Start simulation */

      fprintf(out, "Starting time dependent simulation with delayed neutrons...\n\n");
      
      /* Loop over batches */
      
      for (nb = 0; nb < (long)RDB[DATA_SRC_BATCHES]; nb++)
	{

	  /* Reset time bin index */

	  WDB[DATA_DYN_TB] = 0.0;

	  /* Set time cut-off for first cycle */
	  
	  WDB[DATA_TIME_CUT_TMIN] = RDB[tme];
	  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + 1];

	  /* Initialize values in precursor detectors */
	  /* This initializes the values in buffers for the first */
	  /* time bin, currently this has to be done separately for */
	  /* each batch */

	  InitPrecDet();

	  /* Set simulation normalization and calculate weights */
	  /* of live neutrons and precursors */

	  NormalizePrecDet();

	  /* Do population control for initial source */

	  PrecursorPopControl();

	  /* Handle decay of precursor during interval */
	  /* and calculate weight to emit */

	  DecayMeshPrecDet();
	  
	  /* Stop tracks at outer boundary if time cut-off is used */

	  if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
	    WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

	  /* Reset cycle k-eff and put starting weight */

	  WDB[DATA_CYCLE_KEFF] = 1.0;
	  WDB[DATA_DYN_WGT0] = RDB[DATA_SRC_POP];

	  /* Set batch counters */

	  WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_INTERVAL];
	  WDB[DATA_MICRO_CALC_BATCH_COUNT] = 
	    RDB[DATA_MICRO_CALC_BATCH_SIZE];

	  /* Put cycle index */

	  WDB[DATA_CYCLE_IDX] = (double)nb;

	  /* Get beginning time */

	  t0 = TimerVal(TIMER_TRANSPORT);
	  c0 = TimerCPUVal(TIMER_TRANSPORT);

	  /* Sample delayed neutrons for the interval */

	  SampleDelnu();

	  /* Handle decay of pointwise precursors over interval */

	  DecayPointPrecDet();

	  /* Re-open buffer for writing */
	  /* Was closed in sampledelnu  */

	  WDB[DATA_BUF_REDUCED] = (double)NO;

	  /* Start parallel timer */

	  StartTimer(TIMER_OMP_PARA);

	  /* Get number of source neutrons */

	  nsrc = (long)RDB[DATA_SRC_POP];

	  /* Try to override with number of live neutrons */

	  if ((ptr = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
	    nsrc = (long)RDB[ptr + PRECDET_N_LIVE];
	  
#ifdef DNPRINT
	  fprintf(out, "Sampling %ld live neutrons (printing from transportcycle.c)\n", nsrc);
#endif

	  /* Parallel loop for sampling live source */
	  /* TODO: set random number indexes in sampledelnu */

#ifdef OPEN_MP
#pragma omp parallel private(id, idx, nn) 
#endif
	  {
	    /* Get Open MP thread id */
		
	    id = OMP_THREAD_NUM;
		
#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif	  
	    /* Loop over source neutrons */
	  
	    for (nn = 0; nn < nsrc; nn++)
	      {
		/* Calculate particle index */

		idx = (long)RDB[DATA_NHIST_TOT];
		idx = idx + (long)(nb*RDB[DATA_SRC_POP]) + nn;

		/* Sample source point */

		SampleSrcPoint(id, nn, idx);

	      }
	  }	  

#ifdef DNPRINT
	  fprintf(out, "Moving to transport neutrons (printing from transportcycle.c)\n");
#endif

	  /* Parallel loop over histories */
	  /* Track neutrons */

	  tosimulate = 1;

	  while (tosimulate > 0)
	    {

	      /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
	      {
		/* Get Open MP thread id */
		    
		id = OMP_THREAD_NUM;
		    
		/* Track neutrons */
		    
		Tracking(id);      
	      }

	      /* Even out ques for next generation */

	      tosimulate = ReDistributeQues();

	    }

	  /* Stop parallel timer */

	  StopTimer(TIMER_OMP_PARA);

	  /* Collect data from interval */

	  CollectDynData();

	  /* Get banked precursors */

	  GetBankedPrecursors();
	  
	  /* Reset normalization coefficients if normalization was not */
	  /* set by dynamic source */
	  
	  if ((long)RDB[DATA_PTR_PREC_DET] < VALID_PTR)
	  {
	    WDB[DATA_NORM_COEF_N] = -1.0;
	    WDB[DATA_NORM_COEF_G] = -1.0;
	  }
	  
	  /* Check for remaining time intervals */

	  if (maxt > 1)
	    {
	      /* Fix normalization to first step */
	      /* With dynamic source, normalization is fixed by source */
	      
	      if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
		NormCoef(PARTICLE_TYPE_NEUTRON);
	      if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
		NormCoef(PARTICLE_TYPE_GAMMA);
	      
	      /* Re-open buffer for writing */

	      WDB[DATA_BUF_REDUCED] = (double)NO;
	      
	      /* Loop over remaining time intervals */
	      
	      for (nt = 1; nt < maxt; nt++)
		{
		  
		  /* Re-open buffer for writing */
		  
		  WDB[DATA_BUF_REDUCED] = (double)NO;

		  /* Set time bin index */

		  WDB[DATA_DYN_TB] = (double)nt;

		  /* Set new time cut-off */

		  WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
		  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];

		  /* Do population control on precursors */

		  PrecursorPopControl();

		  /* Count live neutrons and their weight */

		  CountDynSrc();

		  /* Handle decay of precursor during interval */
		  /* and calculate weight to emit from mesh    */

		  DecayMeshPrecDet();

		  /* Resize the number of live neutrons */
		  /* This is the new population control */

		  ResizeDynSrc();

		  /* Sample delayed neutron source for time interval */
		  /* After this we should have RDB[DATA_SRC_POP]     */
		  /* neutrons */
		  
		  SampleDelnu();

		  /* Re-open buffer for writing (was closed in sampledelnu) */
		  
		  WDB[DATA_BUF_REDUCED] = (double)NO;
		  
		  /* Handle decay of pointwise precursors over interval */

		  DecayPointPrecDet();

		  /* Normalize source */
		  /* Should not do population control with delayed neutrons */
		  /* as the number of neutrons is already correct */
		  
		  NormalizeDynSrc();

		  /* Get number of particles to simulate */

		  tosimulate = ReDistributeQues();
		  
#ifdef DNPRINT
	  fprintf(out, "Moving to transport neutrons (printing from transportcycle.c)\n");
#endif

		  /* Start parallel timer */
		  
		  StartTimer(TIMER_OMP_PARA);
		  
		  /* Loop until source is empty */
		  
#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
		  {
		    /* Get Open MP thread id */
		    
		    id = OMP_THREAD_NUM;
		    
		    /* Loop over source (first generation) */
		    
		    while(FromSrc(id) > VALID_PTR)
		      Tracking(id);      

		  }

		  /* Parallel loop over histories */
		  /* Track neutrons */

		  /* Even out ques for second generation */

		  tosimulate = ReDistributeQues();

		  while (tosimulate > 0)
		    {

		      /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
		      {
			/* Get Open MP thread id */
		    
			id = OMP_THREAD_NUM;
		    
			/* Track neutrons */
		    
			Tracking(id);      
		      }

		      /* Even out ques for next generation */

		      tosimulate = ReDistributeQues();

		    }		

		  /* Collect data from interval */

		  CollectDynData();

		  /* Get banked precursors */
		  
		  GetBankedPrecursors();
		  		  
		  /* Stop parallel timer */
		  
		  StopTimer(TIMER_OMP_PARA);
		}
	    }
	      
	  /* Get end time */

	  t0 = TimerVal(TIMER_TRANSPORT) - t0;
	  c0 = TimerCPUVal(TIMER_TRANSPORT) - c0;

	  /* Score time */

	  ptr = (long)RDB[RES_CYCLE_RUNTIME];
	  AddStat(t0, ptr, 0); 
	  AddStat(c0, ptr, 1); 

	  /* Collect and clear buffered results */

	  CollectResults();
	  CollectDet();
	  CollectPrecDet();
	  PoisonEq();
	  CalcMicroGroupXS();
	  ClearBuf();

	  /* Write savesrc data if requested */

	  WriteDynSrc();

	  /* Flush precursor source */

	  FlushPrecSource();

	  /* Flush bank */

	  FlushBank();

	  /* Stop cycle-wise transport timer */

	  StopTimer(TIMER_TRANSPORT_CYCLE);
	  
	  if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
	    t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
	      TimerVal(TIMER_TRANSPORT_CYCLE);
	  else
	    t0 = 0.0;

	  /* CPU usage */

	  ptr = (long)RDB[RES_CPU_USAGE];
	  AddStat(t0, ptr, 0); 

	  /* Print cycle-wise output */

	  PrintCycleOutput();

	  /* Reset and restart cycle-wise transport timer */

	  ResetTimer(TIMER_TRANSPORT_CYCLE);
	  StartTimer(TIMER_TRANSPORT_CYCLE);

	  /* Sort lists */

	  SortAll();

	  /* Print results */
	  
	  if (!((nb + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
	    {
	      MatlabOutput();
	      DetectorOutput();
	      MeshPlotter();
	      PrintCoreDistr();
	      PrintHistoryOutput();
	      PrintPBData();
	      PrintInterfaceOutput();
	      PrintFinix();
	      PrintPrecDet();
	      RROutput();
	      /*
	      GeometryPlotter(NO);
	      */
	      FissMtxOutput();
	      MORAOutput();
	      WriteICMData();
	    }
	}

      /* Update number of neutron histories run */

      WDB[DATA_NHIST_TOT] = WDB[DATA_NHIST_TOT] +
	RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      /***********************************************************************/

      /***** Dynamic external source simulation ******************************/

      /* Reset skip cycles */

      WDB[DATA_CRIT_SKIP] = 0.0;

      /* Set cycle-wise batch size */

      WDB[DATA_CYCLE_BATCH_SIZE] = RDB[DATA_SRC_POP];
      
      /* Set number of source points in batch */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_SRC_POP];

      /* Get number of bins */

      maxt = (long)RDB[DATA_DYN_NB];
      CheckValue(FUNCTION_NAME, "maxt", "", maxt, 1, 1000000000);

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);
      
      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);    

      /* Set time cut-off for first cycle */
	  
      WDB[DATA_TIME_CUT_TMIN] = RDB[tme];
      WDB[DATA_TIME_CUT_TMAX] = RDB[tme + 1];

      /* Reset batch counters */
	      
      WDB[DATA_BATCH_COUNT] = 0.0;
      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;

      /* Check batching interval */

      if ((long)RDB[DATA_SRC_BATCHES] % (long)RDB[DATA_BATCH_INTERVAL])
	Error(0, 
	    "Total number of batches %ld is not a multiple of interval %ld",
	    (long)RDB[DATA_SRC_BATCHES], (long)RDB[DATA_BATCH_INTERVAL]);	       

      /* Reset time bin index */

      WDB[DATA_DYN_TB] = 0.0;

      /* Stop tracks at outer boundary if time cut-off is used */

      if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
	WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;

      /* Reset cycle k-eff and put starting weight */

      WDB[DATA_CYCLE_KEFF] = 1.0;
      WDB[DATA_DYN_WGT0] = RDB[DATA_SRC_POP];

      /* Reset counter for solution relaxation */

      WDB[DATA_SOL_REL_ITER] = 0.0;

      WDB[DATA_SOL_REL_NTOT] = 0.0;
      WDB[DATA_SOL_REL_NCUR] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];
      WDB[DATA_SOL_REL_N1] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

      /* Set number of batches */

      maxb = (long)(RDB[DATA_SRC_BATCHES]/mpitasks);

      /* Start simulation */

      fprintf(out, "Starting time dependent simulation...\n\n");
      
      /* Start active transport timer */

      StartTimer(TIMER_TRANSPORT_ACTIVE);

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      do
	{
	  /* Prepare coupled calculation iteration */

	  ClearInterfaceStat();
	  
	  /* Re-open buffer for writing */
	  
	  WDB[DATA_BUF_REDUCED] = (double)NO;
	  
	  /* Loop over batches */
	  
	  for (nb = 0; nb < maxb ; nb++)
	    {
	      
	      /* Re-open buffer for writing */
	      
	      WDB[DATA_BUF_REDUCED] = (double)NO;
	      
	      /* Put cycle index */
	      
	      WDB[DATA_CYCLE_IDX] = (double)nb;
	      
	      /* Get beginning time */
	      
	      t0 = TimerVal(TIMER_TRANSPORT);
	      c0 = TimerCPUVal(TIMER_TRANSPORT);
	      
	      /* Start parallel timer */
	      
	      StartTimer(TIMER_OMP_PARA);

	      /* Parallel loop for sampling source */

#ifdef OPEN_MP
#pragma omp parallel private(id, idx, nn) 
#endif
	      {
		/* Get Open MP thread id */
		
		id = OMP_THREAD_NUM;
		
#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif	  
		/* Loop over source neutrons */
	  
		for (nn = 0; nn < (long)RDB[DATA_SRC_POP]; nn++)
		  {
		    /* Calculate particle index */

		    idx = (long)RDB[DATA_NHIST_TOT];
		    idx = idx + (long)(nb*RDB[DATA_SRC_POP]) + nn;

		    /* Sample source point */

		    SampleSrcPoint(id, nn, idx);

		  }
	      }	  

	      /* Parallel loop over histories */
	      /* Track neutrons */

	      tosimulate = (long)RDB[DATA_SRC_POP];

	      while (tosimulate > 0)
		{

		  /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
		  {
		    /* Get Open MP thread id */
		    
		    id = OMP_THREAD_NUM;
		    
		    /* Track neutrons */
		    
		    Tracking(id);      
		  }

		  /* Even out ques for next generation */

		  tosimulate = ReDistributeQues();

		}

	      /* Stop parallel timer */

	      StopTimer(TIMER_OMP_PARA);
	      
	      /* Add to batch counters */

	      WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

	      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 
		RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;
	  
	      /* Check batch interval */

	      if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
		{

		  /* Reset normalization coefficients */
		  /* This way we weill get separate coefficients for       */ 
		  /* separate batches. Normalization coefficients will be  */
		  /* fixed to the batch wise average of the first interval */
	      
		  WDB[DATA_NORM_COEF_N] = -1.0;
		  WDB[DATA_NORM_COEF_G] = -1.0;

		  /* Collect and clear buffered results */
		  
		  CollectDynData();
		  CollectResults();
		  CollectDet();
		  PoisonEq();
		  CalcMicroGroupXS();
		  ClearBuf();

		  /* Reset batch counter */
	      
		  WDB[DATA_BATCH_COUNT] = 0.0;
		}
	      
	      /* Store neutrons at the end of the interval */
	      
	      BanksToStore();
	      
	      /* Flush neutrons from bank before next batch */ 
	      
	      FlushBank();

	      /* Stop cycle-wise transport timer */

	      StopTimer(TIMER_TRANSPORT_CYCLE);
	  
	      if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
		t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
		  TimerVal(TIMER_TRANSPORT_CYCLE);
	      else
		t0 = 0.0;

	      /* CPU usage */

	      ptr = (long)RDB[RES_CPU_USAGE];
	      AddStat(t0, ptr, 0); 

	      /* Print cycle output */
	      
	      PrintCycleOutput();

	      /* Sort lists */

	      SortAll();	    
	      
	      /* Reset and restart cycle-wise transport timer */

	      ResetTimer(TIMER_TRANSPORT_CYCLE);
	      StartTimer(TIMER_TRANSPORT_CYCLE);
	      
	    }
	  
	  /************************************/
	  /* Post-batches processing          */
	  /************************************/

	  /* Collect results from MPI tasks */

	  CollectParallelData();	 

	  /* Iterate coupled calculation */
	  
	  IterateCC();
	  
	  /* Update number of neutron histories run */
	  /* For RNG-seed */

	  WDB[DATA_NHIST_TOT] = WDB[DATA_NHIST_TOT] +
	    RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

	}
      while (RDB[DATA_ITERATE] == (double)YES);

      /* Moving on to remaining timesteps */

      /* Set normalization coefficient as the mean of */
      /* first time interval normalization coefficients */
      
      ptr = (long)RDB[RES_NORM_COEF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      WDB[DATA_NORM_COEF_N] = Mean(ptr,0);
      
      /* Check for remaining time intervals */

      if (maxt > 1)
	{
	  /* Re-open buffer for writing */

	  WDB[DATA_BUF_REDUCED] = (double)NO;
	      
	  /* Loop over remaining time intervals */
	      
	  for (nt = 1; nt < maxt; nt++)
	    {

	      /* Move neutrons from EOI stores to BOI stores */

	      MoveStore();

	      /* Reset counter for solution relaxation */

	      WDB[DATA_SOL_REL_ITER] = 0.0;

	      WDB[DATA_SOL_REL_NTOT] = 0.0;
	      WDB[DATA_SOL_REL_NCUR] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];
	      WDB[DATA_SOL_REL_N1] = RDB[DATA_SRC_BATCHES]*RDB[DATA_SRC_POP];

	      if (RDB[DATA_CC_SIG_MODE] != (double)SIG_MODE_NONE)
		{
		  /* Signal external program about moving to next timestep */

		  SignalExternal(SIGUSR2);
	      
		  /* Read updated interfaces */

		  ptr = (long)RDB[DATA_PTR_IFC0];

		  while (ptr > VALID_PTR)
		    {
		      /* Update interface */

		      ReadInterface(ptr, YES);

		      /* Next interface */

		      ptr = NextItem(ptr);

		    }

		  /* Process updated interfaces */

		  ProcessInterface(YES);
		}

	      /* Iteration loop for next timestep */
	      
	      do
		{
		  /* Prepare coupled calculation iteration */

		  ClearInterfaceStat();
		  
		  /* Re-open buffer for writing */
		  
		  WDB[DATA_BUF_REDUCED] = (double)NO;
		  
		  /* Set time bin index */
		  
		  WDB[DATA_DYN_TB] = (double)nt;
		  
		  /* Set new time cut-off */
		  
		  WDB[DATA_TIME_CUT_TMIN] = RDB[tme + nt];
		  WDB[DATA_TIME_CUT_TMAX] = RDB[tme + nt + 1];

		  /* Reset batch counter */
	      
		  WDB[DATA_BATCH_COUNT] = 0.0;
		  
		  for (nb = 0; nb < maxb ; nb++)
		    {
		      /* Re-open buffer for writing */
		      
		      WDB[DATA_BUF_REDUCED] = (double)NO;
		      
		      /* Put cycle index */
		      
		      WDB[DATA_CYCLE_IDX] = (double)nb;
		      
		      /* Normalize source */
		      
		      if (NormalizeDynSrc() < 0)
			break;
		      
		      /* Start parallel timer */
		      
		      StartTimer(TIMER_OMP_PARA);

		      /* Loop until source is empty */
		  
#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
		      {
			/* Get Open MP thread id */
		    
			id = OMP_THREAD_NUM;
		    
			/* Loop over source (first generation) */
		    
			while(FromSrc(id) > VALID_PTR)
			  Tracking(id);      

		      }

		      /* Stop parallel timer */
		      
		      StopTimer(TIMER_OMP_PARA);

		      /* Parallel loop over histories */
		      /* Track neutrons */

		      /* Even out ques for second generation */

		      tosimulate = ReDistributeQues();

		      while (tosimulate > 0)
			{

			  /* Start parallel timer */
		      
			  StartTimer(TIMER_OMP_PARA);

			  /* Track current generation */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
			  {
			    /* Get Open MP thread id */
		    
			    id = OMP_THREAD_NUM;
		    
			    /* Track neutrons */
		    
			    Tracking(id);      
			  }

			  /* Stop parallel timer */
		      
			  StopTimer(TIMER_OMP_PARA);

			  /* Even out ques for next generation */

			  tosimulate = ReDistributeQues();

			}		
		      
		      /************************************/
		      /* Post-batch processing            */
		      /************************************/		

		      /* Add to batch counters */

		      WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

		      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 
			RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;

		      /* Check batch interval */

		      if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
			{

			  /* Collect and clear buffered results */

			  CollectDynData();
			  CollectResults();
			  CollectDet();
			  PoisonEq();
			  CalcMicroGroupXS();
	      
			  /* Clear buffers after each batch */
	      
			  ClearBuf();

			  /* Reset batch counter */
	      
			  WDB[DATA_BATCH_COUNT] = 0.0;
			}

		      /* Store neutrons at the end of the interval */
	      
		      BanksToStore();
	      
		      /* Flush neutrons from bank before next batch */ 
	      
		      FlushBank();

		      /* Stop cycle-wise transport timer */

		      StopTimer(TIMER_TRANSPORT_CYCLE);

		      if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
			t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
			  TimerVal(TIMER_TRANSPORT_CYCLE);
		      else
			t0 = 0.0;

		      /* CPU usage */

		      ptr = (long)RDB[RES_CPU_USAGE];
		      AddStat(t0, ptr, 0); 

		      /* Print cycle output */
	      
		      PrintCycleOutput();

		      /* Sort lists */

		      SortAll();

		      /* Reset and restart cycle-wise transport timer */

		      ResetTimer(TIMER_TRANSPORT_CYCLE);
		      StartTimer(TIMER_TRANSPORT_CYCLE);
		    }

		  /* After all batches have been calculated */

		  /* Collect results from MPI tasks */

		  CollectParallelData();	 

		  /* Iterate coupled codes */

		  IterateCC();
		  
		  /* Check iteration flag */
		  
		}
	      while(RDB[DATA_ITERATE] == (double)YES);

	      /* Print results */
	  
	      if (!((nt + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
		{
		  MatlabOutput();
		  DetectorOutput();
		  MeshPlotter();
		  PrintCoreDistr();
		  PrintHistoryOutput();
		  PrintPBData();
		  PrintInterfaceOutput();
		  RROutput();
		  PrintFinix();
		  /*
		    GeometryPlotter(NO);
		  */
		  FissMtxOutput();
		  MORAOutput();
		  WriteICMData();
		}
	    }
	}

      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /***********************************************************************/
      
      /***** Neutron criticality source simulation ***************************/

      /* Reset simulated batch size */

      WDB[DATA_SIMUL_BATCH_SIZE] = 0.0;

      /* Reset time bin index */

      WDB[DATA_DYN_TB] = 0.0;

      /* Generate initial source */

      if((RDB[DATA_USE_FSP] == (double)NO) || 
	 ((RDB[DATA_BURN_STEP] + RDB[DATA_SOL_REL_ITER] == 0.0)
         && (RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)))
	{

	  fprintf(out, "Sampling initial source...\n");

#ifdef OPEN_MP
#pragma omp parallel private(id, idx) 
#endif
	  {
	    /* Get Open MP thread id */
	
	    id = OMP_THREAD_NUM;
	
#ifdef OPEN_MP
#pragma omp for schedule(dynamic)
#endif	  
	    /* Loop over source neutrons */
	  
	    for (nn = 0; nn < (long)RDB[DATA_CRIT_POP]; nn++)
	      {
		/* Calculate particle index */

		idx = (long)RDB[DATA_NHIST_TOT] + nn;

		/* Sample source point */
		
		SampleSrcPoint(id, nn, idx);
	      }
	  }

	  fprintf(out, "OK.\n\n");

          /* Set number of inactive batches */

	  skip = (long)RDB[DATA_CRIT_SKIP];
	  nb0 = 0;
	}
      else
	{
	  fprintf(out, "Continuing from previous fission source\n\n");	

	  ResizeFissionSrc();

          /* Set number of inactive batches                          */
	  /* by setting the initial inactive batch number            */
	  /* This should make sure that the output and processing    */
	  /* is done after the same number of live batches as before */

	  nb0 = (long)RDB[DATA_CRIT_SKIP] - (long)RDB[DATA_FSP_CRIT_SKIP];

	  skip = (long)RDB[DATA_CRIT_SKIP];
	}

      /* Start cycle-wise transport timer */

      ResetTimer(TIMER_TRANSPORT_CYCLE);
      StartTimer(TIMER_TRANSPORT_CYCLE);

      /* Set number of batches */

#ifdef MPI_MODE1
      
      maxb = (long)RDB[DATA_CRIT_CYCLES];

#else

      maxb = (long)(RDB[DATA_CRIT_CYCLES]/mpitasks);

#endif

      /* Check batching interval */

      if (maxb % (long)RDB[DATA_BATCH_INTERVAL])
	Error(0,
	    "Total number of batches %ld is not a multiple of interval %ld",
	    maxb, (long)RDB[DATA_BATCH_INTERVAL]);

      /* Loop over batches */
      
      for (nb = nb0; nb < maxb + skip; nb++)
	{
	  /* Check number of skip cycles */

	  if (nb == skip)
	    {
	      /* Clear statistics */
	      
	      ClearStat(-1);

	      /* Start active transport timer */

	      StartTimer(TIMER_TRANSPORT_ACTIVE);

	      /* Reset number of active neutron histories */

	      WDB[DATA_NHIST_CYCLE] = 0.0;

	      /* Reset batch counters */
	      
	      WDB[DATA_BATCH_COUNT] = 0.0;
	      WDB[DATA_MICRO_CALC_BATCH_COUNT] = 0.0;
	    }

	  /* Put cycle index */

	  WDB[DATA_CYCLE_IDX] = (double)nb;

	  /* Normalize source */

	  NormalizeCritSrc();

	  /* Clear events */

	  ProcessEvents();

	  /* Get beginning time */

	  t0 = TimerVal(TIMER_TRANSPORT);
	  c0 = TimerCPUVal(TIMER_TRANSPORT);

	  /* Start parallel timer */

	  StartTimer(TIMER_OMP_PARA);

	  /* Loop until source is empty */

#ifdef OPEN_MP
#pragma omp parallel private(id)
#endif
	  {
	    /* Get Open MP thread id */

	    id = OMP_THREAD_NUM;

	    /* Loop over source */

	    while(FromSrc(id) > VALID_PTR)
	      Tracking(id);      
	  }

	  /* Stop parallel timer */

	  StopTimer(TIMER_OMP_PARA);

	  /* Get end time */

	  t0 = TimerVal(TIMER_TRANSPORT) - t0;
	  c0 = TimerCPUVal(TIMER_TRANSPORT) - c0;
	  
	  /* Score time */

	  ptr = (long)RDB[RES_CYCLE_RUNTIME];
	  AddStat(t0, ptr, 0); 
	  AddStat(c0, ptr, 1); 

	  /* Reset normalization coefficients */

	  WDB[DATA_NORM_COEF_N] = -1.0;
	  WDB[DATA_NORM_COEF_G] = -1.0;

	  /* K-eff iteration */

	  IterateKeff();

	  /* Add to batch counters */

	  WDB[DATA_BATCH_COUNT] = RDB[DATA_BATCH_COUNT] + 1.0;

	  WDB[DATA_MICRO_CALC_BATCH_COUNT] = 
	    RDB[DATA_MICRO_CALC_BATCH_COUNT] + 1.0;
	  
	  /* Check batch interval */

	  if ((long)RDB[DATA_BATCH_COUNT] == (long)RDB[DATA_BATCH_INTERVAL])
	    {
	      /* Collect and clear buffered results */

	      CollectResults();
	      CalcMicroGroupXS();
	      CollectPrecDet();
	      CollectDet();
	      PoisonEq();
	      ClearBuf();

	      /* Reset batch counter */
	      
	      WDB[DATA_BATCH_COUNT] = 0.0;
	    }
	  
	  /* Stop cycle-wise transport timer */

	  StopTimer(TIMER_TRANSPORT_CYCLE);
	  
	  if (TimerVal(TIMER_TRANSPORT_CYCLE) > 0.0)
	    t0 = TimerCPUVal(TIMER_TRANSPORT_CYCLE)/
	      TimerVal(TIMER_TRANSPORT_CYCLE);
	  else
	    t0 = 0.0;

	  /* CPU usage */

	  ptr = (long)RDB[RES_CPU_USAGE];
	  AddStat(t0, ptr, 0); 

	  /* Print cycle-wise output */

	  PrintCycleOutput();

	  /* Reset and restart cycle-wise transport timer */

	  ResetTimer(TIMER_TRANSPORT_CYCLE);
	  StartTimer(TIMER_TRANSPORT_CYCLE);

	  /* Sort lists */

	  SortAll();

	  /* Print results */
	  
	  if (!((nb + 1) % (long)RDB[DATA_PRINT_INTERVAL]))
	    {
	      MatlabOutput();
	      DetectorOutput();
	      MeshPlotter();
	      PrintCoreDistr();
	      PrintHistoryOutput();
	      PrintPBData();
	      PrintFinix();
	      PrintPrecDet();
	      PrintInterfaceOutput();
	      RROutput();
	      /*
	      GeometryPlotter(NO);
	      */
	      FissMtxOutput();
 	      MORAOutput();
	      WriteICMData();
	    }
	}

      /* Flush bank if not passing it to next step */
      
      if(RDB[DATA_USE_FSP] == (double)NO)     
	FlushBank();
    
      /* Stop cycle-wise transport timer */

      StopTimer(TIMER_TRANSPORT_CYCLE);

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid mode");
       
  /***************************************************************************/

  /***** Transport cycle completed *******************************************/

  /* Collect results from MPI tasks */

  CollectParallelData();

  /* Put completed flag */

  WDB[DATA_SIMULATION_COMPLETED] = (double)YES;

  /* Dump buffered source points into file */

  WriteSourceFile(-1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -1);

  /* Clear scoring buffer */

  ClearBuf();

  /* Stop transport timers */

  StopTimer(TIMER_TRANSPORT);
  StopTimer(TIMER_TRANSPORT_TOTAL);

  if (((long)RDB[DATA_CRIT_CYCLES] > 0) || ((long)RDB[DATA_SRC_BATCHES] > 0))
    StopTimer(TIMER_TRANSPORT_ACTIVE);

  /* Remember previous value */

  if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
    WDB[DATA_PRED_TRANSPORT_TIME] = TimerVal(TIMER_TRANSPORT);
  else
    WDB[DATA_CORR_TRANSPORT_TIME] = TimerVal(TIMER_TRANSPORT);

  /* Set estimate for coefficient calculation */

  WDB[DATA_COEF_TRANSPORT_TIME] = TimerVal(TIMER_TRANSPORT);

  /* Put keff for burnup iteration */

  if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == YES)
    ptr = (long)RDB[RES_IMP_KEFF];
  else
    ptr = (long)RDB[RES_COL_KEFF];  

  WDB[DATA_BURN_PREV_KEFF] = Mean(ptr, 0);
  WDB[DATA_BURN_PREV_DKEFF] = StdDev(ptr, 0);

  /* Put poison concentrations */

  PutPoisonConc();

  /* Calculate activities (poison concentrations may be updated) */

  CalculateActivities();

  /* Print output */

  MatlabOutput();
  DetectorOutput();
  MeshPlotter();
  PrintCoreDistr();
  PrintHistoryOutput();
  PrintPBData();
  PrintFinix();
  PrintPrecDet();
  PrintInterfaceOutput();
  RROutput();

  GeometryPlotter(NO);

  FissMtxOutput();
  MORAOutput();
  WriteICMData();
  /*
  ARESOutput();
  */

  /* Statistical tests */

  StatTests();

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
