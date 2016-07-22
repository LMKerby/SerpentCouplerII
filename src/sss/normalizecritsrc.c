/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normalizecritsrc.c                             */
/*                                                                           */
/* Created:       2011/03/10 (JLe)                                           */
/* Last modified: 2016/02/17 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Sets up normalized fission source for criticality source     */
/*              simulation                                                   */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormalizeCritSrc:"

/*****************************************************************************/

void NormalizeCritSrc()
{
  long ptr, pos, n, mat, stp, nsrc, nbatch, id, fmx, idx;
  double wgt, w0, keff, kw, P, kp;
  
  /***************************************************************************/

  /***** Move banked neutrons to source **************************************/

  /* Get pointer to source */

  pos = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(pos)", DATA_ARRAY, pos);
  
  /* Check that source is empty */

  if (ListSize(pos) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Reset total weight and source size */
  
  wgt = 0.0;
  nsrc = 0;

  /* Loop over threads */
  
  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      while ((ptr = FromBank(id)) > VALID_PTR)
	{
	  /* Check type */
	  
	  if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
	    Die(FUNCTION_NAME, "Invalid particle type");
	  
	  /* Add to total weight and source size */
	  
	  wgt = wgt + RDB[ptr + PARTICLE_WGT];
	  nsrc = nsrc + 1;

	  /* Put neutron in source */
	  
	  if ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == YES)
	    AddSortItem(DATA_PART_PTR_SOURCE, pos, ptr, PARTICLE_RNG_IDX, 
			SORT_MODE_ASCEND);
	  else
	    AddItem(DATA_PART_PTR_SOURCE, ptr);
	  
	  /* Update position */
	  
	  pos = ptr;
	}
    }

  /* Check weight */

  if (wgt == 0.0)
    Error(0, "Unable to initiate self-sustaining chain reaction");

  /* Get data from MPI parallel tasks */

  Rendezvous(&nsrc, &wgt);

  /* Set batch size */

  if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
    nbatch = (long)RDB[DATA_SIMUL_BATCH_SIZE];
  else
    nbatch = (long)RDB[DATA_CRIT_POP];

  /* Set cycle batch size and weight */

  WDB[DATA_CYCLE_BATCH_SIZE] = (double)nsrc;

  /* Score mean population size */

  ptr = (long)RDB[RES_MEAN_POP_SIZE];
  AddStat(nsrc, ptr, 0); 

  /* Score mean population weight */

  ptr = (long)RDB[RES_MEAN_POP_WGT];
  AddStat(wgt, ptr, 0); 

  /* Get previous cycle-wise k-eff */

  keff = RDB[DATA_CYCLE_KEFF];

  /* Calculate new cycle-wise k-eff */

  keff = keff*wgt/((double)(nbatch));

  /* Put cycle-wise k-eff */
      
  WDB[DATA_CYCLE_KEFF] = keff;

  /* Check if Wieland shift is used */

  if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
    {
      /* Avoid compiler error */

      kw = -1.0;
      kp = -1.0;
      P = -1.0;

      /* Check mode */

      if ((long)RDB[DATA_WIELANDT_MODE] == WIELANDT_MODE_FIX_K)
	{
	  /* Get user-specified k-eff */

	  kw = RDB[DATA_WIELANDT_KEFF];
	  CheckValue(FUNCTION_NAME, "kw", "", kw, ZERO, INFTY);
	  
	  /* Calculate backtransformed k-eff */
	  
	  kp = kw*keff/(kw + keff);
	  CheckValue(FUNCTION_NAME, "kp", "", kp, ZERO, INFTY);

	  /* Calculate probability of banking the neutron */
      
	  P = 1.0 - kp/kw;
	  CheckValue(FUNCTION_NAME, "P", "", P, ZERO, 1.0);
	}
      else if ((long)RDB[DATA_WIELANDT_MODE] == WIELANDT_MODE_FIX_P)
	{
	  /* Get user-specified probability */

	  P = RDB[DATA_WIELANDT_P];
	  CheckValue(FUNCTION_NAME, "P", "", P, ZERO, 1.0);

	  /* Calculate backtransformed k-eff */

	  kp = keff*P;
	  CheckValue(FUNCTION_NAME, "kp", "", kp, ZERO, INFTY);

	  /* Calculate k-eff */

	  kw = kp/(1.0 - P);
	  CheckValue(FUNCTION_NAME, "kw", "", kw, ZERO, INFTY);
	}
      else
	Die(FUNCTION_NAME, "Invalid Wielandt mode");

      /* Store values */

      WDB[DATA_WIELANDT_KEFF] = kw;
      WDB[DATA_WIELANDT_KP] = kp;
      WDB[DATA_WIELANDT_P] = P;

      /* Store kw and P */

      ptr = (long)RDB[RES_WIELANDT_K];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(kw, 1.0, ptr, 0, 0);

      ptr = (long)RDB[RES_WIELANDT_P];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(P, 1.0, ptr, 0, 0);

      /* Store backtransformed keff */

      ptr = (long)RDB[RES_ANA_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(kp, 1.0, ptr, 0, 0);
    }
  else
    {
      /* Store keff */

      ptr = (long)RDB[RES_ANA_KEFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(keff, 1.0, ptr, 0, 0);
    }

  /***************************************************************************/

  /***** Re-normalize source *************************************************/

  /* Reset prompt and delayed weights */

  WDB[DATA_CYCLE_PROMPT_WGT] = 0.0;
  WDB[DATA_CYCLE_DELAYED_WGT] = 0.0;

  /* Pointer to first item after dummy */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  ptr = NextItem(ptr);

  /* Loop over source */

  while (ptr > VALID_PTR)
    {
      /* Normalize weight */

      w0 = RDB[ptr + PARTICLE_WGT]*((double)nbatch)/wgt;

      /* Set value */

      WDB[ptr + PARTICLE_WGT] = w0;

      /* Check delayed neutron group */

      if ((long)RDB[ptr + PARTICLE_DN_GROUP] == 0)
	WDB[DATA_CYCLE_PROMPT_WGT] = RDB[DATA_CYCLE_PROMPT_WGT] + w0;
      else
	WDB[DATA_CYCLE_DELAYED_WGT] = RDB[DATA_CYCLE_DELAYED_WGT] + w0;

      /* Get material pointer (may be null for initial source) */

      mat = (long)RDB[ptr + PARTICLE_PTR_MAT];

      /* Score source rate */
      
      if ((mpiid == 0) || ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO))
	{
	  stp = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
	  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	  AddBuf1D(1.0, w0, stp, 0, 0);

	  /* Score source rate in fissile and non-fissile materials */
	  
	  if (mat > VALID_PTR)
	    {
	      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
		AddBuf1D(1.0, w0, stp, 0, 1);
	      else
		AddBuf1D(1.0, w0, stp, 0, 2);
	    }
	}

      /* Score fission matrix source term */

      if ((idx = (long)RDB[ptr + PARTICLE_FMTX_IDX]) > -1)
	{
	  /* Get pointer */
	  
	  fmx = (long)RDB[DATA_PTR_FMTX];
	  CheckPointer(FUNCTION_NAME, "(fmx)", DATA_ARRAY, fmx);
	  
	  stp = (long)RDB[fmx + FMTX_PTR_SRC];
	  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	  
	  /* Score total */
 
	  AddBuf(1.0, w0, stp, 0, -1, 0, idx);
	  
	  /* Score prompt or delayed */

	  if ((long)RDB[ptr + PARTICLE_DN_GROUP] == 0)
	    AddBuf(1.0, w0, stp, 0, -1, 1, idx);
	  else
	    AddBuf(1.0, w0, stp, 0, -1, 2, idx);
	}
      
      /* Next */
      
      ptr = NextItem(ptr);
    }

  /* Calculate entropies */

  CalculateEntropies();
 
  /***************************************************************************/

  /***** Set particle indexes ************************************************/

  /* Get next rng index */

  n = (long)RDB[DATA_NHIST_TOT];

  /* Reset MPI index counter */

  id = 0;

  /* Reset weight */
  
  wgt = 0.0;

  /* Pointer to source buffer */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get pointer to last particle */

  ptr = LastItem(ptr);

  /* Loop over list and set indexes */

  while(1 != 2)
    {
      /* Break if dummy */

      if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
	break;

      /* Put index */

      WDB[ptr + PARTICLE_RNG_IDX] = (double)(++n);
      WDB[ptr + PARTICLE_HISTORY_IDX] = RDB[ptr + PARTICLE_RNG_IDX];

      /* Add to weight */

      wgt = wgt + RDB[ptr + PARTICLE_WGT];

      /* Put MPI id */

      if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
	WDB[ptr + PARTICLE_MPI_ID] = (double)mpiid;
      else
	WDB[ptr + PARTICLE_MPI_ID] = (double)(id++);

      /* Check id */

      if (id == mpitasks)
	id = 0;

      /* Update number of histories */

      WDB[DATA_NHIST_CYCLE] = RDB[DATA_NHIST_CYCLE] + 1.0;

      /* Next particle */

      ptr = PrevItem(ptr);
    }

  /* Check total weight */

  if(fabs(wgt/((double)nbatch) - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Total weight not preserved");

  /* Put number of histories */

  WDB[DATA_NHIST_TOT] = (double)n;

  /* Check minimum stack size */

  if (RDB[DATA_OMP_MAX_THREADS]*RDB[DATA_PART_MIN_NSTACK]/
      RDB[DATA_PART_ALLOC_N] < 0.2)
    {
      /* Allow memory allocation */
      
      Mem(MEM_ALLOW);

      /* Allocate 20% more particles */

      AllocParticleStack(PARTICLE_TYPE_NEUTRON, 
			 (long)(0.2*RDB[DATA_PART_ALLOC_N]));

      /* Disallow memory allocation */

      Mem(MEM_DENY);
    }

  /* Re-distribute stacks */

  ReDistributeStacks();

  /* Check that source is sorted */

#ifdef DEBUG
  
  /* Check reproducibility */

  if ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == YES)
    {
      /* Pointer to first after dummy */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);

      /* Loop over source */

      while (ptr > VALID_PTR)
	{
	  /* Compare */
	  
	  if ((pos = NextItem(ptr)) > VALID_PTR)
	    if (RDB[pos + PARTICLE_RNG_IDX] >= RDB[ptr + PARTICLE_RNG_IDX])
	      Die(FUNCTION_NAME, "Sorting failed");

	  /* Next */

	  ptr = NextItem(ptr);
	}
    }
 
#endif

  /* Plot source point distribution */

  GeometryPlotter(NO);

  /***************************************************************************/
}

/*****************************************************************************/
