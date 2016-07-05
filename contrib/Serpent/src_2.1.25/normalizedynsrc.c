/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : normalizedynsrc.c                              */
/*                                                                           */
/* Created:       2012/09/23 (JLe)                                           */
/* Last modified: 2016/02/01 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Normalizes source distribution for dynamic criticality       */
/*              source simulation                                            */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NormalizeDynSrc:"

/*****************************************************************************/

long NormalizeDynSrc()
{
  long ptr, pos, part, n, m, nsrc, nbatch, mul, N, np, id, idx;
  long nb, ptr2, i, new, stp, loc0;
  double wgt0, wgt, P;

  /***************************************************************************/

  /***** Move banked neutrons to source **************************************/

  /* Get pointer to source */

  pos = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(pos)", DATA_ARRAY, pos);
  
  /* Check that source is empty */

  if (ListSize(pos) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Reset total weight and source size */

  wgt0 = 0.0;
  nsrc = 0;

  /* In DYN mode get particles from store instead of bank */

  if (RDB[DATA_SIMULATION_MODE] == (double)SIMULATION_MODE_DYN)
    {

      /* Get current batch number */

      nb = (long)RDB[DATA_CYCLE_IDX];

      ptr2 = (long)RDB[DATA_PTR_DYN_PARTCOUNT];

      /* Loop over threads */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	{

	  /* Get particles from BOI store */

	  idx = nb + id*(long)RDB[DATA_SRC_BATCHES];

	  np = (long)RDB[ptr2 + idx];

	  for(i = 0 ; i < np ; i++)
	    {
	      ptr = FromStore(id,0);

	      /* Check type */
          
	      if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
		Die(FUNCTION_NAME, "Invalid particle type");

	      /* Copy particle */
          
	      new = DuplicateParticle(ptr,id);

	      /* Put the original back to BOI store */

	      ToStore(ptr, id, 0);

	      /* Add to total weight and source size */

	      wgt0 = wgt0 + RDB[new + PARTICLE_WGT];
	      nsrc = nsrc + 1;

	      /* Put neutron in source */
          
	      if ((long)RDB[DATA_OPTI_OMP_REPRODUCIBILITY] == YES)
		AddSortItem(DATA_PART_PTR_SOURCE, pos, new, PARTICLE_RNG_IDX, 
			    SORT_MODE_ASCEND);
	      else
		AddItem(DATA_PART_PTR_SOURCE, new);
          
	      /* Update position */
          
	      pos = new;

	    }

	}

    }
  else
    {
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
	  
	      wgt0 = wgt0 + RDB[ptr + PARTICLE_WGT];
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
    }

  /* Check if source is empty */

  if(nsrc == 0)
    {
      Warn(FUNCTION_NAME, "Source is empty");
      return -1;
    }

  /* Check extinction and divergence */

  if (wgt0/RDB[DATA_SRC_POP] < 1E-16)
    {
      /* Reset id */

      id = 0;

      /* Pointer to first item after dummy */
      
      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);
      
      /* Loop over source */
      
      while (ptr > VALID_PTR)
	{
	  /* Get poiner */

	  part = ptr;

	  /* Pointer to Next */
      
	  ptr = NextItem(ptr);

	  /* Remove particle */

	  RemoveItem(part);

	  /* Put particle back to stack */

	  ToStack(part, id++);

	  /* Check id */

	  if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
	    id = 0;
	}

      /* Exit subroutine */
      
      return -1;
    }
  else if (wgt0/RDB[DATA_SRC_POP] > 1E+16)
    Error(0,"Population exceeds maximum at %1.2E seconds, adjust time cut-off",
	  WDB[DATA_TIME_CUT_TMAX]);

  /* Set batch size */

  nbatch = (long)RDB[DATA_SRC_POP];

  /* Try to get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
    {
      /* If precursor detector is in use, population control is done in */
      /* ResizeDynSrc and SampleDelnu */

      /* Set nbatch to equal nsrc to not do population control again */
      /* Due to the stochastic nature of the population control, nsrc might */
      /* not be exactly DATA_SRC_POP */

      nbatch = nsrc;
    }
    
  /***************************************************************************/

  /***** Population control **************************************************/

  /* Reset weight and number of particles */

  wgt = wgt0;
  np = nsrc;

  /* Compare population size to batch size given in input */

  if (nsrc < nbatch)
    {
      /* Calculate multiplication */

      P = ((double)nbatch)/((double)nsrc);
      mul = (long)P;
      P = P - (double)mul;

      /* Pointer to first neutron after dummy */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);

      /* Loop over population and sample particles for duplication */

      for (n = 0; n < nsrc; n++)
	{
	  /* Check pointer */

	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Sample multiplication */

	  if (RandF(0) < P)
	    N = mul;
	  else
	    N = mul - 1;

	  /* Loop over multiplication */

	  for (m = 0; m < N; m++)
	    {
	      /* Duplicate neutron */
	  
	      part = DuplicateParticle(ptr, 0);

	      /* Add to weight and population size */

	      wgt = wgt + RDB[ptr + PARTICLE_WGT];
	      np = np + 1;

	      /* Put particle in source*/
	  
	      AddItem(DATA_PART_PTR_SOURCE, part);
	    }

	  /* Next particle */

	  ptr = NextItem(ptr);
	}
    }
  else if (nsrc > nbatch)
    {
      /* Calculate probability */

      P = 1.0 - ((double)nbatch)/((double)nsrc);

      /* Pointer to first neutron after dummy */

      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
      ptr = NextItem(ptr);

      /* Reset id */

      id = 0;

      /* Loop over population and sample particles for duplication */

      for (n = 0; n < nsrc; n++)
	{
	  /* Check pointer */

	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Get pointer and next */
	  
	  part = ptr;
	  ptr = NextItem(ptr);

	  /* Sample removal */

	  if (RandF(0) < P)
	    {
	      /* Remove particle */

	      RemoveItem(part);

	      /* subtract from weight and population size */
	      
	      wgt = wgt - RDB[part + PARTICLE_WGT];
	      np = np - 1;

	      /* Put particle back to stack */

	      ToStack(part, id++);

	      /* Check id */

	      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
		id = 0;
	    }
	}
    }

  /* Compare population size to minimum and maximum */

  if (nsrc < (long)RDB[DATA_DYN_POP_MIN])
    WDB[DATA_DYN_POP_MIN] = nsrc;
  if (nsrc > (long)RDB[DATA_DYN_POP_MAX])
    WDB[DATA_DYN_POP_MAX] = nsrc;

  /* Check */

  if (np < 1)
    return -1;

  /***************************************************************************/

  /***** Normalize source ****************************************************/

  /* Pointer to first item after dummy */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  ptr = NextItem(ptr);

  /* Normalize weights */

  while(ptr > VALID_PTR)
    {
      /* Normalize */

      WDB[ptr + PARTICLE_WGT] = RDB[ptr + PARTICLE_WGT]*wgt0/wgt;
      
      /* Next */
      
      ptr = NextItem(ptr);
    }

  /***************************************************************************/

  /***** Set particle indexes ************************************************/

  /* Pointer to source buffer */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get next rng index */

  n = (long)RDB[DATA_NHIST_TOT];

  /* Reset total weight and source size */

  wgt = 0.0;
  nsrc = 0;

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

      /* Add to total weight and source size */
      
      wgt = wgt + RDB[ptr + PARTICLE_WGT];
      nsrc = nsrc + 1;

      /* Score initial source rate and source weight for */
      /* interval in case of dynamic mode                */
      /* For MODE_SRC these are scored in samplesrcpoint */      

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
	{

	  /* Check particle type */

	  if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_GAMMA)
	    {
	      /* Score source rate */

	      stp = (long)RDB[RES_TOT_PHOTON_SRCRATE];  
	      CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	      AddBuf1D(1.0, wgt, stp, id, 0);
	    }
	  else
	    {
	      /* Score source rate */

	      stp = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
	      CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	      AddBuf1D(1.0, RDB[ptr + PARTICLE_WGT], stp, 0, 0);
	  
	      /* Score initial source weight */
	  
	      stp = (long)RDB[RES_INI_SRC_WGT];
	      CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	      AddBuf1D(1.0, RDB[ptr + PARTICLE_WGT], stp, 0, 0);
	  
	    }
	}
      /* Next particle */

      ptr = PrevItem(ptr);
    }

  /* Check weight */

  if (fabs(wgt/wgt0 - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Mismatch in weight");

  /* Put number of histories */

  WDB[DATA_NHIST_TOT] = (double)n;

  /* Score mean population size */

  ptr = (long)RDB[RES_MEAN_POP_SIZE];
  AddStat((double)nsrc, ptr, 0); 

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

  /* Exit */

  return 0;

  /***************************************************************************/
}

/*****************************************************************************/
