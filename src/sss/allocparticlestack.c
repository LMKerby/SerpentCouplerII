/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : allocparticlestacs.c                           */
/*                                                                           */
/* Created:       2012/10/17 (JLe)                                           */
/* Last modified: 2016/01/31 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Allocates memory for particle histories (stacks)             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AllocParticleStack:"

/*****************************************************************************/

void AllocParticleStack(long type, long np)
{
  long ptr, loc0, id, n, m;

#ifdef OLD_HIST

  long hst;

#endif

  /* Avoid compiler warning */

  loc0 = -1;

  /* Check type */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      /* Get pointer */

      loc0 = DATA_PART_PTR_NSTACK;

      /* Update stack size */

      WDB[DATA_PART_ALLOC_N] = RDB[DATA_PART_ALLOC_N] + (double)np;

      /* Reset minimum size */

      WDB[DATA_PART_MIN_NSTACK] = RDB[DATA_PART_ALLOC_N];
    }
  else if (type == PARTICLE_TYPE_GAMMA)
    {
      /* Get pointer */

      loc0 = DATA_PART_PTR_GSTACK;

      /* Update stack size */

      WDB[DATA_PART_ALLOC_G] = RDB[DATA_PART_ALLOC_G] + (double)np;

      /* Reset minimum size */

      WDB[DATA_PART_MIN_GSTACK] = RDB[DATA_PART_ALLOC_G];
    }
  else if (type == PARTICLE_TYPE_PRECURSOR)
    {
      /* Get pointer */

      loc0 = DATA_PART_PTR_PSTACK;

      /* Update stack size */

      WDB[DATA_PART_ALLOC_P] = RDB[DATA_PART_ALLOC_P] + (double)np;

      /* Reset minimum size */

      WDB[DATA_PART_MIN_PSTACK] = RDB[DATA_PART_ALLOC_P];
    }
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Reset OpenMP index */
      
  id = 0;
  
  /* Loop over particles */
      
  for (n = 0; n < np; n++)
    {
      /* Allocate memory */
      
      ptr = NewItem(OMPPtr(loc0, id), PARTICLE_BLOCK_SIZE);

      /* Put type */
	  
      WDB[ptr + PARTICLE_TYPE] = (double)type;

      /* Allocate memory for fission progenies */

#ifdef OLD_IFP

      if (type == PARTICLE_TYPE_NEUTRON)
	for (m = 0; m < (long)RDB[DATA_IFP_CHAIN_LENGTH]; m++)
	  NewItem(ptr + PARTICLE_PTR_FISS_PROG, FISS_PROG_BLOCK_SIZE);
#endif	  

      /* Allocate memory for history data (tohon NHIST) */

#ifdef OLD_HIST
	  
      if ((long)RDB[DATA_HIST_LIST_SIZE] > 0)
	{
	  /* Loop over events */

	  for (m = 0; m < (long)RDB[DATA_HIST_LIST_SIZE]; m++)
	    {
	      /* Allocate memory for data */
	      
	      hst = NewItem(ptr + PARTICLE_PTR_HIST, HIST_BLOCK_SIZE);
	      
	      /* Reset weight to indicate unused value */
	      
	      WDB[hst + HIST_WGT] = -1.0;
	    }
	  
	  /* Make list into ring */
	  
	  hst = (long)RDB[ptr + PARTICLE_PTR_HIST];
	  MakeRing(hst);
	}

#endif
      
      /* Update OpenMP id */
      
      if (++id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
	id = 0;
    }

  /* Update memory size */
  
  WDB[DATA_TOT_MISC_BYTES] = RDB[DATA_TOT_MISC_BYTES] + MemCount();
}

/*****************************************************************************/
