/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : redistributestacs.c                            */
/*                                                                           */
/* Created:       2012/10/13 (JLe)                                           */
/* Last modified: 2016/02/16 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Redistributes stacked particles between OpenMp threads       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReDistributeStacks:"

/*****************************************************************************/

void ReDistributeStacks()
{
  long m, n, type, loc0, loc1, ptr, id, id0, id1, min, max, sz, tot, nt;

  /* Get number of OpenMp threads */

  nt = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Reset total counts */

  tot = 0;

  /* Loop over particle types */

  for (m = 0; m < 3; m++)
    {
      /* Get pointer */

      if (m == 0)
	{
	  type = PARTICLE_TYPE_NEUTRON;
	  loc0 = DATA_PART_PTR_NSTACK;
	}
      else if (m == 1)
	{
	  type = PARTICLE_TYPE_GAMMA;
	  loc0 = DATA_PART_PTR_GSTACK;
	}
      else 
	{
	  type = PARTICLE_TYPE_PRECURSOR;
	  loc0 = DATA_PART_PTR_PSTACK;
	}

      /* Check pointer */

      if ((long)RDB[loc0] < VALID_PTR)
	Die(FUNCTION_NAME, "Pointer error");

      /* Loop */

      while (1 != 2)
	{
	  /* Find stack with lowest and highest number of particles */
	  
	  id0 = -1;
	  id1 = -1;
	  
	  min = 100000000000;
	  max = -1;
	  
	  for (id = 0; id < nt; id++)
	    {
	      /* Get list size */
	      
	      loc1 = (long)RDB[OMPPtr(loc0, id)];
	      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	      sz = ListSize(loc1) - 1;
	      
	      /* Compare to minimum */
	      
	      if (sz < min)
		{
		  min = sz;
		  id0 = id;
		}
	      
	      /* Compare to maximum */
	      
	      if (sz > max)
		{
		  max = sz;
		  id1 = id;
		}
	    }
	  
	  /* Check condition */

	  if ((sz = max - min) < 5*nt)
	    break;
	  else
	    {
	      /* Move half of particles to other stack */
	      
	      for (n = 0; n < (long)(0.5*sz); n++)
		{
		  ptr = FromStack(type, id1);
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  ToStack(ptr, id0);
		}
	    }
	}

      /* Add stacks to total count */

      for (id = 0; id < nt; id++)
	{
	  /* Add stack size */
	  
	  ptr = (long)RDB[OMPPtr(loc0, id)];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Do not count precursor stacks */

	  if (m < 2)
	    tot = tot + ListSize(ptr) - 1;
	}
    }

  /* Add ques to total count */

  for (id = 0; id < nt; id++)
    {
      /* Add que size */
      
      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_QUE, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      tot = tot + ListSize(ptr) - 1;
    }

  /* Add stores to total count */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {

      for (id = 0; id < nt; id++)
	{
	  /* Add BOI store size */
      
	  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BOI_STORE, id)];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  tot = tot + ListSize(ptr) - 1;

	  /* Add EOI store size */
      
	  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_EOI_STORE, id)];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  tot = tot + ListSize(ptr) - 1;
	}

    }
  /* Add source to total count */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  tot = tot + ListSize(ptr) - 1;

  sz = tot - (long)RDB[DATA_PART_ALLOC_N] - (long)RDB[DATA_PART_ALLOC_G];

  /* Check count */

  if (sz != 0)
    Die(FUNCTION_NAME, "%ld particles lost", -sz);
}

/*****************************************************************************/
