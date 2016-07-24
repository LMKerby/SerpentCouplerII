#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : flushbank.c                                    */
/*                                                                           */
/* Created:       2012/10/10 (JLe)                                           */
/* Last modified: 2014/02/12 (JLe)                                           */
/* Version:       2.1.17                                                      */
/*                                                                           */
/* Description: Removes all particles from bank                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FlushBank:"

/*****************************************************************************/

void FlushBank()
{
  long ptr, part, id;

  /* Get pointer to source */

  ptr = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Check that source is empty */

  if (ListSize(ptr) != 1)
    Die(FUNCTION_NAME, "Source is not empty");

  /* Loop over threads */

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Loop until bank is empty */

      while ((part = FromBank(id)) > VALID_PTR)
	{
	  /* Put particle back to stack */

	  ToStack(part, id);
	}
    }

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

  /* Re-distribute stacks (does now work in track plot mode) */
  
  if ((long)RDB[DATA_STOP_AFTER_PLOT] != STOP_AFTER_PLOT_TRACKS)
    ReDistributeStacks();
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
