/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : duplicateparticle.c                            */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2013/10/08 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Duplicates particle                                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DuplicateParticle:"

/*****************************************************************************/

long DuplicateParticle(long ptr, long id)
{
  long new;

#ifdef OLD_HIST

  long hst0, hst1;

#endif

#ifdef OLD_IFP

  long prg0, prg1;

#endif

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get particle from stack */

  new = FromStack((long)RDB[ptr + PARTICLE_TYPE], id);

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(new)", DATA_ARRAY, new);

  /* Get history pointers */

#ifdef OLD_HIST

  hst0 = (long)RDB[ptr + PARTICLE_PTR_HIST];
  hst1 = (long)RDB[new + PARTICLE_PTR_HIST];

#endif

  /* Get fission progeny pointers */

#ifdef OLD_IFP

  prg0 = (long)RDB[ptr + PARTICLE_PTR_FISS_PROG];
  prg1 = (long)RDB[new + PARTICLE_PTR_FISS_PROG];

#endif

  /* Copy data */

  memcpy(&WDB[new + LIST_DATA_SIZE], &RDB[ptr + LIST_DATA_SIZE],
	 (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

  /* Restore pointers */

#ifdef OLD_HIST

  WDB[new + PARTICLE_PTR_HIST] = (double)hst1;

#endif

#ifdef OLD_IFP

  WDB[new + PARTICLE_PTR_FISS_PROG] = (double)prg1;

#endif

  /* Copy fission progeny data */

#ifdef OLD_IFP

  while (prg0 > VALID_PTR)
    {
      /* Check second pointer */
      
      CheckPointer(FUNCTION_NAME, "(prg1)", DATA_ARRAY, prg1);

      /* Copy data */

      memcpy(&WDB[prg1 + LIST_DATA_SIZE], &RDB[prg0 + LIST_DATA_SIZE],
	     (FISS_PROG_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

      /* Pointers to next */

      prg0 = NextItem(prg0);
      prg1 = NextItem(prg1);
    }

#endif
  
  /* Check history pointer */

#ifdef OLD_HIST

  if (hst0 > VALID_PTR)
    {
      /* Check second pointer */
      
      CheckPointer(FUNCTION_NAME, "(hst1)", DATA_ARRAY, hst1);

      /* Loop over data */

      do
	{
	  /* Check second pointer */

	  CheckPointer(FUNCTION_NAME, "(hst1)", DATA_ARRAY, hst1);

	  /* Copy data */

	  memcpy(&WDB[hst1 + LIST_DATA_SIZE], &RDB[hst0 + LIST_DATA_SIZE],
		 (HIST_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

	  /* Pointers to previous */

	  hst0 = PrevItem(hst0);
	  hst1 = PrevItem(hst1);
	}
      while (hst0 != (long)RDB[ptr + PARTICLE_PTR_HIST]);
    }

#endif

#ifdef OPEN_MP
#pragma omp critical
#endif
  {
    /* Loop over events and add to counters */

    ptr = (long)RDB[new + PARTICLE_PTR_EVENTS];
    while (ptr > VALID_PTR)
      {
	/* Add to count */

	if ((long)RDB[ptr + EVENT_HIS_COUNT] < 1)
	  Die(FUNCTION_NAME, "WTF? %ld", (long)RDB[ptr + EVENT_HIS_COUNT]);
	else
	  WDB[ptr + EVENT_HIS_COUNT] = RDB[ptr + EVENT_HIS_COUNT] + 1.0;
	
	/* Next */
	
	ptr = NextItem(ptr);
      }
  }

  /* Return pointer */
  
  return new;
}

/*****************************************************************************/
