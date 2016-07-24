#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : movestore.c                                    */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Moves dynamic mode neutrons from EOI store to BOI store      */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveStore:"

/*****************************************************************************/

void MoveStore()
{
  long id, ptr, ptr2, new, np, np3, nb, i;
  long tb, tt, idxb, idxe;

  nb = (long)RDB[DATA_CYCLE_IDX];

  ptr2 = (long)RDB[DATA_PTR_DYN_PARTCOUNT];

  /* Get total number of batches */
  
  tb = (long)RDB[DATA_SRC_BATCHES];

  /* Get total number of threads */

  tt = (long)RDB[DATA_OMP_MAX_THREADS];

  /* Clear second banks and copy third banks to second banks */

  /* Loop over batches */

  for(nb = 0 ; nb < (long)RDB[DATA_SRC_BATCHES] ; nb++)
    {
      /* Loop over threads */

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	{
     
	  /* Calculate index for this thread and batch */

	  idxb =  nb + id*tb;

	  idxe =  nb + id*tb + tb*tt;

	  /* Clear the previous neutrons from the BOI store */

	  np = RDB[ptr2 + idxb];

	  for(i = 0; i < np; i++)
	    {
	      /* Get particle from BOI store */

	      ptr = FromStore(id, 0);

	      /* Put particle to stack */

	      ToStack(ptr, id);

	    }

	  /* Put this batches particles from EOI store to BOI store (and count them) */

	  np = 0;

	  np3 = RDB[ptr2 + idxe];

	  for(i = 0; i < np3; i++)
	    {
	      /* Particle from EOI store */

	      ptr = FromStore(id, 1);

	      /* Check type */
	 
	      if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
		Die(FUNCTION_NAME, "Invalid particle type");

	      /* Duplicate from second stack */
	 
	      new = DuplicateParticle(ptr, id);

	      /* Put duplicate to BOI store */

	      ToStore(new, id, 0);

	      /* Put original to stack */

	      ToStack(ptr, id);

	    }

	  /* Store number of particles in second bank for this thread */

	  WDB[ptr2 + idxb] = np3;
	  /*printf("%ld ", np3);*/

	  /* Third bank should not contain particles for this batch */

	  WDB[ptr2 + idxe] = 0.0;

	}
      /* printf("\n");*/
    }

}
#ifdef __cplusplus 
} 
#endif 
