/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : bankstostore.c                                 */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Stores dynamic mode neutrons in the end of time interval     */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BanksToStore:"

/*****************************************************************************/

void BanksToStore()
{
  long ptr, new, id, np, ptr2, nb, i, tb, tt, idx;

 /* Loop over threads */

  /*printf("Storing neutrons to 3rd bank\n");*/

  /* Get pointer to particle counts */

  ptr2 = (long)RDB[DATA_PTR_DYN_PARTCOUNT];

  nb = (long)RDB[DATA_CYCLE_IDX];

  /* Get total number of batches */
  
  tb = (long)RDB[DATA_SRC_BATCHES];

  /* Get total number of threads */

  tt = (long)RDB[DATA_OMP_MAX_THREADS];

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
   {
     /* printf("np = %f \n", RDB[ptr2 + nb + id*(long)RDB[DATA_SRC_BATCHES] + (long)RDB[DATA_SRC_BATCHES]*(long)RDB[DATA_OMP_MAX_THREADS]]);*/

     /* Clear the previous neutrons from the EOI store */

     /* Calculate index for this thread and batch */

     idx =  nb + id*tb + tb*tt;

     /* Get current particle count in EOI store for this batch */
     
     np = (long)RDB[ptr2 + idx];

     for(i = 0; i < np; i++)
       {
	 ptr = FromStore(id, 1);

	 ToStack(ptr, id);     
       }

     /* Put particles from first bank to (the back of) the EOI store (and count them) */

     np = 0;

     while ((ptr = FromBank(id)) > VALID_PTR)
       {
	 /* Check type */
	 
	 if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
	   Die(FUNCTION_NAME, "Invalid particle type");
	 
	 new = DuplicateParticle(ptr, id);

	 /* Put new to EOI store */

	 ToStore(new, id, 1);

	 /* Put original to stack */
	     
	 ToStack(ptr,id);
	 
	 /* Count particle*/

	 np++;
       }
         
     /* Store number of particles in third bank for this thread */

     WDB[ptr2 + idx] = np;

     /*printf("%ld \n", np);*/

   }
  /*printf("\n");*/
}
