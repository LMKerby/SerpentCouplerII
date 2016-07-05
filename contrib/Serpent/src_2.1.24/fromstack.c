/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fromstack.c                                    */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2014/10/08 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Retrieves neutron / photon from stack                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FromStack:"

/*****************************************************************************/

long FromStack(long type, long id)
{
  long ptr, sz;

#ifdef OLD_HIST

  long hst;

#endif

#ifdef OLD_IFP

  long prg;

#endif

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Avoid compiler warning */

  ptr = -1;

  /* Check type and get pointer to stack */

  if (type == PARTICLE_TYPE_NEUTRON)
    ptr = (long)RDB[OMPPtr(DATA_PART_PTR_NSTACK, id)];
  else if (type == PARTICLE_TYPE_GAMMA)
    ptr = (long)RDB[OMPPtr(DATA_PART_PTR_GSTACK, id)];
  else
    Die(FUNCTION_NAME, "Invalid particle type");
  
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get size */

  sz = ListSize(ptr) - 2;

  /* Check minimum level */
  
#ifdef OPEN_MP
#pragma omp critical
#endif
  {
    /* Compare sum to minimum level */
    
    if ((type == PARTICLE_TYPE_NEUTRON) &&
	(sz < (long)RDB[DATA_PART_MIN_NSTACK]))
      WDB[DATA_PART_MIN_NSTACK] = (double)sz;

    if ((type == PARTICLE_TYPE_GAMMA) &&
	(sz < (long)RDB[DATA_PART_MIN_GSTACK]))
      WDB[DATA_PART_MIN_GSTACK] = (double)sz;
  }

  /* Get pointer to last item */

  ptr = LastItem(ptr);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check type */
  
  if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
    {
      if (type == PARTICLE_TYPE_NEUTRON)
	Error(0, "Insufficient neutron buffer size, increase value of\nparameter \"nbuf\" (currently set to %1.1f)", RDB[DATA_PART_NBUF_FACTOR]);
      else
  	Error(0, "Insufficient photon buffer size, increase value of\nparameter \"gbuf\" (currently set to %1.1f)", RDB[DATA_PART_GBUF_FACTOR]);
    }

  /* Remove particle from stack */
  
  RemoveItem(ptr);

  /* Remember pointers to fission progeny and history data */

#ifdef OLD_IFP

  prg = (long)RDB[ptr + PARTICLE_PTR_FISS_PROG];

#endif

#ifdef OLD_HIST

  hst = (long)RDB[ptr + PARTICLE_PTR_HIST];

#endif
  
  /* Wipe data */

  memset(&WDB[ptr + LIST_DATA_SIZE], 0.0, 
	 (PARTICLE_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

  /* Put type */

  WDB[ptr + PARTICLE_TYPE] = (double)type;

  /* Put pointers */

  WDB[ptr + PARTICLE_PTR_EVENTS] = NULLPTR;

#ifdef OLD_IFP

  WDB[ptr + PARTICLE_PTR_FISS_PROG] = (double)prg;

#endif

#ifdef OLD_HIST

  WDB[ptr + PARTICLE_PTR_HIST] = (double)hst;

  /* Check history pointer */

  if (hst > VALID_PTR)
    {
      /* Loop over points */

      do
	{
	  /* Wipe data */
	  
	  memset(&WDB[hst + LIST_DATA_SIZE], 0.0, 
		 (HIST_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));
	  
	  /* Reset weight to indicate unused value */
		
	  WDB[hst + HIST_WGT] = -1.0;

	  /* Pointers to previous */

	  hst = PrevItem(hst);
	}
      while (hst != (long)RDB[ptr + PARTICLE_PTR_HIST]);
    }

#endif
  
  /* Loop over progeny data */

#ifdef OLD_IFP

  while (prg > VALID_PTR)
    {
      /* Reset data */
      
      WDB[prg + FISS_PROG_DN_GROUP] = 0.0;
      WDB[prg + FISS_PROG_LIFETIME] = 0.0; 
      WDB[prg + FISS_PROG_LAMBDA] = 0.0;
      
      /* Pointer to next */

      prg = NextItem(prg);
    }	  

#endif

  /* Reset ICM data */

  WDB[ptr + PARTICLE_ICM_PTR_ICM] = -1.0;
  WDB[ptr + PARTICLE_ICM_IDX] = -1.0;
  WDB[ptr + PARTICLE_ICM_MUA] = -1.0;
  WDB[ptr + PARTICLE_ICM_MUS] = -1.0;
  WDB[ptr + PARTICLE_ICM_G] = -1.0;
  WDB[ptr + PARTICLE_ICM_WGT] = -1.0;

  /* Reset albedo stuff */

  WDB[ptr + PARTICLE_ALB_PTR_GCU] = -1.0;
  WDB[ptr + PARTICLE_ALB_SURF_IDX] = -1.0;
  WDB[ptr + PARTICLE_ALB_G] = -1.0;
  
  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
