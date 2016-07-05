/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : leak.c                                         */
/*                                                                           */
/* Created:       2011/03/12 (JLe)                                           */
/* Last modified: 2014/04/04 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Handles leakage                                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Leak:"

/*****************************************************************************/

void Leak(long part, double x, double y, double z, double u, double v, 
	  double w, double E, double wgt, long id)
{
  long gcu, ptr, ntot, ng, ncol, type, n0, ma0, ms0, ng0, wgt0, icm;

  /* Check pointer */
      
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Put particle back in stack */

  ToStack(part, id);

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];

  /* Score total leak rate */

  if (type == PARTICLE_TYPE_NEUTRON)
    ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
  else
    ptr = (long)RDB[RES_TOT_PHOTON_LEAKRATE];

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, wgt, ptr, id, 0);

  /* Check active cycle and corrector step */

  if ((RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP]) ||
      (type == PARTICLE_TYPE_GAMMA) ||
      ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
    return;  

  /***************************************************************************/
  
  /***** Group constant generation *******************************************/

#ifdef SERPENT1_GC

  /* Check for group constant calculation */

  if ((long)RDB[DATA_OPTI_GC_CALC] == YES)
    {
      /* Number of groups */

      ntot = (long)RDB[DATA_ERG_FG_NG];

      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      ncol = (long)GetPrivateData(ptr, id);
      
      /* Get pointer to GCU structure (cannot be multiple) */
      
      gcu = TestValuePair(DATA_GCU_PTR_UNI, ncol, id);
      
      /* Get pointer to few-group structure */
  
      ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get few-group and universe indexes */

      if ((gcu > VALID_PTR) && ((ng = GridSearch(ptr, E)) > -1))
	{
	  /* Score few-group leakage */

	  ptr = (long)RDB[gcu + GCU_RES_FG_LEAK];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(1.0, wgt, ptr, id, 0);
	  AddBuf1D(1.0, wgt, ptr, id, ntot - ng);
	}
    }

#else

  /* Avoid compiler warning */

  ncol = -1;
  ng = -1;
  ntot = -1;
  gcu = -1;

#endif

  /***************************************************************************/
  
  /***** ICM leak rate *******************************************************/

  /* Get universe pointer */
  
  if ((long)RDB[DATA_ICM_CALC] == YES)
    if ((icm = (long)RDB[part + PARTICLE_ICM_PTR_ICM]) > VALID_PTR)
      {
	/* Get data (toi nollallinen muuttuja on jo käytössä) */
	
	n0 = (long)RDB[part + PARTICLE_ICM_IDX];
	ma0 = (long)RDB[part + PARTICLE_ICM_MUA];
	ms0 = (long)RDB[part + PARTICLE_ICM_MUS];
	ng0 = (long)RDB[part + PARTICLE_ICM_G];
	wgt0 = RDB[part + PARTICLE_ICM_WGT];
	      
	/* Check */
	
	if (n0 > -1)
	  {
	    ptr = (long)RDB[icm + ICM_RES_LEAK1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    AddBuf(1.0, wgt, ptr, id, -1, n0, ma0, ms0, ng0);
	    
	    ptr = (long)RDB[icm + ICM_RES_LEAK2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    AddBuf(wgt0, wgt, ptr, id, -1, n0, ma0, ms0, ng0);
	  }
      }

  /***************************************************************************/
}

/*****************************************************************************/
