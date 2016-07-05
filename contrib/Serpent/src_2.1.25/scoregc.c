/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoregc.c                                      */
/*                                                                           */
/* Created:       2011/05/04 (JLe)                                           */
/* Last modified: 2014/11/14 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Scores reaction rates needed for group constant generation   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreGC:"

/*****************************************************************************/

void ScoreGC(double flx, double tot, double capt, double fiss, double fissE,
	     double nsf, double ela, double sprod, long mat, long part, 
	     double E, double spd, double wgt, double g, long id)
{
  long gcu, ptr, loc0, ng, ntot, ncol, uni, dng;
  double val, lambda;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /***************************************************************************/

  /***** Meulekamp estimate of beta-eff **************************************/

  /* Get pointer to universe */

  if ((gcu = (long)RDB[part + PARTICLE_PTR_GCU]) > VALID_PTR)
    {
      /* Delayed neutrons */
      
      if ((dng = (long)RDB[part + PARTICLE_DN_GROUP]) > 0)
	{
	  /* Score fission rate */
	  
	  ptr = (long)RDB[gcu + GCU_MEULEKAMP_BETA_EFF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  AddBuf1D(fiss, wgt, ptr, id, 0);
	  AddBuf1D(fiss, wgt, ptr, id, dng);
	  
	  /* Get decay constant */
	  
	  lambda = RDB[part + PARTICLE_DN_LAMBDA];
	  
	  /* Score lambda */
	  
	  ptr = (long)RDB[gcu + GCU_MEULEKAMP_LAMBDA];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  AddBuf1D(fiss*lambda, wgt, ptr, id, 0);
	  AddBuf1D(fiss*lambda, wgt, ptr, id, dng);
	}
      
      /* Score total fission rate */
      
      ptr = (long)RDB[gcu + GCU_MEULEKAMP_TOT_FISS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(fiss, wgt, ptr, id, 0);
    }

  /***************************************************************************/

  /***************************************************************************/

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Check for multiple levels */

  if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
    {
      /* Single level, get pointer */

      if ((gcu = TestValuePair(DATA_GCU_PTR_UNI, ncol, id)) < VALID_PTR)
	return;
    }
  else
    {
      /* Multiple levels, get pointer to list */

      gcu = (long)RDB[DATA_PTR_GCU0];
      CheckPointer(FUNCTION_NAME, "(gcu)", DATA_ARRAY, gcu);
    }

  /****************************************************************************/

  /***** Score data **********************************************************/

  /* Loop over universes */
  
  while (gcu > VALID_PTR)
    {
      /* Check multi-level mode */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == YES)
	{
	  /* Pointer to universe */

	  uni = (long)RDB[gcu + GCU_PTR_UNIV];
	  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

	  /* Check collision */

	  if (TestValuePair(uni + UNIVERSE_COL_COUNT, ncol, id) < 0.0)
	    {
	      /* Next universe */
	      
	      gcu = NextItem(gcu);
	      
	      /* Cycle loop */
	      
	      continue;
	    }
	}
      
      /***********************************************************************/

      /***** MORA data *******************************************************/

      if ((loc0 = (long)RDB[gcu + GCU_PTR_MORA]) > VALID_PTR)
	{
	  /* Get pointer to energy grid */
	  
	  ptr = (long)RDB[loc0 + MORA_PTR_EG];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  /* Total number of groups */
	  
	  ntot = (long)RDB[loc0 + MORA_N_EG];
	  
	  /* Find group */
	  
	  if ((ng = GridSearch(ptr, E)) > -1)
	    {
	      /* Flux */
	      
	      ptr = (long)RDB[loc0 + MORA_PTR_FLX];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(flx, wgt, ptr, id, ng);
	      AddBuf1D(flx, wgt, ptr, id, ntot);
	      
	      /* Score total reaction rate */
	      
	      ptr = (long)RDB[loc0 + MORA_PTR_TOT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(tot, wgt, ptr, id, ng);
	      
	      /* Score fission rate */
	      
	      ptr = (long)RDB[loc0 + MORA_PTR_FISS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(fiss, wgt, ptr, id, ng);
	      
	      /* Score capture rate */
	      
	      ptr = (long)RDB[loc0 + MORA_PTR_CAPT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(capt, wgt, ptr, id, ng);
	      
	      /* Score fission energy production rate */
	      
	      ptr = (long)RDB[loc0 + MORA_PTR_KAPPA];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(fissE, wgt, ptr, id, ng);
	    }
	}
      
      /************************************************************************/

      /***** Infinite spectrum ***********************************************/

#ifdef SERPENT1_GC

      /* Number of groups */
      
      ntot = (long)RDB[DATA_ERG_FG_NG];
      
      /* Get pointer to few-group structure */
      
      ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
    
      /* Get few-group index */
      
      if ((ng = GridSearch(ptr, E)) > -1)
	{
	  /* Convert index */
	  
	  ng = ntot - ng;
	  CheckValue(FUNCTION_NAME, "ng1", "", ng, 0, ntot);

	  /* Flux */
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx, wgt, ptr, id, ng);
	  AddBuf1D(flx, wgt, ptr, id, 0);
	  
	  /* 1/v */
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_RECIPVEL];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx/spd, wgt, ptr, id, ng);
	  AddBuf1D(flx/spd, wgt, ptr, id, 0);
	  
	  /* Check material pointer */
	  
	  if (mat > VALID_PTR)
	    {  
	      /* Score total reaction rate */
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_TOTXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(tot, wgt, ptr, id, ng);
	      AddBuf1D(tot, wgt, ptr, id, 0);
	  
	      /* Score fission rate */
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_FISSXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(fiss, wgt, ptr, id, ng);
	      AddBuf1D(fiss, wgt, ptr, id, 0);
	      
	      /* Score fission neutron production rate */
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_NSF];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(nsf, wgt, ptr, id, ng);
	      AddBuf1D(nsf, wgt, ptr, id, 0);

	      /* Score fission energy production rate */
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_FISSE];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(fissE, wgt, ptr, id, ng);
	      AddBuf1D(fissE, wgt, ptr, id, 0);
	      
	      /* Score capture rate */
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_CAPTXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(capt, wgt, ptr, id, ng);
	      AddBuf1D(capt, wgt, ptr, id, 0);
	      
	      /* Score absorption rate */
	      
	      val = fiss + capt;
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_ABSXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(val, wgt, ptr, id, ng);
	      AddBuf1D(val, wgt, ptr, id, 0);
	      
	      /* Score reduced absorption rate */
	      
	      val = fiss + capt - sprod;
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_RABSXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(val, wgt, ptr, id, ng);
	      AddBuf1D(val, wgt, ptr, id, 0);
	      
	      /* Score elastic scattering rate */
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_ELAXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(ela, wgt, ptr, id, ng);
	      AddBuf1D(ela, wgt, ptr, id, 0);
	      
	      /* Score inelastic scattering rate (avoid round-off to */
	      /* negative) */
	      
	      if ((val = tot - capt - fiss - ela) < 0.0)
		val = 0.0;
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_INLXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(val, wgt, ptr, id, ng);
	      AddBuf1D(val, wgt, ptr, id, 0);
	      
	      /* Score (total) scattering rate */
	      
	      val = tot - capt - fiss;
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_SCATTXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(val, wgt, ptr, id, ng);
	      AddBuf1D(val, wgt, ptr, id, 0);
	      
	      /* Score scattering production rate */
	      
	      val = tot - capt - fiss + sprod;
	      
	      ptr = (long)RDB[gcu + GCU_RES_FG_SCATTPRODXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf1D(val, wgt, ptr, id, ng);
	      AddBuf1D(val, wgt, ptr, id, 0);
	    }
	}
      
#else
      /* Avoid compiler warning */

      val = 0.0;

#endif

      /***********************************************************************/

      /***** Micro-group data ************************************************/

      /* Get pointer to microgroup energy grid */
      
      ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	
      
      /* Number of groups */
      
      ntot = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
      
      /* Get group index */
      
      if ((ng = GridSearch(ptr, E)) > -1)
	{
	  /* Convert index */
	  
	  ng = ntot - ng - 1;
	  CheckValue(FUNCTION_NAME, "ng2", "", ng, 0, ntot - 1);
      
	  /* Put values */
	  
	  ptr = RDB[gcu + GCU_MICRO_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx, id);
	  
	  ptr = RDB[gcu + GCU_MICRO_TOT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*tot, id);
	  
	  ptr = RDB[gcu + GCU_MICRO_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*(capt + fiss), id);
	  
	  ptr = RDB[gcu + GCU_MICRO_FISS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*fiss, id);
	  
	  ptr = RDB[gcu + GCU_MICRO_FISSE];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*fissE, id);

	  ptr = RDB[gcu + GCU_MICRO_NSF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*nsf, id);
	  
	  ptr = RDB[gcu + GCU_MICRO_INV_V];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx/spd, id);
	  
	  /* Flux in fissile materials */
	  
	  if (mat > VALID_PTR)
	    if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT)
	      {
		ptr = RDB[gcu + GCU_MICRO_FISS_FLX];
		CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
		AddPrivateRes(ptr + ng, wgt*flx, id);
	      }
	}
      
      /***********************************************************************/

      /* Next universe */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == NO)
	break;
      else
	gcu = NextItem(gcu);
    }

  /***************************************************************************/
}

/*****************************************************************************/
