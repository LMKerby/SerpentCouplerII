#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iteratekeff.c                                  */
/*                                                                           */
/* Created:       2013/10/02 (JLe)                                           */
/* Last modified: 2014/04/02 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: K-eff iteration                                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IterateKeff:"

/*****************************************************************************/

void IterateKeff()
{
  long ptr, mode, skip, ncyc, idx, fix;
  double nsf, fiss, capt, nuxn, leak, L0, L1, val, f, keff, val0;
  
  /* Check mode */

  if ((mode = (long)RDB[DATA_ITER_MODE]) == ITER_MODE_NONE)
    return;

  /* Get fix mode */

  fix = (long)RDB[DATA_ITER_FIX];

  /* Number of cycles and actual number of skip cycles (setoptimization.c) */

  ncyc = (long)RDB[DATA_ITER_NCYC];
  idx = (long)RDB[DATA_CYCLE_IDX];

  if (fix == YES)
    skip = (long)((RDB[DATA_CRIT_SKIP] - RDB[DATA_ITER_NCYC])/2.0);
  else
    skip = (long)(RDB[DATA_CRIT_SKIP] - RDB[DATA_ITER_NCYC]);

  /* Check cycles */

  if ((idx < skip) || ((fix == YES) && (idx > skip + ncyc)))
    return;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Collect MPI parallel data */

  CollectBuf();

  /* Check mode */

  if (mode == ITER_MODE_ALBEDO)
    {
      /***********************************************************************/

      /***** Albedo iteration ************************************************/

      /* Get k-eff */

      keff = RDB[DATA_ITER_KEFF];
      CheckValue(FUNCTION_NAME, "keff", "", keff, 0.1, 2.5);

      /* Fission nubar */
      
      ptr = (long)RDB[RES_TOT_NSF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      nsf = BufVal(ptr, 0);
      
      /* Fission term */
      
      ptr = (long)RDB[RES_TOT_FISSRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      fiss = BufVal(ptr, 0);
      
      /* Total capture rate */
      
      ptr = (long)RDB[RES_TOT_CAPTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      capt = BufVal(ptr, 0);
      
      /* Scattering production rate */
      
      ptr = (long)RDB[RES_TOT_INLPRODRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      nuxn = BufVal(ptr, 0);

      /* Physical leakage rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      leak = BufVal(ptr, 0);

      /* Get previous albedo leakage rate */
      
      ptr = (long)RDB[RES_ALB_NEUTRON_LEAKRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      L0 = BufVal(ptr, 0);
     
      /* Calculate estimate for new albedo leakage rate */
      
      L1 = nsf/keff - capt - fiss - leak + nuxn;

      /* Avoid compiler warning */

      val = -1.0;

      /* Get previous value */

      if ((val0 = RDB[DATA_ITER_VAL]) < 0.0)
	{
	  /* Not set, use initial guess */

	  ptr = (long)RDB[RES_ANA_KEFF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  if (Mean(ptr, 0) > 1.0)
	    val = 0.999;
	  else
	    val = 1.001;
	}
      else if (L0 != 0.0)
	{
	  /* Calculate new */

	  val = (val0 - 1.0)*L1/L0 + 1.0;

	  /* Add to statistics */

	  ptr = (long)RDB[RES_ITER_VAL];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddStat(val, ptr, 0);
	  
	  /* Fix value for last iteration */
	  
	  if ((fix == YES) && (idx > skip + ncyc))
	    val = Mean(ptr, 0);
	}    
      else
	Die(FUNCTION_NAME, "L0 == 0");

      /* Put value */

      WDB[DATA_ITER_VAL] = val;

      /* Put albedos */

      if ((f = RDB[DATA_ITER_ALB_F1]) > 0.0)
	WDB[DATA_GEOM_ALBEDO1] = val*f;

      if ((f = RDB[DATA_ITER_ALB_F2]) > 0.0)
	WDB[DATA_GEOM_ALBEDO2] = val*f;
      
      if ((f = RDB[DATA_ITER_ALB_F3]) > 0.0)
	WDB[DATA_GEOM_ALBEDO3] = val*f;
      
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid iteration mode");
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
