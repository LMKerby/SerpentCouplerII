/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorepoison.c                                  */
/*                                                                           */
/* Created:       2012/12/04 (JLe)                                           */
/* Last modified: 2014/12/19 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Scores production and absorption rates for fission product   */
/*              poisoning                                                    */
/*                                                                           */
/* Comments: - Toi luuppi fissiolistan yli meinaa sitä että pelkästään       */
/*             alimman energian yieldi otetaan mukaan vaikka rea-rakenteessa */
/*             yieldit on kaikille. Serpent 1 toimii samalla tavalla.        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScorePoison:"

/*****************************************************************************/

void ScorePoison(double flx, long mat, double E, double wgt, double g, long id)
{
  long rea, ptr, ncol, gcu, ntot, ng, uni;
  double adens, xs, Ip, Xep, Pmp, Smp, Ia, Xea, Pma, Sma, Emax, Emin;
  double aI, aXe, aPm, aSm;

  /* Check poison calculation flag */

  if (((long)RDB[DATA_OPTI_POISON_CALC] == NO) &&
      ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] != YES) &&
      ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] != YES))
    return;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Check fissile flag */

  if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT))
    return;

  /***************************************************************************/

  /***** Calculate production rates ******************************************/

  /* Reset cross sections */

  Ip = 0.0;
  Xep = 0.0;
  Pmp = 0.0;
  Smp = 0.0;

  /* Reset densities */

  aI = 0.0;
  aXe = 0.0;
  aPm = 0.0;
  aSm = 0.0;

  /* Pointer to fission list */
  
  ptr = (long)RDB[mat + MATERIAL_PTR_FISS_REA_LIST];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Reset reaction pointer (rewind list) */
	  
  rea = -1;

  /* Loop over reactions */
	  
  while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
    {
      /* Check reaction pointer */
	      
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Get cross section */

      xs = adens*g*MicroXS(rea, E, id);
      
      /* Add to cross sections */

      Ip = Ip + xs*RDB[rea + REACTION_I135_YIELD];
      Xep = Xep + xs*RDB[rea + REACTION_XE135_YIELD];
      Pmp = Pmp + xs*RDB[rea + REACTION_PM149_YIELD];
      Smp = Smp + xs*RDB[rea + REACTION_SM149_YIELD];

      /* Check energy cut-off */
      
      if (E < Emin)
	break;
    }

  /* Get I-135 absorption cross section */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_I135_ISO]) > VALID_PTR)
    { 
      /* Get atomic density */

      aI = RDB[ptr + COMPOSITION_ADENS];

      /* Get pointer to nuclide data */

      ptr = (long)RDB[ptr + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      
      Ia = MicroXS(rea, E, id);

      /* Add capture to isomeric state (for the sake of completeness) */

      if ((rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS_ISO]) > VALID_PTR)
	Ia = Ia + MicroXS(rea, E, id);
    }
  else
    Ia = 0.0;

  /* Get Xe-135 absorption cross section */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_XE135_ISO]) > VALID_PTR)
    {
      /* Get atomic density */

      aXe = RDB[ptr + COMPOSITION_ADENS];

      /* Get pointer to nuclide data */

      ptr = (long)RDB[ptr + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      
      Xea = MicroXS(rea, E, id);

      /* Add capture to isomeric state (for the sake of completeness) */

      if ((rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS_ISO]) > VALID_PTR)
	Xea = Xea + MicroXS(rea, E, id);
    }
  else
    Xea = 0.0;

  /* Get Pm-149 absorption cross section */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_PM149_ISO]) > VALID_PTR)
    {
      /* Get atomic density */

      aPm = RDB[ptr + COMPOSITION_ADENS];

      /* Get pointer to nuclide data */

      ptr = (long)RDB[ptr + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      
      Pma = MicroXS(rea, E, id);

      /* Add capture to isomeric state (for the sake of completeness) */

      if ((rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS_ISO]) > VALID_PTR)
	Pma = Pma + MicroXS(rea, E, id);
    }
  else
    Pma = 0.0;

  /* Get Sm-149 absorption cross section */

  if ((ptr = (long)RDB[mat + MATERIAL_PTR_SM149_ISO]) > VALID_PTR)
    {
      /* Get atomic density */

      aSm = RDB[ptr + COMPOSITION_ADENS];

      /* Get pointer to nuclide data */

      ptr = (long)RDB[ptr + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      
      Sma = MicroXS(rea, E, id);

      /* Add capture to isomeric state (for the sake of completeness) */

      if ((rea = (long)RDB[ptr + NUCLIDE_PTR_NGAMMAXS_ISO]) > VALID_PTR)
	Sma = Sma + MicroXS(rea, E, id);
    }
  else
    Sma = 0.0;

  /***************************************************************************/

  /***** Material-wise rates *************************************************/

  /* Production rates */

  ptr = (long)RDB[mat + MATERIAL_PTR_I135_PROD_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Ip, wgt, ptr, id, 0);

  ptr = (long)RDB[mat + MATERIAL_PTR_XE135_PROD_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Xep, wgt, ptr, id, 0);

  ptr = (long)RDB[mat + MATERIAL_PTR_PM149_PROD_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Pmp, wgt, ptr, id, 0);

  ptr = (long)RDB[mat + MATERIAL_PTR_SM149_PROD_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Smp, wgt, ptr, id, 0);

  /* Microscopic absorption rates */

  ptr = (long)RDB[mat + MATERIAL_PTR_I135_ABS_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Ia, wgt, ptr, id, 0);

  ptr = (long)RDB[mat + MATERIAL_PTR_XE135_ABS_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Xea, wgt, ptr, id, 0);

  ptr = (long)RDB[mat + MATERIAL_PTR_SM149_ABS_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Sma, wgt, ptr, id, 0);

  ptr = (long)RDB[mat + MATERIAL_PTR_PM149_ABS_RATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx*Pma, wgt, ptr, id, 0);

  /***************************************************************************/

  /***** Parameters for group constant generation ****************************/

  /* Check poison calculation */
  
  if ((long)RDB[DATA_OPTI_POISON_CALC] == NO)
    return;

  /* Check that group constants are calculated */

  if ((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check if active cycle */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Check corrector calculation */

 if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    return;

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /***************************************************************************/

  /***** Group constant generation *******************************************/

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
  
      /***** Direct multi-group calculation **********************************/

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
	  
	  /* Production rates */
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_I135_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Ip, wgt, ptr, id, ng);
	  AddBuf1D(flx*Ip, wgt, ptr, id, 0);
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_XE135_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Xep, wgt, ptr, id, ng);
	  AddBuf1D(flx*Xep, wgt, ptr, id, 0);
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_PM149_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Pmp, wgt, ptr, id, ng);
	  AddBuf1D(flx*Pmp, wgt, ptr, id, 0);
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_SM149_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Smp, wgt, ptr, id, ng);
	  AddBuf1D(flx*Smp, wgt, ptr, id, 0);
	  
	  /* Microscopic absorption rates */
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_I135_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Ia, wgt, ptr, id, ng);
	  AddBuf1D(flx*Ia, wgt, ptr, id, 0);
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_XE135_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Xea, wgt, ptr, id, ng);
	  AddBuf1D(flx*Xea, wgt, ptr, id, 0);
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_SM149_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Sma, wgt, ptr, id, ng);
	  AddBuf1D(flx*Sma, wgt, ptr, id, 0);
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_PM149_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx*Pma, wgt, ptr, id, ng);
	  AddBuf1D(flx*Pma, wgt, ptr, id, 0);
	}

#endif
      
      /***********************************************************************/

      /***** Micro-group calculation *****************************************/

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
	  
	  /* Production rates */
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_I135_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Ip, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_XE135_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Xep, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM149_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Pmp, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_SM149_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Smp, id);
	  
	  /* Microscopic absorption rates */
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_I135_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Ia, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_XE135_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Xea, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM149_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Pma, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_SM149_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Sma, id);
	  
	  /* Macroscopic absorption rates */
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_I135_MACRO_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Ia*aI, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_XE135_MACRO_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Xea*aXe, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM149_MACRO_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Pma*aPm, id);
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_SM149_MACRO_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddPrivateRes(ptr + ng, wgt*flx*Sma*aSm, id);
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
