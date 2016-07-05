/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storesimdata.c                                 */
/*                                                                           */
/* Created:       2013/03/01 (JLe)                                           */
/* Last modified: 2014/08/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Stores group constant data for simulator output              */
/*                                                                           */
/* Comments: - Tää pitää fiksata kokonaan                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreSimData:"

/* Local function */

void StoreSimData0(long, long, long, long);

/*****************************************************************************/

void StoreSimData()
{

#ifdef SERPENT1_GC

  long gcu, nfg, ndg, ndfs, ndfc, loc0, ptr, adf, np;

  /* Check option */

  if ((long)RDB[DATA_SIMULATOR_DATA] == NO)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* Allow memory allocation */

  Mem(MEM_ALLOW);

  /* Loop over gc universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Create structure */

      loc0 = NewItem(DATA_PTR_SIM0, SIM_BLOCK_SIZE);

      /* Put gcu pointer */

      WDB[loc0 + SIM_PTR_GCU] = (double)gcu;

      /* Get number of energy groups */

      nfg = (long)RDB[DATA_ERG_FG_NG];
      CheckValue(FUNCTION_NAME, "nfg", "", nfg, 1, 100000);
      WDB[loc0 + SIM_NFG] = (double)nfg;

      /* Get pointer to few-group energy grid */
  
      ptr = (long)RDB[DATA_ERG_FG_PTR_GRID];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to groups */
      
      ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      WDB[loc0 + SIM_PTR_GROUPS] = (double)ptr;

      /* Get number of delayed neutron groups */

      ndg = (long)RDB[DATA_PRECURSOR_GROUPS];
      WDB[loc0 + SIM_NDG] = (double)ndg;

      /* Get pointer to ADF block */

      if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) > VALID_PTR)
	{
	  /* Get number of ADF's */
      
	  ndfs = (long)RDB[adf + ADF_NSURF];
	  WDB[loc0 + SIM_NDFS] = (double)ndfs;

	  ndfc = (long)RDB[adf + ADF_NCORN];
	  WDB[loc0 + SIM_NDFC] = (double)ndfc;

	  /* Pointer to surface */

	  ptr = (long)RDB[adf + ADF_PTR_SURF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Get surface type */

	  WDB[loc0 + SIM_ADF_SURF_TYPE] = RDB[ptr + SURFACE_TYPE];

	  /* Surface data */

	  if (ndfs > 0)
	    {
	      StoreSimData0(loc0 + SIM_ADFS, gcu + GCU_RES_FG_DF_SURF_DF, ndfs, nfg);
	      StoreSimData0(loc0 + SIM_ADF_NET_INCURR, gcu + GCU_RES_FG_DF_SURF_NET_CURR, 
			    ndfs, nfg);
	      StoreSimData0(loc0 + SIM_ADF_SURF_FLX, gcu + GCU_RES_FG_DF_HET_SURF_FLUX, 
			    ndfs, nfg);
	    }
      
	  /* Corner data */
      
	  if (ndfc > 0)
	    {
	      StoreSimData0(loc0 + SIM_ADF_CORN_FLX, gcu + GCU_RES_FG_DF_HET_CORN_FLUX, ndfs, nfg);
	      StoreSimData0(loc0 + SIM_ADFC, gcu + GCU_RES_FG_DF_CORN_DF, ndfs, nfg);
	    }

	  /* Volume flux */

	  StoreSimData0(loc0 + SIM_ADF_VOL_FLX, gcu + GCU_RES_FG_DF_HET_VOL_FLUX, nfg, 1); 
	}

      /* K-inf */
      
      if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == YES)
	StoreSimData0(loc0 + SIM_KINF, RES_IMP_KEFF, 1, 1);
      else
	StoreSimData0(loc0 + SIM_KINF, RES_COL_KEFF, 1, 1);

      /* Total power density */
      
      StoreSimData0(loc0 + SIM_POWDENS, RES_TOT_POWDENS, 1, 1);

      /* Neutron generation time */

      if ((np = (long)RDB[DATA_IFP_CHAIN_LENGTH]) > 0)
	StoreSimData0(loc0 + SIM_NEUTRON_GENTIME, RES_ADJ_IFP_GEN_TIME, 3, np);
      else
	StoreSimData0(loc0 + SIM_NEUTRON_GENTIME, RES_FWD_IMP_GEN_TIME, 1, 1);

      /* Delayed neutron parameters */

      if (ndg > 0)
	{
	  /* Beta-eff's */

	  if ((np = (long)RDB[DATA_IFP_CHAIN_LENGTH]) > 0)
	    StoreSimData0(loc0 + SIM_BETA_EFF, RES_ADJ_IFP_IMP_BETA_EFF, 
			  ndg + 1, np);
	  else
	    StoreSimData0(loc0 + SIM_BETA_EFF, RES_ADJ_MEULEKAMP_BETA_EFF, 
			  ndg + 1, 1);

	  /* Beta-zero's */

	  StoreSimData0(loc0 + SIM_BETA_ZERO, RES_FWD_ANA_BETA_ZERO, 
			ndg + 1, 1);

	  /* Lambda's */

	  if ((np = (long)RDB[DATA_IFP_CHAIN_LENGTH]) > 0)
	    StoreSimData0(loc0 + SIM_LAMBDA, RES_ADJ_IFP_IMP_LAMBDA, 
			  ndg + 1, np);
	  else
	    StoreSimData0(loc0 + SIM_LAMBDA, RES_ADJ_MEULEKAMP_LAMBDA, 
			  ndg + 1, 1);
	}

      /* Inverse velocity */
      
      StoreSimData0(loc0 + SIM_RECIPVEL, gcu + GCU_RES_FG_RECIPVEL, nfg + 1, 1);

      /* Few-group constants in infinite spectrum */
      
      StoreSimData0(loc0 + SIM_INF_FLX, gcu + GCU_INF_FLX,  nfg, 1);
      StoreSimData0(loc0 + SIM_INF_REP_TIME, gcu + GCU_INF_REP_TIME,  nfg, 1);
      StoreSimData0(loc0 + SIM_INF_PROMPT_LIFE, gcu + GCU_INF_PROMPT_LIFE,  nfg, 1);
      StoreSimData0(loc0 + SIM_INF_TOT, gcu + GCU_INF_TOT,  nfg, 1);
      StoreSimData0(loc0 + SIM_INF_CAPT, gcu + GCU_INF_CAPT,  nfg, 1);
      StoreSimData0(loc0 + SIM_INF_FISS, gcu + GCU_INF_FISS,  nfg, 1);
      StoreSimData0(loc0 + SIM_INF_ABS, gcu + GCU_INF_ABS,  nfg, 1);
      StoreSimData0(loc0 + SIM_INF_NSF, gcu + GCU_INF_NSF, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_KAPPA, gcu + GCU_INF_KAPPA, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_INVV, gcu + GCU_INF_INVV, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_NUBAR, gcu + GCU_INF_NUBAR, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_RABSXS, gcu + GCU_INF_RABSXS, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_REMXS, gcu + GCU_INF_REMXS, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_CHIT, gcu + GCU_INF_CHIT, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_CHIP, gcu + GCU_INF_CHIP, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_CHID, gcu + GCU_INF_CHID, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_S0, gcu + GCU_INF_S0, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SP0, gcu + GCU_INF_SP0, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_S1, gcu + GCU_INF_S1, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SP1, gcu + GCU_INF_SP1, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_S2, gcu + GCU_INF_S2, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SP2, gcu + GCU_INF_SP2, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_S3, gcu + GCU_INF_S3, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SP3, gcu + GCU_INF_SP3, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_S4, gcu + GCU_INF_S4, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SP4, gcu + GCU_INF_SP4, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_S5, gcu + GCU_INF_S5, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SP5, gcu + GCU_INF_SP5, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_S6, gcu + GCU_INF_S6, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SP6, gcu + GCU_INF_SP6, nfg, nfg);
      StoreSimData0(loc0 + SIM_INF_SCATT0, gcu + GCU_INF_SCATT0, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATT1, gcu + GCU_INF_SCATT1, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATT2, gcu + GCU_INF_SCATT2, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATT3, gcu + GCU_INF_SCATT3, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATT4, gcu + GCU_INF_SCATT4, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATT5, gcu + GCU_INF_SCATT5, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATT6, gcu + GCU_INF_SCATT6, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATT7, gcu + GCU_INF_SCATT7, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP0, gcu + GCU_INF_SCATTP0, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP1, gcu + GCU_INF_SCATTP1, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP2, gcu + GCU_INF_SCATTP2, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP3, gcu + GCU_INF_SCATTP3, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP4, gcu + GCU_INF_SCATTP4, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP5, gcu + GCU_INF_SCATTP5, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP6, gcu + GCU_INF_SCATTP6, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_SCATTP7, gcu + GCU_INF_SCATTP7, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_I135_YIELD, gcu + GCU_INF_I135_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_XE135_YIELD, gcu + GCU_INF_XE135_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_PM149_YIELD, gcu + GCU_INF_PM149_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_SM149_YIELD, gcu + GCU_INF_SM149_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_I135_ABS, gcu + GCU_INF_I135_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_XE135_ABS, gcu + GCU_INF_XE135_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_PM149_ABS, gcu + GCU_INF_PM149_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_SM149_ABS, gcu + GCU_INF_SM149_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_I135_MACRO_ABS, gcu + GCU_INF_I135_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_XE135_MACRO_ABS, gcu + GCU_INF_XE135_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_PM149_MACRO_ABS, gcu + GCU_INF_PM149_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_SM149_MACRO_ABS, gcu + GCU_INF_SM149_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_INF_TRANSPXS, gcu + GCU_INF_TRANSPXS, nfg, 1);
      StoreSimData0(loc0 + SIM_INF_DIFFCOEF, gcu + GCU_INF_DIFFCOEF, nfg, 1);

      /* Few-group constants in critical spectrum */

      StoreSimData0(loc0 + SIM_B1_FLX, gcu + GCU_B1_FLX,  nfg, 1);
      StoreSimData0(loc0 + SIM_B1_REP_TIME, gcu + GCU_B1_REP_TIME,  nfg, 1);
      StoreSimData0(loc0 + SIM_B1_PROMPT_LIFE, gcu + GCU_B1_PROMPT_LIFE,  nfg, 1);
      StoreSimData0(loc0 + SIM_B1_TOT, gcu + GCU_B1_TOT,  nfg, 1);
      StoreSimData0(loc0 + SIM_B1_CAPT, gcu + GCU_B1_CAPT,  nfg, 1);
      StoreSimData0(loc0 + SIM_B1_FISS, gcu + GCU_B1_FISS,  nfg, 1);
      StoreSimData0(loc0 + SIM_B1_ABS, gcu + GCU_B1_ABS,  nfg, 1);
      StoreSimData0(loc0 + SIM_B1_NSF, gcu + GCU_B1_NSF, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_KAPPA, gcu + GCU_B1_KAPPA, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_INVV, gcu + GCU_B1_INVV, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_NUBAR, gcu + GCU_B1_NUBAR, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_RABSXS, gcu + GCU_B1_RABSXS, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_REMXS, gcu + GCU_B1_REMXS, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_CHIT, gcu + GCU_B1_CHIT, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_CHIP, gcu + GCU_B1_CHIP, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_CHID, gcu + GCU_B1_CHID, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_S0, gcu + GCU_B1_S0, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SP0, gcu + GCU_B1_SP0, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_S1, gcu + GCU_B1_S1, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SP1, gcu + GCU_B1_SP1, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_S2, gcu + GCU_B1_S2, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SP2, gcu + GCU_B1_SP2, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_S3, gcu + GCU_B1_S3, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SP3, gcu + GCU_B1_SP3, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_S4, gcu + GCU_B1_S4, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SP4, gcu + GCU_B1_SP4, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_S5, gcu + GCU_B1_S5, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SP5, gcu + GCU_B1_SP5, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_S6, gcu + GCU_B1_S6, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SP6, gcu + GCU_B1_SP6, nfg, nfg);
      StoreSimData0(loc0 + SIM_B1_SCATT0, gcu + GCU_B1_SCATT0, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATT1, gcu + GCU_B1_SCATT1, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATT2, gcu + GCU_B1_SCATT2, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATT3, gcu + GCU_B1_SCATT3, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATT4, gcu + GCU_B1_SCATT4, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATT5, gcu + GCU_B1_SCATT5, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATT6, gcu + GCU_B1_SCATT6, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATT7, gcu + GCU_B1_SCATT7, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP0, gcu + GCU_B1_SCATTP0, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP1, gcu + GCU_B1_SCATTP1, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP2, gcu + GCU_B1_SCATTP2, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP3, gcu + GCU_B1_SCATTP3, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP4, gcu + GCU_B1_SCATTP4, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP5, gcu + GCU_B1_SCATTP5, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP6, gcu + GCU_B1_SCATTP6, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_SCATTP7, gcu + GCU_B1_SCATTP7, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_I135_YIELD, gcu + GCU_B1_I135_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_XE135_YIELD, gcu + GCU_B1_XE135_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_PM149_YIELD, gcu + GCU_B1_PM149_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_SM149_YIELD, gcu + GCU_B1_SM149_YIELD, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_I135_ABS, gcu + GCU_B1_I135_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_XE135_ABS, gcu + GCU_B1_XE135_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_PM149_ABS, gcu + GCU_B1_PM149_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_SM149_ABS, gcu + GCU_B1_SM149_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_I135_MACRO_ABS, gcu + GCU_B1_I135_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_XE135_MACRO_ABS, gcu + GCU_B1_XE135_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_PM149_MACRO_ABS, gcu + GCU_B1_PM149_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_SM149_MACRO_ABS, gcu + GCU_B1_SM149_MACRO_ABS, nfg, 1);	  
      StoreSimData0(loc0 + SIM_B1_TRANSPXS, gcu + GCU_B1_TRANSPXS, nfg, 1);
      StoreSimData0(loc0 + SIM_B1_DIFFCOEF, gcu + GCU_B1_DIFFCOEF, nfg, 1);
        
      /* Next universe */

      gcu = NextItem(gcu);
    }
  
  /* Disallow memory allocation */

  Mem(MEM_DENY);
}

/*****************************************************************************/

void StoreSimData0(long param, long gc, long sz1, long sz2)
{
  long loc0, ptr, n, m;

  /* Allocate memory for data */

  loc0 = ReallocMem(DATA_ARRAY, sz1*sz2);
  WDB[param] = (double)loc0;

  /* Pointer to data */
      
  if ((ptr = (long)RDB[gc]) < VALID_PTR)
    return;

  /* Check dimension and write data */

  if ((long)RDB[ptr + SCORE_DIM] == 1)
    {
      for (n = 0; n < sz1; n++)
	WDB[loc0++] = Mean(ptr, n);
    }
  else
    {
      for (n = 0; n < sz1; n++)
	for (m = 0; m < sz2; m++)
	  WDB[loc0++] = Mean(ptr, n, m);
    }

#endif
}

/*****************************************************************************/

