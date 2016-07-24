#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processstats.c                                 */
/*                                                                           */
/* Created:       2011/03/10 (JLe)                                           */
/* Last modified: 2015/03/04 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Allocates memory for statistical variables                   */
/*                                                                           */
/* Comments: - Statistics for detectors is handled in processdetectors.c     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessStats:"

/*****************************************************************************/

void ProcessStats()
{
  long loc0, loc1, ptr, np;

  /***************************************************************************/

  /***** Neutron integral parameters *****************************************/

  /* Reaction rates with single bin */

  ptr = NewStat("TOT_INLPRODRATE", 1, 1);
  WDB[RES_TOT_INLPRODRATE] = (double)ptr;  

  ptr = NewStat("TOT_NEUTRON_LEAKRATE", 1, 1);
  WDB[RES_TOT_NEUTRON_LEAKRATE] = (double)ptr;  

  ptr = NewStat("TOT_NEUTRON_CUTRATE", 1, 1);
  WDB[RES_TOT_NEUTRON_CUTRATE] = (double)ptr;  

  ptr = NewStat("ALB_NEUTRON_LEAKRATE", 1, 1);
  WDB[RES_ALB_NEUTRON_LEAKRATE] = (double)ptr;  

  ptr = NewStat("TOT_ELARATE", 1, 1);
  WDB[RES_TOT_ELARATE] = (double)ptr;  

  ptr = NewStat("TOT_NEUTRON_RR", 1, 1);
  WDB[RES_TOT_NEUTRON_RR] = (double)ptr;  
  
  /* Three bins for burnup normalization */

  ptr = NewStat("TOT_NSF", 1, 3);
  WDB[RES_TOT_NSF] = (double)ptr;  

  ptr = NewStat("TOT_NUBAR", 1, 3);
  WDB[RES_TOT_NUBAR] = (double)ptr;  

  ptr = NewStat("TOT_FISSE", 1, 3);
  WDB[RES_TOT_FISSE] = (double)ptr;  

  ptr = NewStat("TOT_FISSRATE", 1, 3);
  WDB[RES_TOT_FISSRATE] = (double)ptr;  

  ptr = NewStat("TOT_CAPTRATE", 1, 3);
  WDB[RES_TOT_CAPTRATE] = (double)ptr;  

  ptr = NewStat("TOT_NEUTRON_SRCRATE", 1, 3);
  WDB[RES_TOT_NEUTRON_SRCRATE] = (double)ptr;  

  ptr = NewStat("TOT_NEUTRON_FLUX", 1, 3);
  WDB[RES_TOT_NEUTRON_FLUX] = (double)ptr;  

  ptr = NewStat("RES_TOT_NEUTRON_POWER", 1, 3);
  WDB[RES_TOT_NEUTRON_POWER] = (double)ptr;

  ptr = NewStat("RES_TOT_POWDENS", 1, 3);
  WDB[RES_TOT_POWDENS] = (double)ptr;

  ptr = NewStat("RES_TOT_GENRATE", 1, 3);
  WDB[RES_TOT_GENRATE] = (double)ptr;

  ptr = NewStat("TOT_ABSRATE", 1, 3);
  WDB[RES_TOT_ABSRATE] = (double)ptr;  

  ptr = NewStat("TOT_NEUTRON_LOSSRATE", 1, 3);
  WDB[RES_TOT_NEUTRON_LOSSRATE] = (double)ptr;  

  /* Source weights */

  ptr = NewStat("INI_SRC_WGT", 1, 1);
  WDB[RES_INI_SRC_WGT] = (double)ptr;  

  ptr = NewStat("NEW_SRC_WGT", 1, MAX_EXT_K_GEN + 1);
  WDB[RES_NEW_SRC_WGT] = (double)ptr;  

  /* Source keff */

  ptr = NewStat("EXT_K", 1, MAX_EXT_K_GEN);
  WDB[RES_EXT_K] = (double)ptr;  

  /* Total 1/v */

  ptr = NewStat("TOT_RECIPVEL", 1, 1);
  WDB[RES_TOT_RECIPVEL] = (double)ptr;  

  /* Normalization coefficient */

  ptr = NewStat("NORM_COEF", 1, 2);
  WDB[RES_NORM_COEF] = (double)ptr;  

  /* Implicit estimator of k-eff and k-inf */

  ptr = NewStat("IMP_KEFF", 1, 1);
  AllocStatHistory(ptr);
  WDB[RES_IMP_KEFF] = (double)ptr;  

  ptr = NewStat("IMP_KINF", 1, 1);
  WDB[RES_IMP_KINF] = (double)ptr;  

  /* Analog estimator of k-eff */

  ptr = NewStat("ANA_KEFF", 1, 3);
  AllocStatHistory(ptr);
  WDB[RES_ANA_KEFF] = (double)ptr;  

  /* Source multiplication */

  ptr = NewStat("SRC_MULT", 1, 1);
  WDB[RES_SRC_MULT] = (double)ptr;  

  /* Number of generations */

  ptr = NewStat("SRC_MEAN_NGEN", 1, 1);
  WDB[RES_MEAN_NGEN] = (double)ptr;  

  /* Prompt fission generation fractions */

  if ((np = (long)RDB[DATA_MAX_PROMPT_CHAIN_LENGTH]) > 0)
    {
      ptr = NewStat("PROMPT_GEN_CUMU", 1, np + 1);
      WDB[RES_PROMPT_GEN_CUMU] = (double)ptr;  

      ptr = NewStat("PROMPT_GEN_POP", 1, np + 1);
      WDB[RES_PROMPT_GEN_POP] = (double)ptr;  

      ptr = NewStat("PROMPT_GEN_TIMES", 1, np + 1);
      WDB[RES_PROMPT_GEN_TIMES] = (double)ptr;  
    }

  ptr = NewStat("PROMPT_CHAIN_LENGTH", 1, 1);
  WDB[RES_PROMPT_CHAIN_LENGTH] = (double)ptr;  

  /* Wielandt method */
  
  if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
    {
      /* K-eff */

      ptr = NewStat("ANA_WIELANDT_K", 1, 1);
      AllocStatHistory(ptr);
      WDB[RES_WIELANDT_K] = (double)ptr;
      
      /* Probability */
      
      ptr = NewStat("COL_WIELANDT_P", 1, 1);
      AllocStatHistory(ptr);
      WDB[RES_WIELANDT_P] = (double)ptr;
    }  

  /* Collision estimator of k-eff */

  ptr = NewStat("COL_KEFF", 1, 1);
  AllocStatHistory(ptr);
  WDB[RES_COL_KEFF] = (double)ptr;  

  /* Iteration eigenvalue */

  ptr = NewStat("ITER_VAL", 1, 1);
  WDB[RES_ITER_VAL] = (double)ptr;  

  /* Albedo */

  ptr = NewStat("GEOM_ALBEDO", 1, 3);
  WDB[RES_GEOM_ALBEDO] = (double)ptr;  

  /* Conversion ratio */

  ptr = NewStat("CONV_RATIO", 1, 3);
  WDB[RES_ANA_CONV_RATIO] = (double)ptr;  

  /* Analog fission and capture fractions */

  ptr = NewStat("FISS_FRAC", 2, 8, 2);
  WDB[RES_ANA_FISS_FRAC] = (double)ptr;  

  ptr = NewStat("CAPT_FRAC", 2, 10, 2);
  WDB[RES_ANA_CAPT_FRAC] = (double)ptr;  

  /***************************************************************************/

  /***** Time constants ******************************************************/

  /* Mean number of collisions per history (total and to fission) */

  ptr = NewStat("ANA_MEAN_NCOL", 1, 2);
  WDB[RES_ANA_MEAN_NCOL] = (double)ptr;  

  /* Forward-weighted time constants */

  ptr = NewStat("FWD_ANA_BETA_ZERO", 1, 9);
  WDB[RES_FWD_ANA_BETA_ZERO] = (double)ptr;  

  ptr = NewStat("FWD_ANA_LAMBDA", 1, 9);
  WDB[RES_FWD_ANA_LAMBDA] = (double)ptr;  

  /* Adjoint-weighted time constants */

  ptr = NewStat("ADJ_MEULEKAMP_BETA_EFF", 1, 9);
  WDB[RES_ADJ_MEULEKAMP_BETA_EFF] = (double)ptr;  

  ptr = NewStat("ADJ_MEULEKAMP_LAMBDA", 1, 9);
  WDB[RES_ADJ_MEULEKAMP_LAMBDA] = (double)ptr;  
 
  ptr = NewStat("ADJ_NAUCHI_GEN_TIME", 1, 3);
  WDB[RES_ADJ_NAUCHI_GEN_TIME] = (double)ptr;  
  
  ptr = NewStat("ADJ_NAUCHI_LIFETIME", 1, 3);
  WDB[RES_ADJ_NAUCHI_LIFETIME] = (double)ptr;  

  ptr = NewStat("ADJ_NAUCHI_BETA_EFF", 1, 9);
  WDB[RES_ADJ_NAUCHI_BETA_EFF] = (double)ptr;  

  ptr = NewStat("ADJ_NAUCHI_LAMBDA", 1, 9);
  WDB[RES_ADJ_NAUCHI_LAMBDA] = (double)ptr;  

  /* Get IFP chain length */

  if ((np = (long)RDB[DATA_IFP_CHAIN_LENGTH]) > 0)
    {
      ptr = NewStat("ADJ_IFP_GEN_TIME", 2, 3, np);
      WDB[RES_ADJ_IFP_GEN_TIME] = (double)ptr;  

      ptr = NewStat("ADJ_IFP_LIFETIME", 2, 3, np);
      WDB[RES_ADJ_IFP_LIFETIME] = (double)ptr;  

      ptr = NewStat("ADJ_IFP_IMP_BETA_EFF", 2, 9, np);
      WDB[RES_ADJ_IFP_IMP_BETA_EFF] = (double)ptr;  

      ptr = NewStat("ADJ_IFP_IMP_LAMBDA", 2, 9, np);
      WDB[RES_ADJ_IFP_IMP_LAMBDA] = (double)ptr;  

      ptr = NewStat("ADJ_IFP_ANA_BETA_EFF", 2, 10, np);
      WDB[RES_ADJ_IFP_ANA_BETA_EFF] = (double)ptr;  

      ptr = NewStat("ADJ_IFP_ANA_LAMBDA", 2, 9, np);
      WDB[RES_ADJ_IFP_ANA_LAMBDA] = (double)ptr;  

      ptr = NewStat("ADJ_IFP_ROSSI_ALPHA", 1, np);
      WDB[RES_ADJ_IFP_ROSSI_ALPHA] = (double)ptr;

      ptr = NewStat("ADJ_PERT_GEN_TIME", 1, 3);
      WDB[RES_ADJ_PERT_GEN_TIME] = (double)ptr;  

      ptr = NewStat("ADJ_PERT_LIFETIME", 1, 3);
      WDB[RES_ADJ_PERT_LIFETIME] = (double)ptr;  

      ptr = NewStat("ADJ_PERT_BETA_EFF", 1, 9);
      WDB[RES_ADJ_PERT_BETA_EFF] = (double)ptr;  

      ptr = NewStat("ADJ_PERT_ROSSI_ALPHA", 1, 1);
      WDB[RES_ADJ_PERT_ROSSI_ALPHA] = (double)ptr;
    }

  /* Misc. analog time constants */

  ptr = NewStat("ANA_DELAYED_EMTIME", 1, 1);
  WDB[RES_ANA_DELAYED_EMTIME] = (double)ptr;  

  ptr = NewStat("ANA_SLOW_TIME", 1, 3);
  WDB[RES_ANA_SLOW_TIME] = (double)ptr;  

  ptr = NewStat("ANA_THERM_TIME", 1, 3);
  WDB[RES_ANA_THERM_TIME] = (double)ptr;  

  ptr = NewStat("ANA_THERM_FRAC", 1, 3);
  WDB[RES_ANA_THERM_FRAC] = (double)ptr;  

  ptr = NewStat("ANA_PHOTON_LIFETIME", 1, 1);
  WDB[RES_ANA_PHOTON_LIFETIME] = (double)ptr;  

  /***************************************************************************/

  /***** Cycle statistics ****************************************************/

  /* Mean population size */

  ptr = NewStat("MEAN_POP_SIZE", 1, 1);
  AllocStatHistory(ptr);
  WDB[RES_MEAN_POP_SIZE] = (double)ptr;  

  /* Mean population weight */

  ptr = NewStat("MEAN_POP_WGT", 1, 1);
  AllocStatHistory(ptr);
  WDB[RES_MEAN_POP_WGT] = (double)ptr;  

  /* Cycle-wise running time */

  ptr = NewStat("TRANSPORT_RUNTIME", 1, 2);
  AllocStatHistory(ptr);
  WDB[RES_CYCLE_RUNTIME] = (double)ptr;  

  /* CPU usage in transport cycle */

  ptr = NewStat("TRANSPORT_CPU_USAGE", 1, 1);
  AllocStatHistory(ptr);
  WDB[RES_CPU_USAGE] = (double)ptr;  

  /***************************************************************************/

  /***** Photon integral parameters ******************************************/

  ptr = NewStat("TOT_PHOTON_LEAKRATE", 1, 1);
  WDB[RES_TOT_PHOTON_LEAKRATE] = (double)ptr;  

  ptr = NewStat("TOT_PHOTON_LOSSRATE", 1, 1);
  WDB[RES_TOT_PHOTON_LOSSRATE] = (double)ptr;  

  ptr = NewStat("TOT_PHOTON_SRCRATE", 1, 1);
  WDB[RES_TOT_PHOTON_SRCRATE] = (double)ptr;  

  ptr = NewStat("TOT_PHOTON_CUTRATE", 1, 1);
  WDB[RES_TOT_PHOTON_CUTRATE] = (double)ptr;  

  ptr = NewStat("TOT_PHOTON_RR", 1, 1);
  WDB[RES_TOT_PHOTON_RR] = (double)ptr;  

  ptr = NewStat("TOT_PHOTON_FLUX", 1, 1);
  WDB[RES_TOT_PHOTON_FLUX] = (double)ptr;  

  ptr = NewStat("TOT_PHOTON_HEATRATE", 1, 1);
  WDB[RES_TOT_PHOTON_HEATRATE] = (double)ptr;  

  ptr = NewStat("PAIRPROD_CAPT_RATE", 1, 1);
  WDB[RES_PAIRPROD_CAPT_RATE] = (double)ptr;  

  ptr = NewStat("PHOTOELE_CAPT_RATE", 1, 1);
  WDB[RES_PHOTOELE_CAPT_RATE] = (double)ptr;  

  /***************************************************************************/

  /***** Dynamic simulation **************************************************/

  /* Check time cut-off */

  if (RDB[DATA_DYN_TMAX] < INFTY)
    {
      /* Number of time bins */

      np = (long)RDB[DATA_DYN_NB];
      
      ptr = NewStat("DYN_POP", 1, np);
      WDB[RES_DYN_POP] = (double)ptr;
      
      ptr = NewStat("DYN_PERIOD", 1, np);
      WDB[RES_DYN_PERIOD] = (double)ptr;
    }

  /***************************************************************************/

  /***** Collsion and reaction sampling **************************************/

  ptr = NewStat("ST_FRAC", 1, 2);
  WDB[RES_ST_TRACK_FRAC] = (double)ptr;

  ptr = NewStat("DT_FRAC", 1, 2);
  WDB[RES_DT_TRACK_FRAC] = (double)ptr;

  ptr = NewStat("DT_EFF", 1, 4);
  WDB[RES_DT_TRACK_EFF] = (double)ptr;

  ptr = NewStat("IFC_COL_EFF", 1, 2);
  WDB[RES_IFC_COL_EFF] = (double)ptr;

  ptr = NewStat("TOT_COL_EFF", 1, 4);
  WDB[RES_TOT_COL_EFF] = (double)ptr;

  ptr = NewStat("REA_SAMPLING_EFF", 1, 4);
  WDB[RES_REA_SAMPLING_EFF] = (double)ptr;

  ptr = NewStat("REA_SAMPLING_FAIL", 1, 2);
  WDB[RES_REA_SAMPLING_FAIL] = (double)ptr;

  ptr = NewStat("TMS_SAMPLING_EFF", 1, 2);
  WDB[RES_TMS_SAMPLING_EFF] = (double)ptr;

  ptr = NewStat("AVG_TRACKS", 1, 2);
  WDB[RES_AVG_TRACKS] = (double)ptr;

  ptr = NewStat("AVG_SURF_CROSS", 1, 2);
  WDB[RES_AVG_SURF_CROSS] = (double)ptr;

  ptr = NewStat("AVG_REAL_COL", 1, 2);
  WDB[RES_AVG_REAL_COL] = (double)ptr;

  ptr = NewStat("AVG_VIRT_COL", 1, 2);
  WDB[RES_AVG_VIRT_COL] = (double)ptr;

  ptr = NewStat("MIN_MACROXS", 1, 2);
  WDB[RES_MIN_MACROXS] = (double)ptr;

  ptr = NewStat("SRC_SAMPLING_EFF", 1, 4);
  WDB[RES_SRC_SAMPLING_EFF] = (double)ptr;

  ptr = NewStat("SRC_MEAN_WGT", 1, 2);
  WDB[RES_SRC_MEAN_WGT] = (double)ptr;

  ptr = NewStat("AVG_TRACK_LOOPS", 1, 4);
  WDB[RES_AVG_TRACK_LOOPS] = (double)ptr;

  /* Store last global stat */

  WDB[DATA_LAST_GLOBAL_STAT] = (double)ptr;
  
  /* TMS failure statistics (0 = total samples, 1 = majorant fail, */
  /* 2 = lower limit fail, 3 = upper limit fail) */

  ptr = NewStat("TMS_FAIL_STAT", 1, 4);
  WDB[RES_TMS_FAIL_STAT] = (double)ptr;

  /***************************************************************************/

  /***** Equilibrium Xe-135 and Sm-149 calculation ***************************/
  
  /* Check equilibrium Xe-135 calculation */

  if ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] > -1)
    {
      ptr = NewStat("I135_EQUIL_CONC", 1, 1);
      WDB[RES_I135_EQUIL_CONC] = (double)ptr;

      ptr = NewStat("XE135_EQUIL_CONC", 1, 1);
      WDB[RES_XE135_EQUIL_CONC] = (double)ptr;

      ptr = NewStat("XE135_ABS_RATE", 1, 1);
      WDB[RES_XE135_ABSRATE] = (double)ptr;

      /* Store last global stat */

      WDB[DATA_LAST_GLOBAL_STAT] = (double)ptr;
    }

  /* Check equilibrium Sm-149 calculation */

  if ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] > -1)
    {
      ptr = NewStat("PM149_EQUIL_CONC", 1, 1);
      WDB[RES_PM149_EQUIL_CONC] = (double)ptr;

      ptr = NewStat("SM149_EQUIL_CONC", 1, 1);
      WDB[RES_SM149_EQUIL_CONC] = (double)ptr;

      ptr = NewStat("SM149_ABS_RATE", 1, 1);
      WDB[RES_SM149_ABSRATE] = (double)ptr;

      /* Store last global stat */

      WDB[DATA_LAST_GLOBAL_STAT] = (double)ptr;
    }

  /***************************************************************************/
  
  /***** Normalization *******************************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_NORM];
  while (loc0 > VALID_PTR)
    {
      /* Allocate memory */

      ptr = NewStat("NORM_FISSRATE", 1, 1);
      WDB[loc0 + NORM_PTR_FISSRATE] = (double)ptr;

      ptr = NewStat("NORM_NSF", 1, 1);
      WDB[loc0 + NORM_PTR_NSF] = (double)ptr;

      ptr = NewStat("NORM_FISSE", 1, 1);
      WDB[loc0 + NORM_PTR_FISSE] = (double)ptr;

      ptr = NewStat("NORM_NEUTRON_FLUX", 1, 1);
      WDB[loc0 + NORM_PTR_NEUTRON_FLUX] = (double)ptr;

      ptr = NewStat("NORM_PHOTON_FLUX", 1, 1);
      WDB[loc0 + NORM_PTR_PHOTON_FLUX] = (double)ptr;

      ptr = NewStat("NORM_PHOTON_HEATRATE", 1, 1);
      WDB[loc0 + NORM_PTR_PHOTON_HEATRATE] = (double)ptr;

      /* Store last global stat */

      WDB[DATA_LAST_GLOBAL_STAT] = (double)ptr;

      /* Next */
      
      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Surface currents for weight window importance ***********************/

  /* Loop over weight windows */

  loc0 = (long)RDB[DATA_PTR_WWD0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over mesh */
      
      loc1 = (long)RDB[loc0 + WWD_PTR_MESH_DATA];
      while (loc1 > VALID_PTR)
	{
	  /* Allocate memory for stats */
	  
	  ptr = NewStat("WWD_MESH_IN_CURR", 1, 1);
	  WDB[loc1 + WWD_MESH_RES_IN_CURR] = (double)ptr;  

	  ptr = NewStat("WWD_MESH_OUT_CURR", 1, 1);
	  WDB[loc1 + WWD_MESH_RES_OUT_CURR] = (double)ptr;  

	  ptr = NewStat("WWD_MESH_SRC_RATE", 1, 1);
	  WDB[loc1 + WWD_MESH_RES_SRC_RATE] = (double)ptr;  

	  /* Next */

	  loc1 = NextItem(loc1);
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
