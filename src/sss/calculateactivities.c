#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculateactivities.c                          */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2015/10/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Calculates activities, decay heat and spontaneous fission    */
/*              rates for materials                                          */
/*                                                                           */
/* Comments: - Fotonitransportmoodissa missä materiaalikoostumukset on       */
/*             annettu neutroni / decay -datana alkuperäiset listat on       */
/*             korvattu replacephotondata.c:ssä. Alkuperäinen lista luetaan  */
/*             erillisestä paikasta.                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateActivities:"

/*****************************************************************************/

void CalculateActivities()
{
  long mat, iso, nuc, ptr;
  double vol, adens, sum1, sum2, A, H, SF, HT, GT, I;

  /* Check if decay data file is given */

  if ((long)RDB[DATA_PTR_DECDATA_FNAME_LIST] < 1)
    return;

  fprintf(out, "Calculating activities...\n");
  
  /***************************************************************************/

  /***** Calculate total activities ******************************************/

  /* Reset global data */

  WDB[DATA_TOT_ACTIVITY] = 0.0;
  WDB[DATA_TOT_SFRATE] = 0.0;
  WDB[DATA_TOT_DECAY_HEAT] = 0.0;
  WDB[DATA_ACT_ACTIVITY] = 0.0;
  WDB[DATA_ACT_DECAY_HEAT] = 0.0;
  WDB[DATA_FP_ACTIVITY] = 0.0;
  WDB[DATA_FP_DECAY_HEAT] = 0.0;
  WDB[DATA_TOT_ING_TOX] = 0.0;
  WDB[DATA_TOT_INH_TOX] = 0.0;

  WDB[DATA_SR90_ACTIVITY] = 0.0;
  WDB[DATA_TE132_ACTIVITY] = 0.0;
  WDB[DATA_I131_ACTIVITY] = 0.0;
  WDB[DATA_I132_ACTIVITY] = 0.0;
  WDB[DATA_CS134_ACTIVITY] = 0.0;
  WDB[DATA_CS137_ACTIVITY] = 0.0;

  WDB[DATA_TOT_PHOTON_SRC_RATE] = 0.0;
  WDB[DATA_TOT_PHOTON_SRC_VOL] = 0.0;
  WDB[DATA_PHOTON_SRC_MAX_I] = 0.0;

  /* Reset material-wise values */
 
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset values */

      WDB[mat + MATERIAL_ACTIVITY] = 0.0;
      WDB[mat + MATERIAL_SFRATE] = 0.0;
      WDB[mat + MATERIAL_DECAY_HEAT] = 0.0;
      WDB[mat + MATERIAL_PHOTON_SRC_RATE] = 0.0;
   
      /* Pointer to next */

      mat = NextItem(mat);
    }

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check physical flag */
      
      if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))
	{
	  /* Pointer to next */

	  mat = NextItem(mat);
 
	  /* Cycle loop */

	  continue;
	}
      
      /* Check divisor type */

      if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
	Die(FUNCTION_NAME, "Divided parent material");
      
      /* Get volume */
      
      vol = RDB[mat + MATERIAL_VOLUME];
      
      /* Check (tällä on korvattu se fail-hässäkkä 20.11.2012 / 2.1.10) */
      
      if ((vol > 0.0) && (vol < 1E+18))
	{
	  /* Loop over composition */
	  
	  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) 
	      < VALID_PTR)
	    iso = (long)RDB[mat + MATERIAL_PTR_COMP];
	  
	  while (iso > VALID_PTR)
	    {
	      /* Pointer to nuclide */
	      
	      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
	      
	      /* Get atomic density */
	      
	      adens = RDB[iso + COMPOSITION_ADENS];
	      
	      /* Calculate activity, decay heat and photon source rate */
	      
	      A = vol*adens*RDB[nuc + NUCLIDE_LAMBDA]*1E+24;
	      H = A*RDB[nuc + NUCLIDE_DECAY_E]*MEV;
	      SF = A*RDB[nuc + NUCLIDE_SFBR];
	      HT = A*RDB[nuc + NUCLIDE_SPEC_INH_TOX];
	      GT = A*RDB[nuc + NUCLIDE_SPEC_ING_TOX];
	      I = A*RDB[nuc + NUCLIDE_SPEC_PHOTON_I];
	      
	      /* Check values */
	      
	      CheckValue(FUNCTION_NAME, "A", "", A, 0.0, 1E+30);
	      CheckValue(FUNCTION_NAME, "H", "", H, 0.0, 1E+30);
	      CheckValue(FUNCTION_NAME, "SF", "", SF, 0.0, A);
	      CheckValue(FUNCTION_NAME, "I", "", I, 0.0, 100.0*A);
	      
	      /* Add to totals */
	      
	      WDB[DATA_TOT_ACTIVITY] = RDB[DATA_TOT_ACTIVITY] + A;
	      WDB[DATA_TOT_DECAY_HEAT] = RDB[DATA_TOT_DECAY_HEAT] + H;
	      WDB[DATA_TOT_SFRATE] = RDB[DATA_TOT_SFRATE] + SF;
	      WDB[DATA_TOT_INH_TOX] = RDB[DATA_TOT_INH_TOX] + HT;
	      WDB[DATA_TOT_ING_TOX] = RDB[DATA_TOT_ING_TOX] + GT;
	      WDB[DATA_TOT_PHOTON_SRC_RATE] = 
		RDB[DATA_TOT_PHOTON_SRC_RATE] + I;
	      
	      /* Add to burnable material values */
	      
	      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
		{
		  WDB[DATA_BURN_DECAY_HEAT] = 
		    RDB[DATA_BURN_DECAY_HEAT] + H;
		  WDB[DATA_BURN_SFRATE] = RDB[DATA_BURN_SFRATE] + SF;
		}
	      
	      /* Add partials */
	      
	      if ((long)RDB[nuc + NUCLIDE_Z] > 89)
		{
		  WDB[DATA_ACT_ACTIVITY] = RDB[DATA_ACT_ACTIVITY] + A;
		  WDB[DATA_ACT_DECAY_HEAT] = RDB[DATA_ACT_DECAY_HEAT] + H;
		}
	      else if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & 
		       NUCLIDE_FLAG_FP)
		{
		  WDB[DATA_FP_ACTIVITY] = RDB[DATA_FP_ACTIVITY] + A;
		  WDB[DATA_FP_DECAY_HEAT] = RDB[DATA_FP_DECAY_HEAT] + H;
		}
	      
	      /* Add isotopic */
	      
	      if ((long)RDB[nuc + NUCLIDE_ZAI] == 380900)
		WDB[DATA_SR90_ACTIVITY] = RDB[DATA_SR90_ACTIVITY] + A;
	      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 521320)
		WDB[DATA_TE132_ACTIVITY] = RDB[DATA_TE132_ACTIVITY] + A;
	      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 531310)
		WDB[DATA_I131_ACTIVITY] = RDB[DATA_I131_ACTIVITY] + A;
	      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 531320)
		WDB[DATA_I132_ACTIVITY] = RDB[DATA_I132_ACTIVITY] + A;
	      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 551340)
		WDB[DATA_CS134_ACTIVITY] = RDB[DATA_CS134_ACTIVITY] + A;
	      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 551370)
		WDB[DATA_CS137_ACTIVITY] = RDB[DATA_CS137_ACTIVITY] + A;
	      
	      /* Add to material data */
	      
	      WDB[mat + MATERIAL_ACTIVITY] = 
		RDB[mat + MATERIAL_ACTIVITY] + A;
	      
	      WDB[mat + MATERIAL_DECAY_HEAT] =
		RDB[mat + MATERIAL_DECAY_HEAT] + H;
	      
	      WDB[mat + MATERIAL_SFRATE] = 
		RDB[mat + MATERIAL_SFRATE] + SF;
	      
	      WDB[mat + MATERIAL_PHOTON_SRC_RATE] = 
		RDB[mat + MATERIAL_PHOTON_SRC_RATE] + I;
	      
	      /* Next isotope */
	      
	      iso = NextItem(iso);
	    }
	  
	  /* Compare to maximum intensity */
	  
	  if (RDB[mat + MATERIAL_PHOTON_SRC_RATE]/vol > 
	      RDB[DATA_PHOTON_SRC_MAX_I])
	    WDB[DATA_PHOTON_SRC_MAX_I] = 
	      RDB[mat + MATERIAL_PHOTON_SRC_RATE]/vol;
	  
	  /* Add to photon source volume */
	  
	  if (RDB[mat + MATERIAL_PHOTON_SRC_RATE] > 0.0)
	    WDB[DATA_TOT_PHOTON_SRC_VOL] = 
	      RDB[DATA_TOT_PHOTON_SRC_VOL] + vol;
	}
      
      /* Add to parent (JLe: Tää lisättiin 1.10.2015 / 2.1.25 että  */
      /* decay gamma sourcen normeerauksen saa kiinnitettyä parent- */
      /* materiaaliin. Voi olla että tää sotkee jotain kokonais-    */
      /* arvoihin liittyviä laskuja) */
      
      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > 0.0)
	{
	  /* Sum only photon source rate for now  */
	  /* (needed in radioactive decay source) */

	  WDB[ptr + MATERIAL_PHOTON_SRC_RATE] = 
	    RDB[ptr + MATERIAL_PHOTON_SRC_RATE] + 
	    RDB[mat + MATERIAL_PHOTON_SRC_RATE];
	}

      /* Next material */

      mat = NextItem(mat);
    }

  fprintf(out, "OK.\n\n");

  /***************************************************************************/

  /***** Xenon entropy *******************************************************/

  /* Calculate mean concentration */

  sum2 = 0.0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check physical flag */
      
      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT)
	{
	  /* Check divisor type */

	  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
	    Die(FUNCTION_NAME, "Divided parent material");
	
	  /* Add Xe-135 concentration */

	  if ((iso = (long)RDB[mat + MATERIAL_PTR_I135_ISO]) > VALID_PTR)
	    sum2 = sum2 + RDB[iso + COMPOSITION_ADENS];
	}

      /* Next material */
      
      mat = NextItem(mat);
    }

  /* Calculate entropies */

  sum1 = 0.0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check physical flag */
      
      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT)
	{
	  /* Check divisor type */

	  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
	    Die(FUNCTION_NAME, "Divided parent material");

	  /* Add Xe-135 concentration */

	  if (sum2 > 0.0)
	    if ((iso = (long)RDB[mat + MATERIAL_PTR_I135_ISO]) > VALID_PTR)
	      {
		/* Get density */
	      
		adens = RDB[iso + COMPOSITION_ADENS];

		/* Add to entropy */
		
		sum1 = sum1 - log2(adens/sum2)/RDB[DATA_N_BURN_MATERIALS];
	      }
	}

      /* Next material */
      
      mat = NextItem(mat);
    }
  
  /* Put value */

  if (sum1 > 0.0)
    WDB[DATA_XENON_ENTROPY] = log2(RDB[DATA_N_BURN_MATERIALS])/sum1;
    
  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
