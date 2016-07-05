/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : identifytfbmat.c                               */
/*                                                                           */
/* Created:       2012/08/16 (JLe)                                           */
/* Last modified: 2012/08/21 (VVa)                                           */
/* Version:       2.1.8                                                      */
/*                                                                           */
/* Description: Identifies material typs used with temperature feedbacks     */
/*                                                                           */
/* Comments: Huomioita materiaalin tunnistukseen (16.8.2012 / 2.1.8 / JLe):  */
/*                                                                           */
/*           - Määrittelemättömän materiaalin tunniste on nyt 0, ei -1       */
/*           - Vesi tunnistetaan vertaamalla O- ja H-atomien suhdetta        */
/*             arvoon 0.5 +/- 0.1                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IdentifyTFBMat:"

/*****************************************************************************/

void IdentifyTFBMat()
{
  long tfb, reg, mat, iso, nuc;
  double adens1, adens2;

  /* Loop over feedbacks */

  tfb = (long)RDB[DATA_PTR_TFB0];
  while (tfb > VALID_PTR)
    {
      /* Loop over regions */

      reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
      while (reg > VALID_PTR)
	{
	  /* Get material pointer */

	  mat = (long)RDB[reg + TFB_REG_PTR_MAT];
	  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

	  /* Check if fissile (fuel) */

	  if (RDB[mat + MATERIAL_INI_FISS_MDENS] > 0.0)
	    {
	      /* Set type */

	      WDB[reg + TFB_REG_MAT_TYPE] = (double)TFB_MAT_TYPE_FUEL;
	  
	      /* Cycle loop */

	      reg = NextItem(reg);
	      continue;
	    }

	  /* Identify zircaloy from Zr isotopes (> 90% of composition) */

	  adens1 = 0.0;
	  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
	  while (iso > VALID_PTR) 
	    {
	      /* Pointer to nuclide */
	      
	      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	      /* Check ZA */

	      if ((RDB[nuc + NUCLIDE_ZA] > 39999) &&
		  (RDB[nuc + NUCLIDE_ZA] < 40999))
		adens1 = adens1 + RDB[iso + COMPOSITION_ADENS];

	      /* Next */

	      iso = NextItem(iso);
	    }

	  /* Check density */

	  if (adens1 > 0.9*RDB[mat + MATERIAL_ADENS])
	    {
	      /* Set type */

	      WDB[reg + TFB_REG_MAT_TYPE] = (double)TFB_MAT_TYPE_ZIRCALOY;

	      /* Cycle loop */

	      reg = NextItem(reg);
	      continue;
	    }

	  /* Identify helium gas (> 90% of composition) */

	  adens1 = 0.0;
	  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
	  while (iso > VALID_PTR) 
	    {
	      /* Pointer to nuclide */
	      
	      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	      /* Check ZA */

	      if ((RDB[nuc + NUCLIDE_ZA] > 1999) &&
		  (RDB[nuc + NUCLIDE_ZA] < 2999))
		adens1 = adens1 + RDB[iso + COMPOSITION_ADENS];

	      /* Next */

	      iso = NextItem(iso);
	    }

	  /* Check density */

	  if (adens1 > 0.9*RDB[mat + MATERIAL_ADENS])
	    {
	      /* Set type */

	      WDB[reg + TFB_REG_MAT_TYPE] = (double)TFB_MAT_TYPE_HELIUM;
	  
	      /* Cycle loop */

	      reg = NextItem(reg);
	      continue;
	    }

	  /* Identify water from ratio of hydrogen and oxygen */

	  adens1 = 0.0;
	  adens2 = 0.0;
	  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
	  while (iso > VALID_PTR) 
	    {
	      /* Pointer to nuclide */
	      
	      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	      /* Check ZA */

	      if ((RDB[nuc + NUCLIDE_ZA] > 999) &&
		  (RDB[nuc + NUCLIDE_ZA] < 1999))
		adens1 = adens1 + RDB[iso + COMPOSITION_ADENS];
	      else if ((RDB[nuc + NUCLIDE_ZA] > 7999) &&
		  (RDB[nuc + NUCLIDE_ZA] < 8999))
		adens2 = adens2 + RDB[iso + COMPOSITION_ADENS];

	      /* Next */

	      iso = NextItem(iso);
	    }
	  
	  /* Check ratio */

	  if ((adens1 > 0.0) && (fabs(adens2/adens1 - 0.5) < 0.1))
	    {
	      /* Set type */

	      WDB[reg + TFB_REG_MAT_TYPE] = (double)TFB_MAT_TYPE_WATER;
	  
	      /* Cycle loop */

	      reg = NextItem(reg);
	      continue;
	    }

	  /* Set type to undefined */

	  WDB[reg + TFB_REG_MAT_TYPE] = (double)TFB_MAT_TYPE_UNDEFINED;

	  /* Next region */

	  reg = NextItem(reg);	  
	}

      /* Next feedback */

      tfb = NextItem(tfb);
    }
}

/*****************************************************************************/
