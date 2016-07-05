/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatemasses.c                              */
/*                                                                           */
/* Created:       2011/05/31 (JLe)                                           */
/* Last modified: 2015/05/02 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Calculates material and fissile masses etc.                  */
/*                                                                           */
/* Comments: - NOTE: tän kanssa pitää olla varovainen, sillä palamat ym.     */
/*             tehdään alkuperäiseen fissiiliin massaan, eli niitä ei        */
/*             päivitetä ensimmäisen askeleen jälkeen. Serpent 1:ssä on      */
/*             optio sille että tuo massa muuttuu, mutta se on nyt jätetty   */
/*             tästä pois.                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateMasses:"

/*****************************************************************************/

void CalculateMasses()
{
  long mat, iso, nuc, ptr;
  double mass, mdens, adens;

  /* Reset total masses */

  WDB[DATA_TOT_FMASS] = 0.0;
  WDB[DATA_TOT_BURN_FMASS] = 0.0;

  /* Reset initial masses */

  if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
    {
      /* Reset totals */

      WDB[DATA_INI_FMASS] = 0.0;
      WDB[DATA_INI_BURN_FMASS] = 0.0;

      /* Loop over materials */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
	{
	  /* Reset material-wise masses */

	  WDB[mat + MATERIAL_INI_FMASS] = 0.0;

	  /* Next */

	  mat = NextItem(mat);
	}
    }

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get material density */
      
      mdens = RDB[mat + MATERIAL_MDENS];

      /* Check value */

      if (mdens < ZERO)
	{
	  /* Check if used in continuous reprocessing */

	  if ((long)RDB[mat + MATERIAL_FLOW_PTR_FIRST] > VALID_PTR)
	    {
	      /* Pointer to next */

	      mat = NextItem(mat);

	      /* Cycle loop */

	      continue;
	    }
	  else      
	    Die(FUNCTION_NAME, "material %s mass density too low (%E)",
		GetText(mat + MATERIAL_PTR_NAME), mdens);
	}

      /* Calculate material mass */
      
      WDB[mat + MATERIAL_MASS] = RDB[mat + MATERIAL_VOLUME]*mdens;

      /* Reset mass */

      mass = 0.0;

      /* Loop over composition and calculate mass of fissile isotopes */

      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
	{
	  /* Get atomic density and pointer to nuclide */

	  adens = RDB[iso + COMPOSITION_ADENS];
	  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	  /* Check fissile flag and add to mass */

	  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
	    mass = mass + adens*RDB[nuc + NUCLIDE_AW]*
	      RDB[mat + MATERIAL_VOLUME]/N_AVOGADRO;

	  /* Next isotope */

	  iso = NextItem(iso);
	}	

      /* Check fissile flag and burnup step */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) &&
	  !((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
	{
	  /* Exlcude divided materials */

	  if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
	    {
	      /* Add to total fissile mass */
	      
	      WDB[DATA_TOT_FMASS] = RDB[DATA_TOT_FMASS] + 1E-3*mass;
	      
	      if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
		WDB[DATA_INI_FMASS] = RDB[DATA_INI_FMASS] + 1E-3*mass;

	      /* Put material-wise mass */

	      if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
		{
		  /* Put material-wise mass */
		  
		  WDB[mat + MATERIAL_INI_FMASS] = 1E-3*mass;

		  /* Add to parent */

		  if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) 
		      > VALID_PTR)
		    WDB[ptr + MATERIAL_INI_FMASS] = 
		      RDB[ptr + MATERIAL_INI_FMASS] + 1E-3*mass;
		}

	      /* Check burn-flag */
	      
	      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
		{
		  /* Add to burnable fissile mass */
		  
		  WDB[DATA_TOT_BURN_FMASS] = RDB[DATA_TOT_BURN_FMASS] +
		    1E-3*mass;
		  
		  if ((long)RDB[DATA_BURN_CALC_INI_MASS] == YES)
		    WDB[DATA_INI_BURN_FMASS] = RDB[DATA_INI_BURN_FMASS] + 
		      1E-3*mass;
		}
	    }
	}
      
      /* Check burn flag and mass */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
	  (RDB[mat + MATERIAL_MASS] <= 0.0) &&
	  ((long)RDB[mat + MATERIAL_VOL_COUNT] > 0))
	Error(mat, "Volume or mass of material \"%s\" must be given", 
	      GetText(mat + MATERIAL_PTR_NAME));
      
      /* Next material */

      mat = NextItem(mat);
    }

  /* Reset initial mass flag */

  WDB[DATA_BURN_CALC_INI_MASS] = (double)NO;
}

/*****************************************************************************/
