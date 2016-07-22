/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : radgammasrc.c                                  */
/*                                                                           */
/* Created:       2012/03/29 (JLe)                                           */
/* Last modified: 2016/04/04 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Gamma source from radioactive decay                          */
/*                                                                           */
/* Comments: - Käytetään nyt noita painoja (jako kokonaislähdetilavuudella)  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RadGammaSrc:"

/*****************************************************************************/

long RadGammaSrc(long src, long mat, double *E, double *wgt, double vol, 
		 long id)
{
  long mat0, mat1, nuc, ptr;
  double I, rnd;

  /* Check source pointer */

  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Check total source rate */
      
  if (RDB[DATA_TOT_PHOTON_SRC_RATE] == 0.0)
    Error(src, "Decay source but no photon emission (check material volumes)");

  /* Check if point is in void */

  if (mat < VALID_PTR)
    return -1;

  /***************************************************************************/
  
  /***** Step 1: Rejection sampling in material ******************************/
  
  /* Check if source material is defined */
  
  if ((mat0 = (long)RDB[src + SRC_PTR_RAD_SRC_MAT]) > VALID_PTR)
    {
      /* Check pointer */
      
      if (mat != mat0)
	{
	  /* Check if material has parent */

	  if ((mat1 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
	    return -1;

	  /* Check parent */

	  if (mat1 != mat0)
	    return -1;
	}
    }
  else if (RDB[mat + MATERIAL_PHOTON_SRC_RATE] == 0.0)
    return -1;

  /* Check volume and intensity */
  
  if (RDB[mat + MATERIAL_VOLUME] == 0.0)
    Die(FUNCTION_NAME, "Volume is zero");
  else if (RDB[DATA_PHOTON_SRC_MAX_I] == 0.0)
    Die(FUNCTION_NAME, "Total intensity is zero");
  else if (RDB[DATA_TOT_PHOTON_SRC_RATE] == 0.0)
    Die(FUNCTION_NAME, "Total source rate is zero");

  /* Check mode */

  if ((long)RDB[src + SRC_RAD_SRC_MODE] == RAD_SRC_MODE_ANALOG)
    {  
      /* Perform rejection sampling on material */
      
      if (RandF(id) > RDB[mat + MATERIAL_PHOTON_SRC_RATE]/
	  RDB[mat + MATERIAL_VOLUME]/RDB[DATA_PHOTON_SRC_MAX_I])
	return -1;
    }
  else if ((long)RDB[src + SRC_RAD_SRC_MODE] == RAD_SRC_MODE_IMPLICIT)
    {
      /* Or adjust weight */
      
      *wgt = *wgt*RDB[mat + MATERIAL_PHOTON_SRC_RATE]/
	RDB[DATA_TOT_PHOTON_SRC_RATE]/RDB[mat + MATERIAL_VOLUME]*
	RDB[DATA_TOT_PHOTON_SRC_VOL];
    }
  else
    Die(FUNCTION_NAME, "Sampling mode is not set");

  /***************************************************************************/
 
  /***** Step 2: Sample nuclide **********************************************/
  
  /* This is used as cut-off in processdecaysrc.c. (could just check the */
  /* pointer, but this is better for debugging) */
  
  if (RDB[mat + MATERIAL_PHOTON_SRC_RATE]/
      RDB[DATA_TOT_PHOTON_SRC_RATE] < 1E-19)
    return -1;
  
  /* Pointer to list */

  ptr = (long)RDB[mat + MATERIAL_PTR_DECAY_SRC];
  CheckPointer(FUNCTION_NAME, "(ptr1)", DATA_ARRAY, ptr);

  /* Sample random number */

  rnd = RandF(id);

  /* Loop over list */

  while (ptr > VALID_PTR)
    {
      if (RDB[ptr + SRC_DECCAY_CUM_P] > rnd)
	break;

      /* Next */

      ptr = NextItem(ptr);
    }

 /* Check pointer */

  if (ptr < VALID_PTR)
    Die(FUNCTION_NAME, "Unable to sample nuclide %s",
	GetText(mat + MATERIAL_PTR_NAME));

  /***************************************************************************/

  /***** Step 3: Sample emission energy **************************************/

  /* Get pointer to nuclide */

  nuc = (long)RDB[ptr + SRC_DECCAY_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Get total intensity */

  I = RDB[nuc + NUCLIDE_SPEC_PHOTON_I];
  CheckValue(FUNCTION_NAME, "I", "", I, ZERO, INFTY);

  /* Sample fraction of intensity */

  I = I*RandF(id);

  /* Loop over line spectra */

  ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_LINE_SPEC];
  while (ptr > VALID_PTR)
    {
      /* Compare to intensity */

      if ((I = I - RDB[ptr + PHOTON_LINE_SPEC_RI]) < 0.0)
	break;
      
      /* Next */

      ptr = NextItem(ptr);
    }

  /* Check pointer */

  if (ptr < VALID_PTR)
    Die(FUNCTION_NAME, "Unable to sample photon line for %s",
	GetText(nuc + NUCLIDE_PTR_NAME));
  
  /* Get energy */

  *E = RDB[ptr + PHOTON_LINE_SPEC_E];

  /* Return material pointer */

  return mat;

  /***************************************************************************/
}
