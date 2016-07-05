/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setnormalization.c                             */
/*                                                                           */
/* Created:       2013/06/24 (JLe)                                           */
/* Last modified: 2013/06/24 (JLe)                                           */
/* Version:       2.1.14                                                     */
/*                                                                           */
/* Description: Sets normalization for depletion interval                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetNormalization:"

/*****************************************************************************/

void SetNormalization(long dep)
{
  long mat0, mat, ptr, norm;
  
  /* Set normalization (NOTE: toi norm-rakenne ei tämän jälkeen  */
  /* välttämättä toimi enää linkitettynä listana, mutta sillä ei */
  /* pitäisi olla väliä koska sen yli ei loopata enää missään.   */
  
  WDB[DATA_PTR_NORM] = RDB[dep + DEP_HIS_PTR_NORM];

  /* Get pointer to normalization */

  if ((norm = (long)RDB[dep + DEP_HIS_PTR_NORM]) < VALID_PTR)
    return;

  /* Get material pointer */

  if ((mat0 = (long)RDB[norm + NORM_PTR_MAT]) < VALID_PTR)
    return;

  /* Reset pointers in materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Reset pointer */

      WDB[mat + MATERIAL_PTR_NORM] = NULLPTR;

      /* Next material */

      mat = NextItem(mat);
    }

  /* Set pointers in materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer */

      if (mat == mat0)
	WDB[mat + MATERIAL_PTR_NORM] = (double)norm;

      /* Check parent */

      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
	if (ptr == mat0)
	  WDB[mat + MATERIAL_PTR_NORM] = (double)norm;

      /* Next material */

      mat = NextItem(mat);
    }
}

/*****************************************************************************/
