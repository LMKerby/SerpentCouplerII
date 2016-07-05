/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoretemp.c                                    */
/*                                                                           */
/* Created:       2012/01/24 (JLe)                                           */
/* Last modified: 2012/11/04 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Scores temperatures in feedback mode                         */
/*                                                                           */
/* Comments: - Tohon tulee joku ihan ihme painotus                           */
/*           - Jos tota fluxia ei tarvii muuhun kuin tähän, niin sen voisi   */
/*             tallettaa samaan muuttujaan                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreTemp:"

/*****************************************************************************/

void ScoreTemp(long mat, double flx, long id, double wgt)
{
  long ptr, ncol, tfb, nst, reg, i, tb, nreg;
  long T;

  /* Check if feedback is in use */

  if ((long)RDB[DATA_USE_TFB] == YES)
    {
      /* Get collision number */
      
      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      ncol = (long)GetPrivateData(ptr, id);

      /* Get time bin index */
		  
      tb = (long)RDB[DATA_DYN_TB];
      
      /* Loop over temperature feedbacks */
      
      tfb = (long)RDB[DATA_PTR_TFB0];
      while (tfb > VALID_PTR)
	{
	  /* Pointer to nest */
	  
	  nst = (long)RDB[tfb + TFB_PTR_NST];
	  CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);
	  
	  /* Get pointer to region */
	  
	  if ((reg = (long)TestValuePair(nst + NEST_PTR_COL_REG, ncol, id))
	      > VALID_PTR)
	    {
	      /* Pointer to feedback region */

	      if ((reg = (long)RDB[reg + NEST_REG_PTR_TFB_REG]) > VALID_PTR)
		{
		  /* Get index and number of regions */

		  i = (long)RDB[reg + TFB_REG_IDX];
		  nreg = (long)RDB[tfb + TFB_N_REG];

		  /* Get temperature */

		  T = GetTemp(mat, id);

		  /* Score flux-averaged temperature */
		  
		  ptr = (long)RDB[tfb + TFB_PTR_MEAN_FTEMP];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  AddBuf(T*flx, wgt, ptr, id, -1, i, tb);
		  AddBuf(flx, wgt, ptr, id, -1, i + nreg, tb);

		  /* Exit subroutine */

		  return;
		}
	    }

	  /* Next feedback */

	  tfb = NextItem(tfb);
	}
    }
}

/*****************************************************************************/
