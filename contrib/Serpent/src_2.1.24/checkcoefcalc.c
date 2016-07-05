/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkcoefcalc.c                                */
/*                                                                           */
/* Created:       2014/08/18 (JLe)                                           */
/* Last modified: 2014/08/28 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Checks that coefficient matrix and branches are correctly    */
/*              defined before running the calculation.                      */
/*                                                                           */
/* Comments: - The purpose is to stop the calculation early, rather than     */
/*             terminating with an error later on.                           */
/*                                                                           */
/*           - TODO: Palamapisteiden ja universumien testaamiseen pitäis     */
/*                   keksiä jotain.                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckCoefCalc:"

/*****************************************************************************/

void CheckCoefCalc()
{
  long loc0, loc1, loc2, ns, n, ptr, bra, mat;

  /* Loop over coefficient matrixes */

  loc0 = (long)RDB[DATA_PTR_COEF0];
  while (loc0 > VALID_PTR)
    {
      /* Loop over branch matrix */
      
      loc1 = (long)RDB[loc0 + COEF_PTR_MTX];
      while (loc1 > VALID_PTR)
	{
	  /* Number of branches */

	  if ((ns = (long)RDB[loc1 + COEF_MTX_N_BRA]) < 1)
	    Die(FUNCTION_NAME, "Number of branches shouldn't be zero");
	  
	  /* Loop over branches */
	  
	  loc2 = (long)RDB[loc1 + COEF_MTX_PTR_BRA];
	  for (n = 0; n < ns; n++)
	    {
	      /***************************************************************/

	      /***** Check that branches are defined *************************/

	      /* Find match */

	      bra = (long)RDB[DATA_PTR_BRA0];
	      while (bra > VALID_PTR)
		{
		  /* Compare */

		  if (CompareStr(bra + DEP_BRA_PTR_NAME, loc2))
		    break;
		  
		  /* Next */

		  bra = NextItem(bra);
		}

	      /* Check */

	      if (bra < VALID_PTR)
		Error(loc0, "Branch %s is not defined", GetText(loc2));

	      /***************************************************************/

	      /***** Check that materials are defined ************************/

	      ptr = (long)RDB[bra + DEP_BRA_PTR_REPLACE_MAT];
	      while (ptr > VALID_PTR)
		{
		  /* Find first material */
		    
		  mat = (long)RDB[DATA_PTR_M0];
		  while (mat > VALID_PTR)
		    {
		      /* Compare name */
		      
		      if (CompareStr(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT1, 
				     mat + MATERIAL_PTR_NAME))
			break;
		      
		      /* Next material */
		      
		      mat = NextItem(mat);
		    }

		  /* Check pointer */
		  
		  if (mat < VALID_PTR)
		    Error(bra, "Material %s is not defined", 
			  GetText(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT1));
		  
		  /* Find second material */
		  
		  mat = (long)RDB[DATA_PTR_M0];
		  while (mat > VALID_PTR)
		    {
		      /* Compare name */
		      
		      if (CompareStr(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT2, 
				     mat + MATERIAL_PTR_NAME))
			break;
		      
		      /* Next material */
		      
		      mat = NextItem(mat);
		    }
		  
		  /* Check pointer */
		  
		  if (mat < VALID_PTR)
		    Error(bra, "Material %s is not defined", 
			  GetText(ptr + DEP_BRA_REPLACE_MAT_PTR_MAT2));
		  
		  /* Set flag */

		  SetOption(mat + MATERIAL_OPTIONS, OPT_REPLACED_MAT);

		  /* Next */

		  ptr = NextItem(ptr);
		}

	      /***************************************************************/

	      /* Next branch in matrix */

	      loc2++;
	    }
	  
	  /* Pointer to next */

	  loc1 = NextItem(loc1);
	}
      
      /* Next coefficient matrix */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
