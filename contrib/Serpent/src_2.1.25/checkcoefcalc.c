/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkcoefcalc.c                                */
/*                                                                           */
/* Created:       2014/08/18 (JLe)                                           */
/* Last modified: 2016/02/20 (JLe)                                           */
/* Version:       2.1.25                                                     */
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
  long loc0, loc1, loc2, ns, n, ptr, bra, mat, tra1, tra2, sab1, sab2;

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

	      /***** Check transformations ***********************************/

	      ptr = (long)RDB[bra + DEP_BRA_PTR_TRANS];
	      while (ptr > VALID_PTR)
		{
		  /* Reset pointer */

		  tra1 = -1;

		  /* Find transformation */
		  
		  tra2 = (long)RDB[DATA_PTR_TR0];
		  while (tra2 > VALID_PTR)
		    {
		      /* Compare names */
		      
		      if (CompareStr(ptr + DEP_BRA_TRANS_PTR_TRANS, 
				     tra2 + TRANS_PTR_NAME))
			{
			  /* Set pointer if not defined */
			  
			  if (tra1 < VALID_PTR)
			    tra1 = tra2;
			  else
			    Error(bra, "Multiple transformations named %s",
				  GetText(ptr + DEP_BRA_TRANS_PTR_TRANS));
			}
		      
		      /* Next */
		      
		      tra2 = NextItem(tra2);
		    }

		  /* Check pointer */

		  if (tra1 < VALID_PTR)
		    Error(bra, "Transformation %s is not defined", 
			  GetText(ptr + DEP_BRA_TRANS_PTR_TRANS));
 
		  /* Next */

		  ptr = NextItem(ptr);
		}
 
	      /***************************************************************/

	      /***** Checks for state-point variation ************************/

	      ptr = (long)RDB[bra + DEP_BRA_PTR_STP];
	      while (ptr > VALID_PTR)
		{
		  /* Find material */

		  mat = (long)RDB[DATA_PTR_M0];
		  while (mat > VALID_PTR)
		    {
		      /* Compare name */
		      
		      if (CompareStr(ptr + DEP_BRA_STP_PTR_MAT, 
				     mat + MATERIAL_PTR_NAME))
			break;
		      
		      /* Next material */
		      
		      mat = NextItem(mat);
		    }
		  
		  /* Check pointer */
		  
		  if (mat < VALID_PTR)
		    Error(bra, "Material %s is not defined", 
			  GetText(ptr + DEP_BRA_STP_PTR_MAT));

		  /* Check if material has S(a,b) data */

		  sab2 = (long)RDB[mat + MATERIAL_PTR_SAB];
		  while (sab2 > VALID_PTR)
		    {
		      /* Find match */

		      sab1 = (long)RDB[ptr + DEP_BRA_STP_PTR_SAB];
		      while (sab1 > VALID_PTR)
			{
			  /* Compare */

			  if (CompareStr(sab1 + DEP_BRA_STP_SAB_PTR_THERM,
					 sab2 + THERM_PTR_ALIAS))
			    break;
			  
			  /* Next */
			  
			  sab1 = NextItem(sab1);
			}
		      
		      /* Check */
		      
		      if (sab1 < VALID_PTR)
			Error(bra, "Missing entry for %s",
			      GetText(sab2 + THERM_PTR_ALIAS));

		      /* Next */
		      
		      sab2 = NextItem(sab2);
		    }

		  /* Loop over data */

		  sab1 = (long)RDB[ptr + DEP_BRA_STP_PTR_SAB];
		  while (sab1 > VALID_PTR)
		    {
		      
		      /* Check that names match */

		      sab2 = (long)RDB[mat + MATERIAL_PTR_SAB];
		      while (sab2 > VALID_PTR)
			{
			  /* Compare */
			  
			  if (CompareStr(sab1 + DEP_BRA_STP_SAB_PTR_THERM,
					 sab2 + THERM_PTR_ALIAS))
			    break;
			  
			  /* Next */
			  
			  sab2 = NextItem(sab2);
			}

		      /* Check pointer */
		      
		      if (sab2 < VALID_PTR)
			Error(bra, 
			      "Material %s has no thermal scattering data %s",
			      GetText(mat + MATERIAL_PTR_NAME),
			      GetText(sab1 + DEP_BRA_STP_SAB_PTR_THERM));
		      
		      
		      /* Next */
		      
		      sab1 = NextItem(sab1);		  
		    }

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
