/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setcoefcalc.c                                  */
/*                                                                           */
/* Created:       2014/04/15 (JLe)                                           */
/* Last modified: 2014/08/30 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Sets depletion branches for coefficient calculations         */
/*                                                                           */
/* Comments: - Called right after ReadInput() to invoke changes before       */
/*             anything else is done.                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetCoefCalc:"

/*****************************************************************************/

long SetCoefCalc(long idx)
{
  long loc0, loc1, loc2, ptr, n, i, m, nc;
  double x;

  /* Check pointer to data */

  if ((loc0 = (long)RDB[DATA_PTR_COEF0]) > VALID_PTR)
    {
      /* Set flag */

      WDB[DATA_MORE_COEF_CALC] = (double)YES;

      /* Check burnup points */

      while (loc0 > VALID_PTR)
	{
	  /* Pointer to burnup points */
	  
	  loc1 = (long)RDB[loc0 + COEF_PTR_BU_PTS];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  
	  /* Loop over points and check values */
	  
	  for (n = 0; n < (long)RDB[loc0 + COEF_N_BU]; n++)
	    {
	      /* Get value */
	      
	      x = atof(GetText(loc1++));
	      
	      /* Check value */
	      
	      if ((x != 0.0) && ((long)RDB[DATA_BURN_PTR_DEP] < VALID_PTR))
		Error(loc0, 
		      "Non-zero restart points without burnup calculation");
	    }
	
	  /* Next */

	  loc0 = NextItem(loc0);
	}

      /* Start from first branch if no burnup calculation */

      if ((((long)RDB[DATA_BURN_PTR_DEP]) < VALID_PTR) && (idx < 0))
	idx = 0;
    }

  /* Check index */

  if (idx < 0)
    return idx;
    
  fprintf(out, "Re-configuring model for coefficient calculation...\n");

  /* Put index */

  WDB[DATA_COEF_CALC_IDX] = (double)(idx + 1);

  /* Loop over coefficient matrixes */

  loc0 = (long)RDB[DATA_PTR_COEF0];
  while (loc0 > VALID_PTR)
    {
      /* Check index */

      if (idx > (long)RDB[loc0 + COEF_N_TOT] - 1)
	{
	  /* Update */

	  idx = idx - (long)RDB[loc0 + COEF_N_TOT];

	  /* Add to current run index */

	  WDB[DATA_COEF_CALC_RUN_IDX] = RDB[DATA_COEF_CALC_RUN_IDX] +
	    (long)(RDB[loc0 + COEF_N_TOT]*RDB[loc0 + COEF_N_BU]);

	  /* Next matrix */

	  loc0 = NextItem(loc0);

	  /* Cycle loop */

	  continue;
	}

      /* Loop over branch matrix */

      for (n = 0; n < (long)RDB[loc0 + COEF_N_TOT]; n++)
	{
	  /* Reset cumulative */

	  nc = 0;

	  /* Pointer to matrix */

	  loc1 = (long)RDB[loc0 + COEF_PTR_MTX];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  loc1 = LastItem(loc1);

	  /* Loop over matrix */

	  while (loc1 > VALID_PTR)
	    {
	      /* Calculate index */

	      m = (long)(((double)n)/RDB[loc1 + COEF_MTX_N_CUMU]);

	      /* Check index and invoke branch */

	      if (n == idx)
		{
		  /* Pointer to branches */

		  ptr = (long)RDB[loc1 + COEF_MTX_PTR_BRA];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  
		  /* Loop to value */
		  
		  i = m - (long)(((double)nc)/RDB[loc1 + COEF_MTX_N_CUMU]); 
		  while (i-- > 0)
		    ptr++;

		  /* Pointer to branches */

		  if ((loc2 = (long)RDB[DATA_PTR_BRA0]) < VALID_PTR)
		    Error(loc0, "No branches defined");

		  /* Loop over branches and find match */

		  while (loc2 > VALID_PTR)
		    {
		      /* Check */
		      
		      if (CompareStr(loc2 + DEP_BRA_PTR_NAME, ptr))
			break;

		      /* Next */
		      
		      loc2 = NextItem(loc2);
		    }

		  /* Check pointer, nominal and invoke branch */
		  
		  if (loc2 < VALID_PTR)
		    Die(FUNCTION_NAME, "Branch \"%s\" not defined", 
			GetText(ptr));
		  else if (loc2 > VALID_PTR)
		    InvokeBranch(loc2);

		  /* Put pointer to branch name for output */

		  WDB[loc1 + COEF_MTX_PTR_BRA] = (double)ptr;

		  /* Put pointer to variable list */

		  if (loc2 > VALID_PTR)
		    WDB[loc1 + COEF_MTX_PTR_VAR] = RDB[loc2 + DEP_BRA_PTR_VAR];
		}
	    
	      /* Update cumulative */

	      nc = ((long)RDB[loc1 + COEF_MTX_N_CUMU])*m;

	      /* Pointer to next */

	      loc1 = PrevItem(loc1);
	    }

	  /* Check index exit subroutine */

	  if (n == idx)
	    {
	      /* Add to current run index */

	      WDB[DATA_COEF_CALC_RUN_IDX] = RDB[DATA_COEF_CALC_RUN_IDX] +
		(double)(nc*RDB[loc0 + COEF_N_BU]);

	      /* Check if last */

	      if ((NextItem(loc0) < VALID_PTR) &&
		  (idx == (long)RDB[loc0 + COEF_N_TOT] - 1))
		WDB[DATA_MORE_COEF_CALC] = (double)NO;
	      else
		WDB[DATA_MORE_COEF_CALC] = (double)YES;

	      /* Set pointer */

	      WDB[DATA_PTR_COEF0] = (double)loc0;

	      /* Exit OK */

	      fprintf(out, "OK.\n\n");

	      /* Exit subroutine */

	      return idx;
	    }
	}

      /* Next coefficient matrix */
      
      loc0 = NextItem(loc0);
    }

  /* Something wrong here */

  Die(FUNCTION_NAME, "Shouldn't be here");

  /* Avoid compiler warning */

  return -1;
}

/*****************************************************************************/
