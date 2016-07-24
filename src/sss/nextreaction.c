#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nextreaction.c                                 */
/*                                                                           */
/* Created:       2011/11/02 (JLe)                                           */
/* Last modified: 2012/05/28 (JLe)                                           */
/* Version:       2.1.6                                                      */
/*                                                                           */
/* Description: Returns pointers and data for next reaction                  */
/*                                                                           */
/* Comments: Changed completely 28.5.2012 (2.1.6)                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NextReaction:"

/*****************************************************************************/

long NextReaction(long loc0, long *rea, double *adens, double *Emin, 
		  double *Emax, long id)
{
  long mat, loc1, next, iso, idx;

  /* Check reaction list pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Pointer to material */

  mat = (long)RDB[loc0 + RLS_PTR_MAT];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Pointer to composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

  /* Get pointer to first or next reaction in list */

  if (*rea < VALID_PTR)
    loc1 = (long)RDB[loc0 + RLS_PTR_REA0];
  else
    loc1 = TestValuePair(loc0 + RLS_PTR_NEXT, 0.0, id);
    
  /* Check pointer and put data */
      
  if (loc1 > VALID_PTR)
    {
      /* Check cut-off */

      if ((long)RDB[loc1 + RLS_DATA_CUT] == YES)
	return -1;

      /* Reaction pointer and minimum and maximum energy */

      *rea = (long)RDB[loc1 + RLS_DATA_PTR_REA];
      *Emin = RDB[loc1 + RLS_DATA_EMIN];
      *Emax = RDB[loc1 + RLS_DATA_EMAX];
      
      /* Composition index */

      idx = (long)RDB[loc1 + RLS_DATA_COMP_IDX];

      /* Pointer to composition */

      iso = ListPtr(iso, idx);
      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

      /* Atomic density */

      *adens = RDB[iso + COMPOSITION_ADENS];

      /* Pointer to next in list */
      
      next = NextItem(loc1);

      /* Update pointer */

      StoreValuePair(loc0 + RLS_PTR_NEXT, 0.0, (double)next, id);
    }

  /* Return pointer */
  
  return loc1;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
