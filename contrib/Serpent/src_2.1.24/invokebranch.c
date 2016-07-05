/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : invokebranch.c                                 */
/*                                                                           */
/* Created:       2014/04/08 (JLe)                                           */
/* Last modified: 2015/03/10 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Invokes depletion branch in automated burnup sequence for    */
/*              group constant generation                                    */
/*                                                                           */
/* Comments: - Called right after ReadInput() to invoke changes before       */
/*             anything else is done.                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InvokeBranch:"

/*****************************************************************************/

void InvokeBranch(long loc0)
{
  long loc1, loc2, mat, cell, lat, nst, ptr;
  double sum;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /***************************************************************************/

  /***** Replace materials ***************************************************/

  /* Loop over branches */

  loc1 = (long)RDB[loc0 + DEP_BRA_PTR_REPLACE_MAT];
  while (loc1 > VALID_PTR)
    {
      /* Find first material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
	{
	  /* Compare name */
	  
	  if (CompareStr(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT1, 
			 mat + MATERIAL_PTR_NAME))
	    break;

	  /* Next material */

	  mat = NextItem(mat);
	}

      /* Check pointer */

      if (mat < VALID_PTR)
	Die(FUNCTION_NAME, "Material %s is not defined", 
	      GetText(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT1));

      /* Remove original */

      RemoveItem(mat);

      /* Find second material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
	{
	  /* Compare name */
	  
	  if (CompareStr(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT2, 
			 mat + MATERIAL_PTR_NAME))
	    break;

	  /* Next material */

	  mat = NextItem(mat);
	}

      /* Check pointer */

      if (mat < VALID_PTR)
	Die(FUNCTION_NAME, "Material %s is not defined", 
	    GetText(loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT2));

      /* Override name */

      WDB[mat + MATERIAL_PTR_NAME] = RDB[loc1 + DEP_BRA_REPLACE_MAT_PTR_MAT1];

      /* Reset replaced flag */

      ResetOption(mat + MATERIAL_OPTIONS, OPT_REPLACED_MAT);

      /* Next */

      loc1 = NextItem(loc1);
    }

  /***************************************************************************/
  
  /***** Replace universes ***************************************************/

  /* Loop over branches */

  loc1 = (long)RDB[loc0 + DEP_BRA_PTR_REPLACE_UNI];
  while (loc1 > VALID_PTR)
    {
      /* Loop over cells */

      cell = (long)RDB[DATA_PTR_C0];
      while (cell > VALID_PTR)
	{
	  /* Compare with first and second universe */

	  if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1, 
			 cell + CELL_PTR_UNI))	  
	    {
	      /* Copy pointer */

	      ptr = cell;

	      /* Pointer to next */
	  
	      cell = NextItem(cell);

	      /* Remove cell from list */

	      RemoveItem(ptr);

	      /* Cycle loop */

	      continue;
	    }
	  else if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2, 
			      cell + CELL_PTR_UNI))	  
	    {
	      /* Replace pointer */

	      WDB[cell + CELL_PTR_UNI] = 
		RDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1];
	    }

	  /* Next cell */
	  
	  cell = NextItem(cell);
	}

      /* Loop over lattices */

      lat = (long)RDB[DATA_PTR_L0];
      while (lat > VALID_PTR)
	{
	  /* Compare with first and second universe */

	  if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1, 
			 lat + LAT_PTR_NAME))	  
	    {
	      /* Copy pointer */

	      ptr = lat;

	      /* Pointer to next */
	  
	      lat = NextItem(lat);

	      /* remove lattice from list */

	      RemoveItem(ptr);

	      /* Cycle loop */

	      continue;
	    }
	  else if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2, 
			      lat + LAT_PTR_NAME))	  
	    {
	      /* Replace pointer */

	      WDB[lat + LAT_PTR_NAME] = 
		RDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1];
	    }

	  /* next lattice */
	  
	  lat = NextItem(lat);
	}

      /* Loop over nests */

      nst = (long)RDB[DATA_PTR_NST0];
      while (nst > VALID_PTR)
	{
	  /* Compare with first and second universe */

	  if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1, 
			 nst + NEST_PTR_NAME))	  
	    {
	      /* Copy pointer */

	      ptr = nst;

	      /* Pointer to next */
	  
	      nst = NextItem(nst);

	      /* Remove nest from list */

	      RemoveItem(ptr);

	      /* Cycle loop */

	      continue;
	    }
	  else if (CompareStr(loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI2, 
			      nst + NEST_PTR_NAME))	  
	    {
	      /* Replace pointer */

	      WDB[nst + NEST_PTR_NAME] = 
		RDB[loc1 + DEP_BRA_REPLACE_UNI_PTR_UNI1];
	    }

	  /* next nest */
	  
	  nst = NextItem(nst);
	}

      /* Next */

      loc1 = NextItem(loc1);
    }

  /***************************************************************************/

  /***** Changes in material states ******************************************/

  /* Loop over branches */

  loc1 = (long)RDB[loc0 + DEP_BRA_PTR_STP];
  while (loc1 > VALID_PTR)
    {
      /* Find material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
	{
	  /* Compare name */

	  if (CompareStr(loc1 + DEP_BRA_STP_PTR_MAT, mat + MATERIAL_PTR_NAME))
	    break;

	  /* Next material */

	  mat = NextItem(mat);
	}

      /* Check pointer */

      if (mat < VALID_PTR)
	Error(loc0, "Material %s is not defined", 
	      GetText(loc1 + DEP_BRA_STP_PTR_MAT));

      /* Override density and temperature */

      WDB[mat + MATERIAL_ADENS] = RDB[loc1 + DEP_BRA_STP_DENSITY];

      if (RDB[loc1 + DEP_BRA_STP_TEMP] > 0.0)
	WDB[mat + MATERIAL_COEF_TEMP] = RDB[loc1 + DEP_BRA_STP_TEMP];

      /* Check if sum is given */

      if (RDB[mat + MATERIAL_ADENS] == -INFTY)
	{
	  /* Loop over composition and calculate sum */

	  sum = 0.0;
	  ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
	  while (ptr > VALID_PTR)
	    {
	      /* Add to sum */

	      sum = sum + RDB[ptr + COMPOSITION_ADENS];

	      /* Next */
	      
	      ptr = NextItem(ptr);
	    }

	  /* Put value */

	  WDB[mat + MATERIAL_ADENS] = sum;
	}

      /* Override S(a,b) data */

      if ((loc2 = (long)RDB[loc1 + DEP_BRA_STP_PTR_SAB]) > VALID_PTR)
	{
	  /* Reset pointer */

	  WDB[mat + MATERIAL_PTR_SAB] = NULLPTR;

	  /* Loop over data */
	  
	  while (loc2 > VALID_PTR)
	    {
	      /* Create new item */

	      ptr = NewItem(mat + MATERIAL_PTR_SAB, THERM_BLOCK_SIZE);

	      /* Put name and ZA */
		  
	      WDB[ptr + THERM_PTR_ALIAS] = RDB[loc2 + DEP_BRA_STP_SAB_PTR_LIB];
	      WDB[ptr + THERM_ZA] = RDB[loc2 + DEP_BRA_STP_SAB_ZA];

	      /* Pointer to next */

	      loc2 = NextItem(loc2);
	    }
	}

      /* Next */

      loc1 = NextItem(loc1);
    }

  /***************************************************************************/

  /***** Other stuff *********************************************************/

  /* Replace universes for group constant generation */
  
  if ((loc1 = (long)RDB[loc0 + DEP_BRA_PTR_GCU]) > VALID_PTR)
    {
      /* Get pointer to existing, and check that list is single-valued */

      if ((ptr = (long)RDB[DATA_PTR_GCU0]) > VALID_PTR)
	if (ListSize(ptr) > 1)
	  Error(loc0, "Replaced GCU-list must be single-valued");

      /* Replace pointer */
      
      WDB[DATA_PTR_GCU0] = (double)loc1;

      /* Loop over ADF's */

      ptr = (long)RDB[DATA_PTR_ADF0];
      while (ptr > VALID_PTR)
	{
	  /* Replace universe name */

	  WDB[ptr + ADF_PTR_GCU] = RDB[loc1 + GCU_PTR_UNIV];

	  /* Pointer to next */

	  ptr = NextItem(ptr);
	}

      /* Loop over pin power distributions */

      ptr = (long)RDB[DATA_PTR_PPW0];
      while (ptr > VALID_PTR)
	{
	  /* Replace universe name */

	  WDB[ptr + PPW_PTR_GCU] = RDB[loc1 + GCU_PTR_UNIV];

	  /* Pointer to next */

	  ptr = NextItem(ptr);
	}

      /* Loop over ALB's */

      ptr = (long)RDB[DATA_PTR_ALB0];
      while (ptr > VALID_PTR)
	{
	  /* Replace universe name */

	  WDB[ptr + ALB_PTR_GCU] = RDB[loc1 + GCU_PTR_UNIV];

	  /* Pointer to next */

	  ptr = NextItem(ptr);
	}	
    }

  /***************************************************************************/
}

/*****************************************************************************/
