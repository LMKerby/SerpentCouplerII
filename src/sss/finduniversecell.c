/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : finduniversecell.c                             */
/*                                                                           */
/* Created:       2010/10/13 (JLe)                                           */
/* Last modified: 2013/10/15 (JLe)                                           */
/* Version:       2.1.16                                                     */
/*                                                                           */
/* Description: Find cell in universe                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindUniverseCell:"

/*****************************************************************************/

long FindUniverseCell(long uni, double x, double y, double z, long *ridx, 
		      long id)
{
  long cell, lst, ptr, msh, found, loc0, loc1, loc2, n;

  /* Check universe pointer */

  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Check plotter mode */
  
  if (((long)RDB[DATA_PLOTTER_MODE] == NO) || 
      ((long)RDB[DATA_QUICK_PLOT_MODE] == YES))
    {
      /* Check previous */

      ptr = RDB[uni + UNIVERSE_PTR_PREV_REG];
      if ((lst = GetPrivateData(ptr, id)) > VALID_PTR)
	{
	  /* Get cell pointer */
	  
	  cell = (long)RDB[lst + CELL_LIST_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);      

	  /* Test cell */
	  
	  if (InCell(cell, x, y, z, NO, id) == YES)
	    {
	      /* Put region index */
	      
	      *ridx = (long)RDB[lst + CELL_LIST_REG_IDX];

	      /* Pointer to counter */

	      loc2 = (long)RDB[lst + CELL_LIST_PTR_COUNT];
	      CheckPointer(FUNCTION_NAME, "(loc2)", PRIVA_ARRAY, loc2);

	      /* Add counter */

	      AddPrivateData(loc2, 1, id);
	      
	      /* Return cell pointer */
	      
	      return cell;
	    }
	}
    }

  /* Check for search mesh */
  
  if ((msh = (long)RDB[uni + UNIVERSE_PTR_CELL_MESH]) > VALID_PTR)
    {
      /* Pointer to mesh */

      msh = (long)RDB[msh + CELL_MESH_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get pointer from search mesh */

      if ((ptr = MeshPtr(msh, x, y, z)) > VALID_PTR)
	ptr = (long)RDB[ptr];
      else
	ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
    }
  else
    ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
    
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
      
  /* Reset flag and cell pointer */

  found = NO;
  cell = GEOM_ERROR_NO_CELL;
  lst = -1;

  /* Loop over cells */
  
  n = 0;

  while ((loc0 = ListPtr(ptr, n++)) > VALID_PTR)
    {
      /* Pointer to cell */

      loc1 = (long)RDB[loc0 + CELL_LIST_PTR_CELL];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Test cell */

      if (InCell(loc1, x, y, z, NO, id) == YES)
	{
	  /* Put region index */

	  *ridx = (long)RDB[loc0 + CELL_LIST_REG_IDX];

	  /* Check plotter mode */

	  if ((long)RDB[DATA_PLOTTER_MODE] == NO)
	    {
	      /* Pointer to counter */

	      loc2 = (long)RDB[loc0 + CELL_LIST_PTR_COUNT];
	      CheckPointer(FUNCTION_NAME, "(loc2)", PRIVA_ARRAY, loc2);
	   
	      /* Add counter */

	      AddPrivateData(loc2, 1, id);

	      /* Put previous pointer */
	      
	      ptr = RDB[uni + UNIVERSE_PTR_PREV_REG];
	      PutPrivateData(ptr, loc0, id);
	      
	      /* Return cell pointer */

	      return loc1;
	    }
	  else if (found == YES)
	    {
	      /* Duplicate definition of geometry region */

	      return GEOM_ERROR_MULTIPLE_CELLS;
	    }
	  else
	    {
	      /* Put pointer */

	      cell = loc1;
	      lst = loc0;

	      /* Put found flag */

	      found = YES;
	    }
	}
    }

  /* Put previous pointer */
  
  ptr = RDB[uni + UNIVERSE_PTR_PREV_REG];
  PutPrivateData(ptr, lst, id);

  /* Return cell pointer */

  return cell;
}

/*****************************************************************************/
