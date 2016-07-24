#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findtetcell.c                                  */
/*                                                                           */
/* Created:       2012/09/11 (JLe)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Finds tetrahedral mesh cell                                  */
/*                                                                           */
/* Comments: -NOTE: Return value is pointer to tet cell, not geometry cell.  */
/*            (tuonne neighbour cell -rakenteeseen tallennetaan tässä        */
/*             moodissa tet-cellin pointteri, muulloin geometriacellin       */
/*             pointteri.)                                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindTetCell:"

#define FAST_MODE

/*****************************************************************************/

long FindTetCell(long ifc, double x, double y, double z, long id)
{
  long msh, lst, loc0, ptr;
  /*
  long nc, cell, loc1;
  double dx, dy, dz, r2, min;
  */

  /* Check pointer */
  
  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

#ifdef FAST_MODE

  /***************************************************************************/

  /***** Try previous cell ***************************************************/

  /* Check previous cell */

  ptr = (long)RDB[ifc + IFC_PTR_PREV_CELL];
  CheckPointer(FUNCTION_NAME, "(ptr1)", PRIVA_ARRAY, ptr);

  if ((loc0 = GetPrivateData(ptr, id)) > VALID_PTR)
    {

      /* Test cell */

      if (InTetCell(ifc, loc0, x, y, z, YES, id) == YES)
	{
	  /* Return pointer */

	  return loc0;
	}
    }

  /***************************************************************************/

  /***** Try neighbour search ************************************************/
      
  /* Check if previous cell was given (NOTE: Tää on nyt disabloitu kun  */
  /* se ei toimi kovin hyvin) */
  /* NOTE: Tää pitää kirjoittaa uudestaan uusien tetrojen vuoksi -Ville */

#ifdef taaeitoiminyt
  if (1 == 2)
  if (loc0 > VALID_PTR)
    {
      /* Reset minimum distance */

      min = INFTY;
      
      /* Reset counter */
      
      nc = 0;
      
      /* Loop */
      
      while (1 != 2)
	{
	  /* Check pointer */
	  
	  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
	  
	  /* Pointer to geometry cell */
	  
	  cell = (long)RDB[loc0 + IFC_TET_MSH_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
	  
	  /* Test cell */
	  
	  if (InCell(cell, x, y, z, YES, id) == YES)
	    {
	      /* Store pointer */
	      
	      ptr = (long)RDB[ifc + IFC_PTR_PREV_CELL];
	      CheckPointer(FUNCTION_NAME, "(ptr2)", PRIVA_ARRAY, ptr);
	      PutPrivateData(ptr, loc0, id);
	      
	      /* Return pointer */
	      
	      return loc0;
	    }
	  
	  /* Reset pointer */
	  
	  loc0 = -1;
	  
	  /* Pointer to intersection list */
	  
	  loc1 = (long)RDB[cell + CELL_PTR_SURF_INSC];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  
	  /* Loop over boundary surfaces using intersection list */
	  
	  while (loc1 > VALID_PTR)
	    {
	      /* Get pointer to neighbour cell */
	      
	      if ((ptr = (long)RDB[loc1 + CELL_INSC_PTR_NEXT_TET_CELL]) 
		  > VALID_PTR)
		{
		  /* Pointer to geometry cell */
		  
		  cell = (long)RDB[ptr + IFC_TET_MSH_PTR_CELL];
		  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
		  
		  /* Calculate square distance */
		  
		  dx = (x - RDB[ptr + CELL_CENTER_X]);
		  dy = (y - RDB[ptr + CELL_CENTER_Y]);
		  dz = (z - RDB[ptr + CELL_CENTER_Z]);
		  
		  /* Compare to minimum */
		  
		  if ((r2 = dx*dx + dy*dy + dz*dz) == 0.0)
		    Die(FUNCTION_NAME, "Error in vector lenght");
		  else if (r2 < min)
		    {
		      /* Update pointer and minimum */
		      
		      loc0 = ptr;
		      min = r2;		      
		    }
		}
	      
	      /* Next boundary surface */
	      
	      loc1 = NextItem(loc1);
	    }
	  
	  /* Check if minimum was found */
	  
	  if (loc0 < VALID_PTR)
	    {
	      /* No neighbour cell (point may be outside), break loop */
	      
	      break;
	    }
	  
	  /* Check maximum number of iterations */
	  
	  if (nc++ > 10)
	    break;
	}
    }
#endif

  /***************************************************************************/

  /***** Try list search *****************************************************/

  /* Get pointer to search mesh */

  msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
  
  /* Get pointer */
  
  if ((lst = MeshPtr(msh, x, y, z)) > VALID_PTR)
    lst = (long)RDB[lst];
  else
    return NULLPTR;

  /* Check pointer */
  
  if (lst < VALID_PTR)
    return NULLPTR;

  /* Loop over content */
      
  while (lst > VALID_PTR)
    {
      /* Pointer to cell */
	  
      loc0 = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
	  
      /* Check with bounding box */
	  
      if ((x < RDB[loc0 + IFC_TET_MSH_XMIN]) ||
	  (x > RDB[loc0 + IFC_TET_MSH_XMAX]) ||
	  (y < RDB[loc0 + IFC_TET_MSH_YMIN]) ||
	  (y > RDB[loc0 + IFC_TET_MSH_YMAX]) ||
	  (z < RDB[loc0 + IFC_TET_MSH_ZMIN]) ||
	  (z > RDB[loc0 + IFC_TET_MSH_ZMAX]))
	{
	  /* Cannot be inside, pointer to next */
	      
	  lst = NextItem(lst);
	  
	  /* Cycle loop */
	  
	  continue;
	}	      
      
      /* Test cell */

      if (InTetCell(ifc, loc0, x, y, z, YES, id) == YES)
	{
	  /* Check plotter mode */

	  if ((long)RDB[DATA_PLOTTER_MODE] == NO)
	    {
	      /* Pointer to counter */
	      
	      ptr = (long)RDB[lst + SEARCH_MESH_PTR_CELL_COUNT];
	      CheckPointer(FUNCTION_NAME, "(ptr3)", PRIVA_ARRAY, ptr);
	      		  
	      /* Add counter */
	      
	      AddPrivateData(ptr, 1, id);
	    }
	      
	  /* Store pointer */

	  ptr = (long)RDB[ifc + IFC_PTR_PREV_CELL];
	  CheckPointer(FUNCTION_NAME, "(ptr4)", PRIVA_ARRAY, ptr);
	  PutPrivateData(ptr, loc0, id);

	  /* Return pointer */

	  return loc0;
	}

      /* Next */
      
      lst = NextItem(lst);
      
      /***********************************************************************/
    }

  /***************************************************************************/

#else 

  /***************************************************************************/
  
  /***** Safe mode (test all cells) ******************************************/

  /* Loop over all tet cells */
      
  loc0 = (long)RDB[ifc + IFC_PTR_TET_MSH];
  while (loc0 > VALID_PTR)
    {
      /* Check with bounding box */

      if ((x < RDB[loc0 + IFC_TET_MSH_XMIN]) ||
	  (x > RDB[loc0 + IFC_TET_MSH_XMAX]) ||
	  (y < RDB[loc0 + IFC_TET_MSH_YMIN]) ||
	  (y > RDB[loc0 + IFC_TET_MSH_YMAX]) ||
	  (z < RDB[loc0 + IFC_TET_MSH_ZMIN]) ||
	  (z > RDB[loc0 + IFC_TET_MSH_ZMAX]))
	{
	  /* Cannot be inside, pointer to next */
	  
	  loc0 = NextItem(loc0);
	  
	  /* Cycle loop */
	  
	  continue;
	}	      
      
      /* Test */
      
      if (InTetCell(ifc, loc0, x, y, z, YES, id) == YES)
	return loc0;
      
      /* Next tet cell */
      
      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

#endif

  /* Not in any, return null */

  return NULLPTR;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
