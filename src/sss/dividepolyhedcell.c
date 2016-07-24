#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dividepolyhedcell.c                            */
/*                                                                           */
/* Created:       2014/02/24 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Splits a polyhedral cell into tetrahedrons by adding points  */
/*              to face centerpoints as well as to the centerpoint of        */
/*              face centerpoints ("cell centerpoint")                       */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DividePolyhedCell:"

/*****************************************************************************/

void DividePolyhedCell(long ifc, long idx, long *sdone)
{
  long surf, np, i, j, n, cgns;
  long facelist, surflist, sidelist, side;
  long ptr, *newcells, *tempcells;
  long ndone, nsides, nc;

  /* Initialize number of created subcells */

  ndone = 0;

  /* Get pointer to interface surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  nc = 0;

  /* Get pointer to parent list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PRNTS];
  CheckPointer(FUNCTION_NAME, "(cgnslist)", DATA_ARRAY, cgns);

  /* Get pointer to cell to be divided */

  cgns = ListPtr(cgns, idx);
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get number of faces */

  nsides = (long)RDB[cgns + IFC_TET_MSH_NF];

  /* Get pointer to face list */
  
  ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, ptr);

  /* Loop over faces to get number of sides and points */

  for (j = 0; j < nsides; j++)
    {
      /* Get index of face */
      
      n = (long)RDB[ptr + j];

      /* Get pointer to face */

      surf = ListPtr(surflist, n);

      /* Get number of points for face */

      np = (long)RDB[surf + SURFACE_N_PARAMS];

      /* Add to number of cells to be created */

      nc += np;
    }

  /*printf("Dividing cell %ld to %ld cells (%ld faces)\n", idx, nc, nsides);*/

  /* Allocate memory for temporary cell list */

  tempcells = (long *)Mem(MEM_ALLOC, nc, sizeof(long));

  /* Get pointer to face list */
  
  facelist = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, facelist);

  /* Get pointer to side list */
  
  sidelist = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, sidelist);

  /* Loop over faces to get number of sides and points */

  for (j = 0; j < nsides; j++)
    {
      /* Get index of face */
      
      n = (long)RDB[facelist + j];

      /* Get pointer to face */

      surf = ListPtr(surflist, n);

      /* Get side for this cell */

      side = (long)RDB[sidelist + j];

      /* Get points on the face */

      np = (long)RDB[surf + SURFACE_N_PARAMS];

      /* Allocate memory for a list of new face cells */
      /* DividePolyhedFace will use it to return a list of created cgns-cells*/

      newcells = (long *)Mem(MEM_ALLOC, np, sizeof(long));

      /* Divide face into np cells */

      DividePolyhedFace(ifc, idx, n, side, np, newcells, tempcells, 
			ndone, sdone);
      
      /* Copy pointers of new face cells to list of subcells for this cell */

      for (i = 0; i < np; i++)
	tempcells[ndone + i] = newcells[i];

      ndone += np;

      /* Free memory */

      Mem(MEM_FREE,newcells);
    }

  /* Copy parent cell pointer to new cells. Ville: Tarkistaisitko että tämä */
  /* menee oikein? Onko kaikki alkuperäisestä cellistä luodut cellit tuossa */
  /* vektorissa? */

  for (i = 0; i < nc; i++)
    {
      /* Get pointer to new cell structure */

      ptr = tempcells[i];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Put parent pointer */
      
      if (ptr != cgns)
	WDB[ptr + IFC_TET_MSH_PTR_PARENT] = (double)cgns;

    }

  /* Free memory */
    
  Mem(MEM_FREE,tempcells);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
