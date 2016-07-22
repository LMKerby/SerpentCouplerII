/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fixpolyhedmesh.c                               */
/*                                                                           */
/* Created:       2014/02/24 (VVa)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Splits polyhedral meshes to tetrahedrons                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FixPolyhedMesh:"

/*****************************************************************************/

void FixPolyhedMesh(long ifc)
{
  long cgns, loc1, sdone;
  long ptr, ptr1, surflist, facepts, cellpts, surf;
  long solist, snlist;
  long nf, nc, np, i, j, n;  
  double newfaces;

  /* Get number of parent cells */

  nc = (long)RDB[ifc + IFC_NC_PRNTS];

  /* Allocate memory for cell centerpoints */

  cellpts = ReallocMem(DATA_ARRAY, nc*3);

  /* Store cell centerpoint list */

  WDB[ifc + IFC_PTR_PRNT_CELL_CP_LIST] = (double)cellpts;

  /* Get number of parent faces */

  nf = (long)RDB[ifc + IFC_NF_PRNTS];

  /* Allocate memory for face centerpoints */

  facepts = ReallocMem(DATA_ARRAY, nf*3);

  /* Store face centerpoint list */

  WDB[ifc + IFC_PTR_PRNT_FACE_CP_LIST] = (double)facepts;

  /* Get pointer to surfacelist */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PRNTS];
  CheckPointer(FUNCTION_NAME, "(cgns00)", DATA_ARRAY, cgns);  

  /* Calculate cell and face centerpoints */

  CalculateTetCenter(cgns, surflist, cellpts, facepts);

  /******* Calculate number of new faces ********************/

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PRNTS];
  CheckPointer(FUNCTION_NAME, "(cgns0)", DATA_ARRAY, cgns);

  newfaces = 0;

  /* Loop over cells to calculate number of new faces */

  for(i = 0; i < nc; i++)
    {
      /* Get pointer to cell */

      cgns = ListPtr(cgns, i);
      CheckPointer(FUNCTION_NAME, "(cgns1)", DATA_ARRAY, cgns);

      /* Get number of faces */

      nf = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Get pointer to face list */

      ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, ptr);

      /* Get pointer to side list */

      ptr1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
      CheckPointer(FUNCTION_NAME, "(PTR1)", DATA_ARRAY, ptr1);

      /* Loop over faces */

      for (j = 0; j < nf; j++)
	{
	  /* Get index of face */

	  n = (long)RDB[ptr + j];

	  /* Get pointer to surface */

	  surf = ListPtr(surflist, n);
	  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

	  /* Get number of points on the face */

	  np=(long)RDB[surf + SURFACE_N_PARAMS];
	  
	  /* Add to number of external faces */

	  if ((long)RDB[ptr1 + j] == -1)
	    {

	      /* This cell owns the child faces */
	      /* Sides out of the cell */

	      newfaces += np;

	    }

	  /* The sides between the cells on this face*/

	  newfaces += np;

	  /* The sides between this faces cells and other faces in this cell */

	  newfaces += np/2.0;

	}

    }

  /* Check that there are no half faces */

  if(newfaces - (long)newfaces != 0)
    Die(FUNCTION_NAME, "Something wrong with number of new faces");

  /* Allocate memory for new owners and neighbors list    */
  /* The number of child faces has to be calculated first */
   
  /* Allocate memory for pointers */
      
  solist = ReallocMem(DATA_ARRAY, (long)newfaces);

  WDB[ifc + IFC_PTR_OWNR_LIST] = (double)solist;

  snlist = ReallocMem(DATA_ARRAY, (long)newfaces);

  WDB[ifc + IFC_PTR_NBR_LIST] = (double)snlist;

  /* Loop over nbr/ownrlists to reset cells */

  for (i = 0; i < (long)newfaces; i++)
    {
      /* Create new surface */

      surf = NewItem(ifc + IFC_PTR_SURF_LIST, SURFACE_BLOCK_SIZE);

      /* Allocate memory for parameters */
		  
      ptr = ReallocMem(DATA_ARRAY, 3);
      WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

      /* Store number of parameters */

      WDB[surf + SURFACE_N_PARAMS] = (double)3.0;

    }

  surf = (long)RDB[ifc + IFC_PTR_SURF_LIST];

  /* Close new surface list */

  CloseList(surf);

  /* Get pointer to cells */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PRNTS];
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Reset number of created surfaces (needed in one of the subroutines) */

  sdone = 0;

  /* Loop over parents and divide the cells */

  for(i = 0; i < nc; i++)
    DividePolyhedCell(ifc, i, &sdone);
    
  /* Get pointer to new cells */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Close list */

  CloseList(cgns);

  /* Store number of child cells */

  WDB[ifc + IFC_NC] = (double)ListSize(cgns);

  /* Loop over new cell list to calculate cell centerpoints */

  nc = (long)RDB[ifc + IFC_NC];

  /*********** Remove minus signs from face list ***********/

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get pointer to surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Loop over cells to create centerpoints */

  for(i = 0; i < nc; i++)
    {

      /* Get pointer to cell */

      cgns = ListPtr(cgns, i);
      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      /* Get number of faces */

      nf = (long)RDB[cgns + IFC_TET_MSH_NF];

      if(nf != 4)
	Die(FUNCTION_NAME, 
	    "More than four faces on a tet-cell?");

      /* Get pointer to face list */
      
      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      for (j = 0; j < nf; j++)
	{

	  /* Get index of face */

	  n = (long)RDB[loc1 + j];

	  /* Remove minus sign from face index */

	  if (n <= 0)
	    WDB[loc1 + j] = (double)(-n);
	  else
	    Die(FUNCTION_NAME, 
		"Positive face index on child cell");

	}

    }

  /* Allocate memory for cell centerpoints */

  cellpts = ReallocMem(DATA_ARRAY, nc*3);

  /* Store cell centerpoint list */

  WDB[ifc + IFC_PTR_CELL_CP_LIST] = (double)cellpts;

  /* Get pointer to tet list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get pointer to surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Calculate cell centerpoints for new cells */

  CalculateTetCenter(cgns, surflist, cellpts, -1);

}

/*****************************************************************************/
