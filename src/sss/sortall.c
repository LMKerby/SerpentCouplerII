#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sortall.c                                      */
/*                                                                           */
/* Created:       2011/03/01 (JLe)                                           */
/* Last modified: 2015/10/31 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Sorts reaction, surface etc. lists to speed up calculation   */
/*                                                                           */
/* Comments: - Reaction list sorting doesn't work in reproduceable MPI mode  */
/*           - Sorting the tet mesh lists is expensive if the search         */
/*             mesh is large (not in use at the moment).                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SortAll:"

/*****************************************************************************/

void SortAll()
{
  long mat, cell, nuc, lst, idx, uni, icm, loc0, msh, nx, ny, nz, i, j, k;
  
  /* Get cycle index */

  idx = (long)RDB[DATA_CYCLE_IDX];

  /* Check condition */
  
  /*
  if (((idx > 20) && (idx % 20)) || ((idx > 5) && (idx % 5)))
    return;
  */

  if (idx != (long)RDB[DATA_SORT_COUNT])
    return;
  else
    WDB[DATA_SORT_COUNT] = RDB[DATA_SORT_COUNT]*3.0;

  /***************************************************************************/

  /***** Reaction lists in nuclides and materials ****************************/

  /* Check MPI option */

  if ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO)
    {
      /* Loop over nuclides */
      
      nuc = (long)RDB[DATA_PTR_NUC0];
      while (nuc > VALID_PTR)
	{
	  /* Check mode */
	  
	  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
	    {
	      /* Pointer to partial reaction list */
	      
	      lst = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
	      CheckPointer(FUNCTION_NAME, "(lst1)", DATA_ARRAY, lst);
	      
	      /* Sort list */
	      
	      SortList(lst, RLS_DATA_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
	    }
	  
	  /* Next nuclide */
	  
	  nuc = NextItem(nuc);
	}
      
      /* Loop over materials */
      
      mat = (long)RDB[DATA_PTR_M0];
      while(mat > VALID_PTR)
	{
	  /* Sort material-wise list of totals */
	  
	  if ((lst = (long)RDB[mat + MATERIAL_PTR_TOT_REA_LIST]) 
	      > VALID_PTR)
	    SortList(lst, RLS_DATA_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);

	  if ((lst = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANT_LIST]) 
	      > VALID_PTR)
	    SortList(lst, RLS_DATA_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);

	  /* Next material */
	  
	  mat = NextItem(mat);
	}
    }

  /***************************************************************************/

  /***** Intersection lists for cells ****************************************/

  /* Loop over cells */

  cell = (long)RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Pointer to intersection list */
      
      if ((lst = (long)RDB[cell + CELL_PTR_SURF_INSC]) > VALID_PTR)
	{
	  /* Sort list */
	  
	  SortList(lst, CELL_INSC_PTR_OUT_COUNT, SORT_MODE_DESCEND_PRIVA);
	}

      /* Next cell */

      cell = NextItem(cell);
    }

  /***************************************************************************/

  /***** Cell lists for universes ********************************************/

  /* Loop over universes */

  uni = (long)RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Pointer to cell list (pointteri on PRIVA-blokkiin) */
      
      if ((lst = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST]) > 0)
	{
	  /* Sort list */

	  SortList(lst, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
	}

      /* Next universe */

      uni = NextItem(uni);
    }

  /***************************************************************************/

  /***** Cell search meshes **************************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_CSM0];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to mesh */

      msh = (long)RDB[loc0 + CELL_MESH_PTR_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];
      
      /* Check size */

      CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);
      
      /* Loop over mesh */
      
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	  for (k = 0; k < nz; k++)
	    {
	      /* Get pointer to list */
	      
	      lst = ReadMeshPtr(msh, i, j, k);
	      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

	      lst = (long)RDB[lst];
	      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);
  
	      /* Sort list */
	      
	      SortList(lst, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
	    }

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** STL solid search meshes *********************************************/

  /* Loop over structures */

  loc0 = (long)RDB[DATA_PTR_STL0];
  while (loc0 > VALID_PTR)
    {
      /* Get pointer to mesh */

      msh = (long)RDB[loc0 + STL_PTR_SOLID_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];
      
      /* Check size */

      CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);
      
      /* Loop over mesh */
      
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	  for (k = 0; k < nz; k++)
	    {
	      /* Get pointer to list */
	      
	      lst = ReadMeshPtr(msh, i, j, k);
	      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

	      /* Sort list */

	      if ((lst = (long)RDB[lst]) > VALID_PTR)
		SortList(lst, CELL_LIST_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
	    }

      /* Pointer to next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Tet mesh search lists ***********************************************/

  /* Loop over interfaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];
  while (loc0 > VALID_PTR)  
    {
      /* Check type */

      if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH)
	{
	  /* Get pointer to search mesh */
      
	  msh = (long)RDB[loc0 + IFC_PTR_SEARCH_MESH];
	  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
	  
	  /* Get size */

	  nx = (long)RDB[msh + MESH_N0];
	  ny = (long)RDB[msh + MESH_N1];
	  nz = (long)RDB[msh + MESH_N2];

	  /* Check size */

	  CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);

	  /* Loop over mesh */

	  for (i = 0; i < nx; i++)
	    for (j = 0; j < ny; j++)
	      for (k = 0; k < nz; k++)
		{
		  /* Get pointer to list */
	  
		  lst = ReadMeshPtr(msh, i, j, k);
		  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

		  /* Check content and sort list */
	
		  if ((lst = (long)RDB[lst]) > VALID_PTR)
		    if (ListSize(lst) > 10000)
		      SortList(lst, SEARCH_MESH_PTR_CELL_COUNT, 
			       SORT_MODE_DESCEND_PRIVA);
		}
	}

      /* Next interface */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Unstructured mesh search lists **************************************/

#ifdef mmmmmmmmmmmmmmmmmmmmmm

  /* Loop over geometries */

  loc0 = (long)RDB[DATA_PTR_UMSH0];
  while (loc0 > VALID_PTR)  
    {
      /* Get pointer to search mesh */
      
      msh = (long)RDB[loc0 + UMSH_PTR_SEARCH_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
	  
      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];
      
      /* Check size */
      
      CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);
      
      /* Loop over mesh */
      
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	  for (k = 0; k < nz; k++)
	    {
	      /* Get pointer to list */
	      
	      lst = ReadMeshPtr(msh, i, j, k);
	      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);
	      
	      /* Check content and sort list */
	      
	      if ((lst = (long)RDB[lst]) > VALID_PTR)
		if (ListSize(lst) > 10000)
		  SortList(lst, SEARCH_MESH_PTR_CELL_COUNT, 
			   SORT_MODE_DESCEND_PRIVA);
	    }
      
      /* Next geometry */

      loc0 = NextItem(loc0);
    }

#endif

  /***************************************************************************/

  /***** Stochastic geometry search lists ************************************/
  
  /* Loop over geometries */

  loc0 = (long)RDB[DATA_PTR_PB0];
  while (loc0 > VALID_PTR)  
    {
      /* Pointer to search mesh */

      msh = (long)RDB[loc0 + PBED_PTR_SEARCH_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
      
      /* Get size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];
      
      /* Check size */
      
      CheckValue(FUNCTION_NAME, "size", "", nx*ny*nz, 1, 1000000000);
      
      /* Loop over mesh */
      
      for (i = 0; i < nx; i++)
	for (j = 0; j < ny; j++)
	  for (k = 0; k < nz; k++)
	    {
	      /* Get pointer to list */
	      
	      lst = ReadMeshPtr(msh, i, j, k);
	      CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

	      /* Check content and sort list */
	      
	      if ((lst = (long)RDB[lst]) > VALID_PTR)
		SortList(lst, SEARCH_MESH_PTR_CELL_COUNT, 
			 SORT_MODE_DESCEND_PRIVA);
	    }
      
      /* Next geometry */
      
      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Misc. stuff *********************************************************/

  /* ICM break list */

  if ((icm = (long)RDB[DATA_PTR_ICM0]) > VALID_PTR)
    SortList(icm, ICM_BREAK_PTR_COUNT, SORT_MODE_DESCEND_PRIVA);
  
  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
