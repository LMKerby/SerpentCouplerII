#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processvr.c                                    */
/*                                                                           */
/* Created:       2011/05/15 (JLe)                                           */
/* Last modified: 2015/10/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Processes some variance reduction stuff                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessVR:"

/*****************************************************************************/

void ProcessVR()
{
  long wwd, msh, loc0, loc1, ptr, nx, ny, nz, lat, i, j, k;
  double lims[6], p, q, x, y, z, min, max;

  /***************************************************************************/
  
  /****** Uniform fission source method **************************************/
  
  if ((long)RDB[DATA_UFS_MODE] != UFS_MODE_NONE)
    {
      /* Check pointer to lattice */

      if (((long)RDB[DATA_UFS_PTR_LAT]) > VALID_PTR)
	{
	  /* Find lattice */

	  lat = (long)RDB[DATA_PTR_L0];
	  while (lat > VALID_PTR)
	    {
	      /* Compare names */

	      if (CompareStr(lat + LAT_PTR_NAME, DATA_UFS_PTR_LAT))
		break;

	      /* Next lattice */

	      lat = NextItem(lat);
	    }

	  /* Check */

	  if (lat < VALID_PTR)
	    Error(0, "Lattice \"%s\" used for UFS not found in geometry",
		  GetText(DATA_UFS_PTR_LAT));
	  else
	    WDB[DATA_UFS_PTR_LAT] = (double)lat;
	  
	  /* Check type */

	  if (((long)RDB[lat + LAT_TYPE] != LAT_TYPE_S) && 
	      ((long)RDB[lat + LAT_TYPE] != LAT_TYPE_HX) && 
	      ((long)RDB[lat + LAT_TYPE] != LAT_TYPE_HY))
	    Error(0, "Lattice \"%s\" is wrong type for UFS", 
		  GetText(lat + LAT_PTR_NAME));

	  /* Put x- and y-bins */

	  WDB[DATA_UFS_NX] = RDB[lat + LAT_NX];
	  WDB[DATA_UFS_NY] = RDB[lat + LAT_NY];
	}
      
      /* Set boundaries if not set in input */
      
      if (RDB[DATA_UFS_XMIN] == -INFTY)
	WDB[DATA_UFS_XMIN] = RDB[DATA_GEOM_MINX];
      
      if (RDB[DATA_UFS_XMAX] == INFTY)
	WDB[DATA_UFS_XMAX] = RDB[DATA_GEOM_MAXX];
      
      if (RDB[DATA_UFS_YMIN] == -INFTY)
	WDB[DATA_UFS_YMIN] = RDB[DATA_GEOM_MINY];
      
      if (RDB[DATA_UFS_YMAX] == INFTY)
	WDB[DATA_UFS_YMAX] = RDB[DATA_GEOM_MAXY];
      
      if (RDB[DATA_UFS_ZMIN] == -INFTY)
	WDB[DATA_UFS_ZMIN] = RDB[DATA_GEOM_MINZ];
      
      if (RDB[DATA_UFS_ZMAX] == INFTY)
	WDB[DATA_UFS_ZMAX] = RDB[DATA_GEOM_MAXZ];
      
      /* Get mesh size */
      
      nx = (long)RDB[DATA_UFS_NX];
      ny = (long)RDB[DATA_UFS_NY];
      nz = (long)RDB[DATA_UFS_NZ];

      /* Check values */
      
      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 10000);
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 10000);
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 10000);
      
      /* Put mesh variables */
      
      lims[0] = RDB[DATA_UFS_XMIN];
      lims[1] = RDB[DATA_UFS_XMAX];
      lims[2] = RDB[DATA_UFS_YMIN];
      lims[3] = RDB[DATA_UFS_YMAX];
      lims[4] = RDB[DATA_UFS_ZMIN];
      lims[5] = RDB[DATA_UFS_ZMAX];
      
      /* Create mesh */
      
      ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_RES, -1,
		       nx, ny, nz, lims, -1);
      
      /* Put pointer */
      
      WDB[DATA_UFS_PTR_SRC_MESH] = (double)ptr;	      
      
      /* Allocate memory for factors */
      
      ptr = ReallocMem(DATA_ARRAY, nx*ny*nz);
      WDB[DATA_UFS_PTR_FACTORS] = (double)ptr;
    }

  /***************************************************************************/

  /***** Process weight windows **********************************************/

  /* Loop over weight windows */

  wwd = (long)RDB[DATA_PTR_WWD0];
  while (wwd > VALID_PTR)
    {
      /* Check normalization */

      if (RDB[wwd + WWD_NORM_FACT] > 0.0)
	{
	  /* Get point */
	  
	  x = RDB[wwd + WWD_NORM_X];
	  y = RDB[wwd + WWD_NORM_Y];
	  z = RDB[wwd + WWD_NORM_Z];

	  /* Pointer to mesh */
      
	  msh = (long)RDB[wwd + WWD_PTR_MESH];
	  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
	  	  
	  /* Pointer to data */
	  
	  if ((ptr = MeshPtr(msh, x, y, z)) < VALID_PTR)
	    Error(wwd, "Normalization point is outside mesh");
	  
	  ptr = (long)RDB[ptr];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  /* Get value */
	  
	  p = RDB[ptr + WWD_MESH_IMP];
	  
	  /* Get exponential */
	  
	  if ((q = RDB[wwd + WWD_POW]) > 0.0)
	    p = powl(p, q);
	  
	  /* Calculate normalization factor */
	  
	  WDB[wwd + WWD_NORM_FACT] = RDB[wwd + WWD_NORM_FACT]/p;
	}

      /* Get pointer to mesh */
      
      if ((msh = (long)RDB[wwd + WWD_PTR_MESH]) > VALID_PTR)
	{
	  /* Get size */

	  nx = (long)RDB[msh + MESH_N0];
	  ny = (long)RDB[msh + MESH_N1];
	  nz = (long)RDB[msh + MESH_N2];
	  
	  /* Reset minimum and maximum */
	  
	  min = INFTY;
	  max = -INFTY;
	  
	  /* Loop over mesh */
	  
	  for (i = 0; i < nx; i++)
	    for (j = 0; j < ny; j++)
	      for (k = 0; k < nz; k++)
		{
		  /* Pointer to data */
		  
		  loc0 = ReadMeshPtr(msh, i, j, k);
		  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
		  
		  loc0 = (long)RDB[loc0];
		  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
		  
		  /* Get value */
		  
		  p = RDB[loc0 + WWD_MESH_IMP];
		  
		  /* Compare to minimum and maximum */
		  
		  if (p > max)
		    max = p;
		  if (p < min)
		    min = p;
		  
		  /* Neighbour cells */
		  
		  if (i > 0)
		    {
		      /* Pointer to data */
		      
		      loc1 = ReadMeshPtr(msh, i - 1, j, k);
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      loc1 = (long)RDB[loc1];
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      /* Create current */
		      
		      ptr = NewItem(loc0 + WWD_MESH_PTR_CURR,
				    WWD_MESH_CURR_BLOCK_SIZE);
		      
		      /* Put pointer */
		      
		      WDB[ptr + WWD_MESH_CURR_PTR_NEIGHBOUR] = (double)loc1;
		    }
		  
		  if (i < nx - 1)
		    {
		      /* Pointer to data */
		      
		      loc1 = ReadMeshPtr(msh, i + 1, j, k);
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      loc1 = (long)RDB[loc1];
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      /* Create current */
		      
		      ptr = NewItem(loc0 + WWD_MESH_PTR_CURR,
				    WWD_MESH_CURR_BLOCK_SIZE);
		      
		      /* Put pointer */
		      
		      WDB[ptr + WWD_MESH_CURR_PTR_NEIGHBOUR] = (double)loc1;
		    }
		  
		  if (j > 0)
		    {
		      /* Pointer to data */
		      
		      loc1 = ReadMeshPtr(msh, i, j - 1, k);
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      loc1 = (long)RDB[loc1];
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      /* Create current */
		      
		      ptr = NewItem(loc0 + WWD_MESH_PTR_CURR,
				    WWD_MESH_CURR_BLOCK_SIZE);
		      
		      /* Put pointer */
		      
		      WDB[ptr + WWD_MESH_CURR_PTR_NEIGHBOUR] = (double)loc1;
		    }
		  
		  if (j < ny - 1)
		    {
		      /* Pointer to data */
		      
		      loc1 = ReadMeshPtr(msh, i, j + 1, k);
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      loc1 = (long)RDB[loc1];
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      /* Create current */
		      
		      ptr = NewItem(loc0 + WWD_MESH_PTR_CURR,
				    WWD_MESH_CURR_BLOCK_SIZE);
		      
		      /* Put pointer */
		      
		      WDB[ptr + WWD_MESH_CURR_PTR_NEIGHBOUR] = (double)loc1;
		    }
		  
		  if (k > 0)
		    {
		      /* Pointer to data */
		      
		      loc1 = ReadMeshPtr(msh, i, j, k - 1);
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      loc1 = (long)RDB[loc1];
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      /* Create current */
		      
		      ptr = NewItem(loc0 + WWD_MESH_PTR_CURR,
				    WWD_MESH_CURR_BLOCK_SIZE);
		      
		      /* Put pointer */
		      
		      WDB[ptr + WWD_MESH_CURR_PTR_NEIGHBOUR] = (double)loc1;
		    }
		  
		  if (k < nz - 1)
		    {
		      /* Pointer to data */
		      
		      loc1 = ReadMeshPtr(msh, i, j, k + 1);
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      loc1 = (long)RDB[loc1];
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      /* Create current */
		      
		      ptr = NewItem(loc0 + WWD_MESH_PTR_CURR,
				    WWD_MESH_CURR_BLOCK_SIZE);
		      
		      /* Put pointer */
		      
		      WDB[ptr + WWD_MESH_CURR_PTR_NEIGHBOUR] = (double)loc1;
		    }
		}
	  
	  /* Put values */

	  WDB[wwd + WWD_MESH_MIN] = min;
	  WDB[wwd + WWD_MESH_MAX] = max;
	}

      /* Next */

      wwd = NextItem(wwd);
    }

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
