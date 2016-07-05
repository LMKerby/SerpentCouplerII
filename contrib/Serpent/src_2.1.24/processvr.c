/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processvr.c                                    */
/*                                                                           */
/* Created:       2011/05/15 (JLe)                                           */
/* Last modified: 2014/05/23 (JLe)                                           */
/* Version:       2.1.21                                                     */
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
  long ptr, nx, ny, nz, lat;
  double lims[6];

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
      
      ptr = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_DATA, -1,
		       nx, ny, nz, lims, -1);
      
      /* Put pointer */
      
      WDB[DATA_UFS_PTR_SRC_MESH] = (double)ptr;	      
      
      /* Allocate memory for factors */
      
      ptr = ReallocMem(DATA_ARRAY, nx*ny*nz);
      WDB[DATA_UFS_PTR_FACTORS] = (double)ptr;
    }

  /***************************************************************************/
}

/*****************************************************************************/
