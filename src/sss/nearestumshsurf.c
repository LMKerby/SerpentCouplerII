#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nearestumshsurf.c                              */
/*                                                                           */
/* Created:       2013/11/23 (JLe)                                           */
/* Last modified: 2015/04/29 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Calculates distance to nearest boundary when not inside      */
/*              mesh cell.                                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NearestUMSHSurf:"

/*****************************************************************************/

double NearestUMSHSurf(long ifc, double x, double y, double z, 
		       double u, double v, double w, long id)
{
  long loc0, msh, lst, surf, ptr, out;
  long i, j, k, n, nf, np, pt;
  long surflist, facelist;
  double xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, l, min, d;
  double x1, y1, z1, params[9];

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Get pointer to search mesh */

  msh = (long)RDB[ifc + IFC_PTR_SEARCH_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get pointer to interface surfaces */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /***************************************************************************/

  /***** Distance to outer boundaries ****************************************/

  /* Get mesh boundaries */

  xmin = RDB[msh + MESH_MIN0];
  xmax = RDB[msh + MESH_MAX0];
  ymin = RDB[msh + MESH_MIN1];
  ymax = RDB[msh + MESH_MAX1];
  zmin = RDB[msh + MESH_MIN2];
  zmax = RDB[msh + MESH_MAX2];
  
  /* Check that particle is in mesh */

  if ((x < xmin) || (x > xmax) || (y < ymin) || (y > ymax) || (z < zmin) ||
      (z > zmax))
    {
      /* Reset minimum distance */

      min = INFTY;

      /* Particle is not in mesh. Calculate distance to outer boundaries */

      if (u != 0.0)
	{
	  if ((l = -(x - xmin)/u) > 0.0)
	    if (l < min)
	      min = l;

	  if ((l = -(x - xmax)/u) > 0.0)
	    if (l < min)
	      min = l;
	}

      if (v != 0.0)
	{
	  if ((l = -(y - ymin)/v) > 0.0)
	    if (l < min)
	      min = l;
	  
	  if ((l = -(y - ymax)/v) > 0.0)
	    if (l < min)
	      min = l;
	}

      if (w != 0.0)
	{
	  if ((l = -(z - zmin)/w) > 0.0)
	    if (l < min)
	      min = l;
	  
	  if ((l = -(z - zmax)/w) > 0.0)
	    if (l < min)
	      min = l;
	}

      /* Check value */

      CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);

      /* Do zero cut-off */

      if (min < ZERO)
	min = ZERO;

      /* Return minimum distance */

      return min;
    }

  /***************************************************************************/

  /***** Point is in search mesh *********************************************/

  /* Distance to search mesh boundaries */

  min = NearestMeshBoundary(msh, x, y, z, u, v, w, NULL);

  /* Get pointer to search mesh (NOTE: tossa luupataan */
  /* toisen kerran adaptiivisen meshin yli) */

  lst = MeshPtr(msh, x, y, z);
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);
  lst = (long)RDB[lst];

  /* Loop over content */

  while (lst > VALID_PTR)
    {
      /* Pointer to tet cell */
      
      loc0 = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /***********************************************************************/

      /***** Distance to bounding box ****************************************/

      /* Reset out flag */

      out = NO;

      /* Get box boundaries */

      xmin = RDB[loc0 + IFC_TET_MSH_XMIN];
      xmax = RDB[loc0 + IFC_TET_MSH_XMAX];
      ymin = RDB[loc0 + IFC_TET_MSH_YMIN];
      ymax = RDB[loc0 + IFC_TET_MSH_YMAX];
      zmin = RDB[loc0 + IFC_TET_MSH_ZMIN];
      zmax = RDB[loc0 + IFC_TET_MSH_ZMAX];
      
      /* Check with bounding box */
		
      if ((dx = x - xmin) < 0.0)
	{
	  if (u != 0.0)
	    if (((d = -dx/u) > 0.0) && (d < min))
	      {
		/* Move to position */

		x1 = x + (d + EXTRAP_L)*u;
		y1 = y + (d + EXTRAP_L)*v;
		z1 = z + (d + EXTRAP_L)*w;
		
		/* Check if inside */

		if ((x1 > xmin) && (x1 < xmax) && (y1 > ymin) && (y1 < ymax) &&
		    (z1 > zmin) && (z1 < zmax))
		  min = d;
	      }

	  out = YES;
	}
      else if ((dx = x - xmax) > 0.0)
	{
	  if (u != 0.0)
	    if (((d = -dx/u) > 0.0) && (d < min))
	      {
		/* Move to position */

		x1 = x + (d + EXTRAP_L)*u;
		y1 = y + (d + EXTRAP_L)*v;
		z1 = z + (d + EXTRAP_L)*w;
		
		/* Check if inside */

		if ((x1 > xmin) && (x1 < xmax) && (y1 > ymin) && (y1 < ymax) &&
		    (z1 > zmin) && (z1 < zmax))
		  min = d;
	      }

	  out = YES;
	}

      if ((dy = y - ymin) < 0.0)
	{
	  if (v != 0.0)
	    if (((d = -dy/v) > 0.0) && (d < min))
	      {
		/* Move to position */

		x1 = x + (d + EXTRAP_L)*u;
		y1 = y + (d + EXTRAP_L)*v;
		z1 = z + (d + EXTRAP_L)*w;
		
		/* Check if inside */

		if ((x1 > xmin) && (x1 < xmax) && (y1 > ymin) && (y1 < ymax) &&
		    (z1 > zmin) && (z1 < zmax))
		  min = d;
	      }

	  out = YES;
	}
      else if ((dy = y - ymax) > 0.0)
	{
	  if (v != 0.0)
	    if (((d = -dy/v) > 0.0) && (d < min))
	      {
		/* Move to position */

		x1 = x + (d + EXTRAP_L)*u;
		y1 = y + (d + EXTRAP_L)*v;
		z1 = z + (d + EXTRAP_L)*w;
		
		/* Check if inside */

		if ((x1 > xmin) && (x1 < xmax) && (y1 > ymin) && (y1 < ymax) &&
		    (z1 > zmin) && (z1 < zmax))
		  min = d;
	      }

	  out = YES;
	}

      if ((dz = z - zmin) < 0.0)
	{
	  if (w != 0.0)
	    if (((d = -dz/w) > 0.0) && (d < min))
	      {
		/* Move to position */

		x1 = x + (d + EXTRAP_L)*u;
		y1 = y + (d + EXTRAP_L)*v;
		z1 = z + (d + EXTRAP_L)*w;
		
		/* Check if inside */

		if ((x1 > xmin) && (x1 < xmax) && (y1 > ymin) && (y1 < ymax) &&
		    (z1 > zmin) && (z1 < zmax))
		  min = d;
	      }

	  out = YES;
	}
      else if ((dz = z - zmax) > 0.0)
	{
	  if (w != 0.0)
	    if (((d = -dz/w) > 0.0) && (d < min))
	      {
		/* Move to position */

		x1 = x + (d + EXTRAP_L)*u;
		y1 = y + (d + EXTRAP_L)*v;
		z1 = z + (d + EXTRAP_L)*w;
		
		/* Check if inside */

		if ((x1 > xmin) && (x1 < xmax) && (y1 > ymin) && (y1 < ymax) &&
		    (z1 > zmin) && (z1 < zmax))
		  min = d;
	      }

	  out = YES;
	}

      /***********************************************************************/

      /***** Distance to cell faces ******************************************/

      /* Check if inside bounding box */

      if (out == NO)
	{
	  /* Get pointer to cell's face list */

	  facelist = (long)RDB[loc0 + IFC_TET_MSH_PTR_FACES];
	  CheckPointer(FUNCTION_NAME, "(facelist)", DATA_ARRAY, facelist);

	  /* Loop over cell faces */

	  nf = (long)RDB[loc0 + IFC_TET_MSH_NF];
	  CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 4);

	  for (i = 0; i < nf; i++)
	    {

	      /* Get index of face */

	      n = (long)RDB[facelist + i];

	      /* Get pointer to face surface */

	      surf = ListPtr(surflist, n);

	      /* Get pointer to surface parameters */

	      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
	      /* Get number of points on the face */
      
	      np = (long)RDB[surf + SURFACE_N_PARAMS];
	      CheckValue(FUNCTION_NAME, "np", "", np, 3, 3);

	      /* Copy points to params */

	      k = 0;

	      for (j = 0; j < np; j++)
		{
		  /* Get pointer to beginning of point coordinates */

		  pt = (long)RDB[ptr + j];
		
		  /* Store coordinates to params */

		  params[k++] = RDB[pt + 0];
		  params[k++] = RDB[pt + 1];
		  params[k++] = RDB[pt + 2];
	
		}

	      /* Get distance */
	      
	      d = SurfaceDistance(surf, params, SURF_PLANE, 9, 
				  x, y, z, u, v, w, id);
	      
	      /* Compare to minimum */
	      
	      if (d < min)
		min = d;

	      /* Next face */

	    }

	}

      /* Next */

      lst = NextItem(lst);
    }

  /* Check minimum distance */

  CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);

  /* Do zero cut-off */

  if (min < ZERO)
    min = ZERO;
  
  /* Return minimum */

  return min;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
