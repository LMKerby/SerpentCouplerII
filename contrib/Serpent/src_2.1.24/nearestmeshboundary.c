/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nearestmeshboundary.c                          */
/*                                                                           */
/* Created:       2014/03/01 (JLe)                                           */
/* Last modified: 2014/12/26 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description:  Calculates distance nearest mesh boundary                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "NearestMeshBoundary:"

/*****************************************************************************/

double NearestMeshBoundary(long msh, double x, double y, double z,
			   double u, double v, double w, long *fail)
{
  long ptr, loc0, idx, nx, ny, nz, i, j, k, n;
  double xmin, xmax, ymin, ymax, zmin, zmax, min, l, px, py, pz, dx, dy, dz;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check type */

  if ((long)RDB[msh + MESH_TYPE] != MESH_TYPE_ADAPTIVE)
    Die(FUNCTION_NAME, "Invalid mesh type");

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_PTR)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Reset minimum distance */

  min = INFTY;

  /* Get mesh index */

  if ((idx = MeshIndex(msh, x, y, z)) < 0)
    {
      /* Point is outside the mesh, get boundaries */

      xmin = RDB[msh + MESH_MIN0];
      xmax = RDB[msh + MESH_MAX0];
      ymin = RDB[msh + MESH_MIN1];
      ymax = RDB[msh + MESH_MAX1];
      zmin = RDB[msh + MESH_MIN2];
      zmax = RDB[msh + MESH_MAX2];

      /* Check boundaries to be sure */

      if ((x > xmin) && (x < xmax) && (y > ymin) && (y < ymax) &&
	  (z > zmin) && (z < zmax))
	Warn(FUNCTION_NAME, "Point is inside mesh");

      /* Calculate distance to outer boundaries */

      if (u != 0.0)
	{
	  if ((l = -(x - xmin)/u) >= 0.0)
	    if (l < min)
	      min = l;
	  
	  if ((l = -(x - xmax)/u) >= 0.0)
	    if (l < min)
	      min = l;
	}
      
      if (v != 0.0)
	{
	  if ((l = -(y - ymin)/v) >= 0.0)
	    if (l < min)
	      min = l;
	  
	  if ((l = -(y - ymax)/v) >= 0.0)
	    if (l < min)
	      min = l;
	}
      
      if (w != 0.0)
	{      
	  if ((l = -(z - zmin)/w) >= 0.0)
	    if (l < min)
	      min = l;
	  
	  if ((l = -(z - zmax)/w) >= 0.0)
	    if (l < min)
	      min = l;
	}

      /* Check value */

      if (min < 0.0)
	Warn(FUNCTION_NAME, "min = %E\n", min);
      else
	{
	  /* Check upper limit */

	  CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);
	}

      /* Return distance */

      return min;
    }
    
  /* Get direct pointer to data */
  
  ptr = (long)RDB[msh + MESH_PTR_PTR] + idx;
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Check pointer */
  
  if ((loc0 = (long)RDB[ptr]) == NULLPTR)
    {
      min = NearestMeshBoundary(-loc0, x, y, z, u, v, w, fail);

      /* Tää voi olla ongelma */
      /*
      Die(FUNCTION_NAME, "WTF?");
      */
    }
  else if (loc0 < -VALID_PTR)
    {
      /* Pointer to new mesh, call recursively */

      min = NearestMeshBoundary(-loc0, x, y, z, u, v, w, fail);
    }
  else
    {
      /* Get mesh boundaries */
      
      xmin = RDB[msh + MESH_MIN0];
      xmax = RDB[msh + MESH_MAX0];
      ymin = RDB[msh + MESH_MIN1];
      ymax = RDB[msh + MESH_MAX1];
      zmin = RDB[msh + MESH_MIN2];
      zmax = RDB[msh + MESH_MAX2];

      /* Get mesh size */

      nx = (long)RDB[msh + MESH_N0];
      ny = (long)RDB[msh + MESH_N1];
      nz = (long)RDB[msh + MESH_N2];
      
      /* Calculate mesh pitch */
      
      px = (xmax - xmin)/((double)nx);
      py = (ymax - ymin)/((double)ny);
      pz = (zmax - zmin)/((double)nz);
      
      CheckValue(FUNCTION_NAME, "px", "", px, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "py", "", py, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "pz", "", pz, ZERO, INFTY);

      /* Calculate mesh indexes */
      
      i = (long)((x - xmin)/px);
      j = (long)((y - ymin)/py);
      k = (long)((z - zmin)/pz);
      
      /* Check limits */
      
      CheckValue(FUNCTION_NAME, "i", "", i, 0, nx - 1);
      CheckValue(FUNCTION_NAME, "j", "", j, 0, ny - 1);
      CheckValue(FUNCTION_NAME, "k", "", k, 0, nz - 1);
      
      /* Calculate distance to mesh walls */
      
      dx = x - xmin - (double)i*px;
      dy = y - ymin - (double)j*py;
      dz = z - zmin - (double)k*pz;
      
      /* Check values (NOTE: toi nolla pitää sallia tässä ja ottaa alemmissa */
      /* if-lauseissa huomioon yhtäsuuruutena, sillä desimaalipyöristyksen   */
      /* takia piste voi osua suoraan pinnalle) */

      CheckValue(FUNCTION_NAME, "dx", "", dx, -1E-12, px);
      CheckValue(FUNCTION_NAME, "dy", "", dy, -1E-12, py);
      CheckValue(FUNCTION_NAME, "dz", "", dz, -1E-12, pz);

      /* Calculate distance to mesh cell walls */

      if (u != 0.0)
	{
	  if ((l = -dx/u) >= 0.0)
	    {
	      if (l < min)
		min = l;
	    }
	  else if ((l = -(dx - px)/u) > 0.0)
	    {
	      if (l < min)
		min = l;
	    }
	  else
	    Warn(FUNCTION_NAME, "dx = %E, u = %E, l = %E", dx, u, l);
	}

      if (v != 0.0)
	{
	  if ((l = -dy/v) >= 0.0)
	    {
	      if (l < min)
		min = l;
	    }
	  else if ((l = -(dy - py)/v) > 0.0)
	    {
	      if (l < min)
		min = l;
	    }
	  else
	    Warn(FUNCTION_NAME, "dy = %E, v = %E, l = %E", dy, v, l);
	}

      if (w != 0.0)
	{
	  if ((l = -dz/w) >= 0.0)
	    {
	      if (l < min)
		min = l;
	    }
	  else if ((l = -(dz - pz)/w) > 0.0)
	    {
	      if (l < min)
		min = l;
	    }
	  else
	    Warn(FUNCTION_NAME, "dz = %E, w = %E, l = %E", dz, w, l);
	}

      /* Check fail condition for STLRayTest() */
  
      if (fail != NULL)
	{
	  /* Count vertices */
	  
	  n = 0;
	  
	  if (fabs(x + min*u - xmin) < 1E-6)
	    n++;
	  else if (fabs(x + min*u - xmax) < 1E-6)
	    n++;
	  
	  if (fabs(y + min*v - ymin) < 1E-6)
	    n++;
	  else if (fabs(y + min*v - ymax) < 1E-6)
	    n++;
	  
	  if (fabs(z + min*w - zmin) < 1E-6)
	    n++;
	  else if (fabs(z + min*w - zmax) < 1E-6)
	    n++;

	  /* Check and set fail flag */

	  if (n > 1)
	    *fail = YES;
	  else
	    *fail = NO;
	}
    }
  
  /* Check value */

  CheckValue(FUNCTION_NAME, "min", "", min, 0.0, INFTY);

  /* Return distance */
  
  return min;
}

/*****************************************************************************/
