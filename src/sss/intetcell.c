/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : intetcell.c                                    */
/*                                                                           */
/* Created:       2010/10/12 (JLe)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Checks if point is inside a tet cell                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InTetCell:"

/*****************************************************************************/

long InTetCell(long ifc, long cgns, double x, double y, double z, long on, long id)
{
  long n, ptr, surf, side, facelist, surflist, sidelist, retval;
  long pt0, pt1, pt2, nf  ;
  long i;
  double x0, x1, x2, y0, y1, y2, z0, z1, z2, d;

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get pointer to interface surfaces */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to cell's face list */

  facelist = (long)RDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(facelist)", DATA_ARRAY, facelist);

  /* Get pointer to cell's side list */

  sidelist = (long)RDB[cgns + IFC_TET_MSH_PTR_SIDES];
  CheckPointer(FUNCTION_NAME, "(sidelist)", DATA_ARRAY, sidelist);

  /* Loop over cell faces */

  nf = (long)RDB[cgns + IFC_TET_MSH_NF];
  /*
  CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 4);
  */
  for (i = 0; i < nf; i++)
    {

      /* Get index of face */

      n = (long)RDB[facelist + i];

      /* Get side for this face */

      side = (long)RDB[sidelist + i];

      /* Get pointer to face surface */

      surf = ListPtr(surflist, n);

      /*      printf("face %ld, side %ld\n", n, side);*/

      /* Get pointer to surface parameters */

      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get pointers to points */

      pt0 = (long)RDB[ptr + 0];
      pt1 = (long)RDB[ptr + 1];
      pt2 = (long)RDB[ptr + 2];

      /* Three points, calculate vectors */

      x0 = -RDB[pt1 + 0] + x;
      y0 = -RDB[pt1 + 1] + y;
      z0 = -RDB[pt1 + 2] + z;

      x1 = -RDB[pt0 + 0] + RDB[pt1 + 0];
      y1 = -RDB[pt0 + 1] + RDB[pt1 + 1];
      z1 = -RDB[pt0 + 2] + RDB[pt1 + 2];

      x2 = -RDB[pt1 + 0] + RDB[pt2 + 0];
      y2 = -RDB[pt1 + 1] + RDB[pt2 + 1];
      z2 = -RDB[pt1 + 2] + RDB[pt2 + 2];

      /* Scalar triple product */

      d = x0*(y1*z2 - y2*z1) - y0*(x1*z2 - x2*z1) + z0*(x1*y2 - x2*y1);

      /* Check if surface itself is included as inside (NOTE: this is */
      /* used with tet-mesh structures to avoid undefined regions on  */
      /* cell boundaries, otherwise the boundary is included in the   */
      /* outside of the surface) */
      
      if (on == YES)
	{
	  /* Check */
	  
	  if (d <= 0.0)
	    retval = YES; 
	  else
	    retval = NO;
	}
      else
	{
	  /* Check */
	  
	  if (d < 0.0)
	    retval = YES; 
	  else
	    retval = NO;
	}

      /* Test surface */

      if (retval == YES)
	side = -side;

      /* Check result */

      if (side < 0)
	return NO;

    }

  /* Point is inside all surfaces */

  return YES;

}

/*****************************************************************************/
