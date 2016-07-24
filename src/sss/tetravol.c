#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : tetravol.c                                     */
/*                                                                           */
/* Created:       2015/02/19 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Calculates volume of a tetrahedral unstructured mesh cell    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TetraVol:"

/*****************************************************************************/

double TetraVol(long ifc, long cgns)
{
  long surf, ptr, n, p0, p1, p2, p3, nf;
  long surflist, facelist;
  double vol, a1, a2, a3, b1, b2, b3;
  double c1, c2, c3;

  /* Check cell pointer */

  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get pointer to interface surfaces */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to cell's face list */

  facelist = (long)RDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(facelist)", DATA_ARRAY, facelist);

  /* Loop over cell faces */

  nf = (long)RDB[cgns + IFC_TET_MSH_NF];
  CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 4);

  /* Get index of face */

  n = (long)RDB[facelist + 0];

  /* Get pointer to face surface */

  surf = ListPtr(surflist, n);

  /* Get pointer to surface parameters */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
  /* Get pointers to points */

  p0 = (long)RDB[ptr + 0];
  p1 = (long)RDB[ptr + 1];
  p2 = (long)RDB[ptr + 2];

  /* Get fourth point from second face */

  /* Get index of face */

  n = (long)RDB[facelist + 1];

  /* Get pointer to face surface */

  surf = ListPtr(surflist, n);

  /* Get pointer to surface parameters */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get point candidate */

  p3 = (long)RDB[ptr + 0];

  /* Check if it is one of the first three */

  if ((p3 == p0) || (p3 == p1) || (p3 == p2))
    {
      /* Try next point */

      p3 = (long)RDB[ptr + 1];

      if ((p3 == p0) || (p3 == p1) || (p3 == p2))
	{
	  /* Try next point */

	  p3 = (long)RDB[ptr + 2];

	  if ((p3 == p0) || (p3 == p1) || (p3 == p2))
	    Die(FUNCTION_NAME, "Identical faces");
	  
	}

    }

  /* Calculate vector components */
  
  a1 = RDB[p1+0] - RDB[p0+0];
  a2 = RDB[p1+1] - RDB[p0+1];
  a3 = RDB[p1+2] - RDB[p0+2];
  
  b1 = RDB[p2+0] - RDB[p0+0];
  b2 = RDB[p2+1] - RDB[p0+1];
  b3 = RDB[p2+2] - RDB[p0+2];
  
  c1 = RDB[p3+0] - RDB[p0+0];
  c2 = RDB[p3+1] - RDB[p0+1];
  c3 = RDB[p3+2] - RDB[p0+2];

   /* Calculate vector triple product and volume */
  
  vol = a1*(b2*c3 - b3*c2) - a2*(b1*c3 - b3*c1) + a3*(b1*c2 - b2*c1);
  vol = fabs(vol)/6.0;

  /* Check */
  
  if (vol == 0.0)
    Warn(FUNCTION_NAME, "Zero volume in cell %ld", 
	(long)RDB[cgns + IFC_TET_MSH_IDX]);
  else if (vol > 1E+9)
    Warn(FUNCTION_NAME, "Suspiciouslylarge volume %1.5E in cell %ld", vol, 
	(long)RDB[cgns + IFC_TET_MSH_IDX]);
 
  /* Return volume */

  return vol;
 
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
