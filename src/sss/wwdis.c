#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : wwdis.c                                        */
/*                                                                           */
/* Created:       2015/10/02 (JLe)                                           */
/* Last modified: 2015/10/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Calculates distance to weight window boundary                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WWDis:"

/*****************************************************************************/

double WWDis(double x, double y, double z, double u, double v, double w)
{
  long wwd, msh, fail;
  double min, l;

  /* Check if weight windows are used */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
    Die(FUNCTION_NAME, "Weight windows not in use");

  /* Reset minimum */

  min = INFTY;

  /* Pointer to weight window structure */
  
  wwd = (long)RDB[DATA_PTR_WWD0];
  CheckPointer(FUNCTION_NAME, "(wwd)", DATA_ARRAY, wwd);

  /* Loop over structures */

  while (wwd > VALID_PTR)
    {
      /* Pointer to mesh */
      
      if ((msh = (long)RDB[wwd + WWD_PTR_MESH]) < VALID_PTR)
	{
	  /* Pointer to next */

	  wwd = NextItem(wwd);

	  /* Cycle loop */

	  continue;
	}

      /* Get distance */
      
      if ((l = NearestMeshBoundary(msh, x, y, z, u, v, w, &fail)) > 0.0)
	if (l < min)
	  min = l;

      /* Pointer to next */

      wwd = NextItem(wwd);
    }

  /* Return distance */

  return min;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
