#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorewwdcurr.c                                 */
/*                                                                           */
/* Created:       2015/10/18 (JLe)                                           */
/* Last modified: 2015/10/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Scores interface currents for weight window importances      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreWWDCurr:"

/*****************************************************************************/

void ScoreWWDCurr(double x1, double y1, double z1, double u, double v, 
		  double w, long src, long id)
{
  long wwd, msh, loc0, loc1, ptr;
  double x0, y0, z0;

  /* Check if weight windows are used */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
    {
      if (src == NO)
	Die(FUNCTION_NAME, "Weight windows not in use");
      else
	return;
    }

  /* Point before crossing */

  x0 = x1 - 2.0*EXTRAP_L*u;
  y0 = y1 - 2.0*EXTRAP_L*v;
  z0 = z1 - 2.0*EXTRAP_L*w;
  
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

      /* Check source mode */

      if (src == YES)
	{
	  /* Score source rate (number of particles) */

	  if ((loc1 = MeshPtr(msh, x1, y1, z1)) > VALID_PTR)
	    {
	      /* Pointer to structure */
	  
	      loc1 = (long)RDB[loc1];
	      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	      
	      /* Score */
	      
	      ptr = (long)RDB[loc1 + WWD_MESH_RES_SRC_RATE];
	      AddBuf1D(1.0, 1.0, ptr, id, 0);
	    }
	  
	  /* Pointer to next */

	  wwd = NextItem(wwd);

	  /* Cycle loop */

	  continue;
	}
      
      /* Score outward current (number of particles) */

      if ((loc0 = MeshPtr(msh, x0, y0, z0)) > VALID_PTR)
	{
	  /* Pointer to structure */
	  
	  loc0 = (long)RDB[loc0];
	  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

	  /* Score */
	  
	  ptr = (long)RDB[loc0 + WWD_MESH_RES_OUT_CURR];
	  AddBuf1D(1.0, 1.0, ptr, id, 0);
	}

      /* Score inward current (number of particles) */

      if ((loc1 = MeshPtr(msh, x1, y1, z1)) > VALID_PTR)
	{
	  /* Pointer to structure */
	  
	  loc1 = (long)RDB[loc1];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  /* Score */
	  
	  ptr = (long)RDB[loc1 + WWD_MESH_RES_IN_CURR];
	  AddBuf1D(1.0, 1.0, ptr, id, 0);
	}

      /* Check */

      if ((loc0 > VALID_PTR) && (loc1 > VALID_PTR) && (loc0 != loc1))
	{
	  /* Loop over currents */
	  
	  ptr = (long)RDB[loc0 + WWD_MESH_PTR_CURR];
	  while (ptr > VALID_PTR)
	    {
	      /* Compare pointers */
	      
	      if ((long)RDB[ptr + WWD_MESH_CURR_PTR_NEIGHBOUR] == loc1)
		break;
	      
	      /* Pointer to next */
	      
	      ptr = NextItem(ptr);
	    }
	  
	  /* Check pointer */
	  
	  if (ptr < VALID_PTR)
	    Warn(FUNCTION_NAME, "Neighbour cell not found (%E %E %E)",
		 x1, y1, z1);
	}
      
      /* Pointer to next */

      wwd = NextItem(wwd);
    }
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
