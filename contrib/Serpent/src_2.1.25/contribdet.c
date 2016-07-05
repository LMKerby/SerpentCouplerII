/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : contribdet.c                                   */
/*                                                                           */
/* Created:       2015/10/02 (JLe)                                           */
/* Last modified: 2015/10/02 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Scores detector contributions (for the lack of better term)  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ContribDet:"

/*****************************************************************************/

void ContribDet(long part, long det1, long ptr, long idx0, long rbin0, 
		double val0, double wgt0, double tcross, long id)
{
  long evn, mat1, rbin1, loc1, idx1;
  double wgt1, flx1, x1, y1, z1, E1, t1, g1, f1, val1;

  /* Check pointers */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  CheckPointer(FUNCTION_NAME, "(det1)", DATA_ARRAY, det1);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Loop over events */
  
  evn = (long)RDB[part + PARTICLE_PTR_EVENTS];
  while (evn > VALID_PTR)
    {
      /* Get weight, flux and time */
      
      wgt1 = RDB[evn + EVENT_WGT];
      flx1 = RDB[evn + EVENT_FLX];
      t1 = RDB[evn + EVENT_T];
      
      /* Check */
      
      if ((flx1*wgt1 < ZERO) || ((tcross > 0.0) && (t1 > tcross)))
	{
	  /* Pointer to next */
	  
	  evn = NextItem(evn);
	  
	  /* Cycle loop */
	  
	  continue;
	}
      
      /* Get particle coordinates */
      
      x1 = RDB[evn + EVENT_X];
      y1 = RDB[evn + EVENT_Y];
      z1 = RDB[evn + EVENT_Z];
      E1 = RDB[evn + EVENT_E];
      mat1 = RDB[evn + EVENT_PTR_MAT];

      /* Check flux, weight, energy and time */

      CheckValue(FUNCTION_NAME, "flx1", "", flx1, 0.0, INFTY);
      CheckValue(FUNCTION_NAME, "wgt1", "", wgt1, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "E1", "", E1, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "t1", "", t1, ZERO, INFTY);
      
      /* Reset response index */
      
      rbin1 = 0;
      
      /* Get pointer to response functions */
      
      loc1 = (long)RDB[det1 + DET_PTR_RBINS];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
      
      /* Loop over responses */
      
      while (loc1 > VALID_PTR)
	{
	  /* Get bin index */
	  
	  if ((idx1 = DetBin(det1, mat1, x1, y1, z1, E1, t1, id)) < 0)
	    break;
 
	  /* Get density factor */

	  g1 = DensityFactor(mat1, x1, y1, z1, t1, id);
	  CheckValue(FUNCTION_NAME, "g", "", g1, 0.0, 1.0);
	  
	  /* Get response function */
	  
	  f1 = DetResponse(det1, loc1, mat1, E1, g1, id);
	  
	  /* Calculate value */
	  
	  val1 = f1*flx1;
	  
	  /* Score */
	  
	  AddBuf(val0*val1, wgt0, ptr, id, -1, idx0, idx1, rbin0, rbin1);

	  /* Update response index */
	  
	  rbin1++;
		      
	  /* Next response */
	  
	  loc1 = NextItem(loc1);
	}

      /* Next event */
		  
      evn = NextItem(evn);
    }
}

/*****************************************************************************/
