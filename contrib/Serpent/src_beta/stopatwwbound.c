/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stopatwwbound.c                                */
/*                                                                           */
/* Created:       2015/10/02 (JLe)                                           */
/* Last modified: 2015/10/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Stops track at weight window boundary                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StopAtWWBound:"

/*****************************************************************************/

long StopAtWWBound(long trk, double *x, double *y, double *z, double u, 
		   double v, double w, double spd, double *t, double *lmax,
		   long *cell, long id)
{
  double x0, y0, z0, d;

  /* Check if weight windows are used */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
    return trk;

  /* Move to initial position */

  x0 = *x - *lmax*u;
  y0 = *y - *lmax*v;
  z0 = *z - *lmax*w;

  /* Distance to boundary */

  if ((d = WWDis(x0, y0, z0, u, v, w)) > *lmax)
    return trk;

  /* Extrapolate */

  d = d + EXTRAP_L;

  /* Adjust time and distence */

  *t = *t - (*lmax - d)/spd;
  *lmax = d;

  /* Move over boundary */

  *x = x0 + d*u;
  *y = y0 + d*v;
  *z = z0 + d*w;

  /* Get new location */
	      
  *cell = WhereAmI(*x, *y, *z, u, v, w, id);
  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

  /* Score current */

  ScoreWWDCurr(*x, *y, *z, u, v, w, NO, id);

  /* Exit subroutine */

  return TRACK_END_WWIN;
}

/*****************************************************************************/
