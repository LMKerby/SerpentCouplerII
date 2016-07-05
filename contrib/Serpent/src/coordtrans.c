/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : coordtrans.c                                   */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2015/09/15 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Transfers and rotates coordinates and direction cosines      */
/*                                                                           */
/* Comments: - Converted from rotatecoord.c in Serpent 1.1.14                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CoordTrans:"

/*****************************************************************************/

void CoordTrans(long ptr, double *x, double *y, double *z, double *u, 
		double *v, double *w, long id)
{
  long lvl;
  double x0, y0, z0, u0, v0, w0;

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.0, 1.0);

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Loop over transformations */

  while (ptr > VALID_PTR)
    {
      /* Check level pointer */
  
      if ((lvl = (long)RDB[ptr + TRANS_PTR_LVL]) > VALID_PTR)
	{
	  /* Get coordinates and direction cosines */
	  
	  *x = GetPrivateData(lvl + LVL_PRIV_X, id);
	  *y = GetPrivateData(lvl + LVL_PRIV_Y, id);
	  *z = GetPrivateData(lvl + LVL_PRIV_Z, id);
	  *u = GetPrivateData(lvl + LVL_PRIV_U, id);
	  *v = GetPrivateData(lvl + LVL_PRIV_V, id);
	  *w = GetPrivateData(lvl + LVL_PRIV_W, id);
	}
      
      /* Check rotations */
      
      if ((long)RDB[ptr + TRANS_ROT] == YES)
	{
	  /* Copy values to temporary variables */
	  
	  x0 = *x;
	  y0 = *y;
	  z0 = *z;
	  u0 = *u;
	  v0 = *v;
	  w0 = *w;
	  
	  /* Do matrix multiplication */
	  
	  *x = x0*RDB[ptr + TRANS_RX1] + y0*RDB[ptr + TRANS_RX2] 
	    + z0*RDB[ptr + TRANS_RX3];
	  *y = x0*RDB[ptr + TRANS_RX4] + y0*RDB[ptr + TRANS_RX5] 
	    + z0*RDB[ptr + TRANS_RX6];
	  *z = x0*RDB[ptr + TRANS_RX7] + y0*RDB[ptr + TRANS_RX8] 
	    + z0*RDB[ptr + TRANS_RX9];
	  *u = u0*RDB[ptr + TRANS_RX1] + v0*RDB[ptr + TRANS_RX2] 
	    + w0*RDB[ptr + TRANS_RX3];
	  *v = u0*RDB[ptr + TRANS_RX4] + v0*RDB[ptr + TRANS_RX5] 
	    + w0*RDB[ptr + TRANS_RX6];
	  *w = u0*RDB[ptr + TRANS_RX7] + v0*RDB[ptr + TRANS_RX8] 
	    + w0*RDB[ptr + TRANS_RX9];
	  
	  /* Copy values to temporary variables */
	  
	  x0 = *x;
	  y0 = *y;
	  z0 = *z;
	  u0 = *u;
	  v0 = *v;
	  w0 = *w;
	  
	  /* Do matrix multiplication */
	  
	  *x = x0*RDB[ptr + TRANS_RY1] + y0*RDB[ptr + TRANS_RY2] 
	    + z0*RDB[ptr + TRANS_RY3];
	  *y = x0*RDB[ptr + TRANS_RY4] + y0*RDB[ptr + TRANS_RY5] 
	    + z0*RDB[ptr + TRANS_RY6];
	  *z = x0*RDB[ptr + TRANS_RY7] + y0*RDB[ptr + TRANS_RY8] 
	    + z0*RDB[ptr + TRANS_RY9];
	  *u = u0*RDB[ptr + TRANS_RY1] + v0*RDB[ptr + TRANS_RY2] 
	    + w0*RDB[ptr + TRANS_RY3];
	  *v = u0*RDB[ptr + TRANS_RY4] + v0*RDB[ptr + TRANS_RY5] 
	    + w0*RDB[ptr + TRANS_RY6];
	  *w = u0*RDB[ptr + TRANS_RY7] + v0*RDB[ptr + TRANS_RY8] 
	    + w0*RDB[ptr + TRANS_RY9];
	  
	  /* Copy values to temporary variables */
	  
	  x0 = *x;
	  y0 = *y;
	  z0 = *z;
	  u0 = *u;
	  v0 = *v;
	  w0 = *w;
	  
	  /* Do matrix multiplication */
	  
	  *x = x0*RDB[ptr + TRANS_RZ1] + y0*RDB[ptr + TRANS_RZ2] 
	    + z0*RDB[ptr + TRANS_RZ3];
	  *y = x0*RDB[ptr + TRANS_RZ4] + y0*RDB[ptr + TRANS_RZ5] 
	    + z0*RDB[ptr + TRANS_RZ6];
	  *z = x0*RDB[ptr + TRANS_RZ7] + y0*RDB[ptr + TRANS_RZ8] 
	    + z0*RDB[ptr + TRANS_RZ9];
	  *u = u0*RDB[ptr + TRANS_RZ1] + v0*RDB[ptr + TRANS_RZ2] 
	    + w0*RDB[ptr + TRANS_RZ3];
	  *v = u0*RDB[ptr + TRANS_RZ4] + v0*RDB[ptr + TRANS_RZ5] 
	    + w0*RDB[ptr + TRANS_RZ6];
	  *w = u0*RDB[ptr + TRANS_RZ7] + v0*RDB[ptr + TRANS_RZ8] 
	    + w0*RDB[ptr + TRANS_RZ9];
	}

      /* Transfer coordinates */
      
      *x = *x - RDB[ptr + TRANS_X0];
      *y = *y - RDB[ptr + TRANS_Y0];
      *z = *z - RDB[ptr + TRANS_Z0];

      /* Next */

      ptr = NextItem(ptr);
    }

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "", *u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", *v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", *w, -1.0, 1.0);
}

/*****************************************************************************/
