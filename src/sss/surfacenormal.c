#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : surfacenormal.c                                */
/*                                                                           */
/* Created:       2014/02/07 (JLe)                                           */
/* Last modified: 2014/12/01 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Calculates normal vector for surface at (x,y,z)              */
/*                                                                           */
/* Comments: - Not called from anywhere                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SurfaceNormal:"

/*****************************************************************************/

void SurfaceNormal(long surf, double x, double y, double z, 
		   double *u, double *v, double *w, long id)
{
  double d, A, B, C, x1, y1, z1, x2, y2, z2, x3, y3, z3;
  long ptr, type, np, n;

  Die(FUNCTION_NAME, "Not used");

  /* Check surface pointer */

  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Do coordinate transformation (use dummies for direction cosines) */

  x1 = 1.0;
  y1 = 0.0;
  z1 = 0.0;

  if ((ptr = (long)RDB[surf + SURFACE_PTR_TRANS]) > VALID_PTR)
    CoordTrans(ptr, &x, &y, &z, &x1, &y1, &z1, id);

  /* Get surface type */

  type = (long)RDB[surf + SURFACE_TYPE];

  /* Pointer to parameter list */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  
  /* Check number of parameters */

#ifdef DEBUG

  if((type != SURF_INF) && (ptr < VALID_PTR))
    Die(FUNCTION_NAME, "Surface %s has no parameters",
	GetText(surf + SURFACE_PTR_NAME));
#endif

  /* Get number of parameters */

  np = (long)RDB[surf + SURFACE_N_PARAMS];
  n = 0;

  /***************************************************************************/

  switch(type)
    {
      /**********************************************************************/

    case SURF_PX:
      {
	/* Put vector */

	*u = 1.0;
	*v = 0.0;
	*w = 0.0;

	/* Break case */

	break;
      }

      /**********************************************************************/

    case SURF_PY:
      {
	/* Put vector */

	*u = 0.0;
	*v = 1.0;
	*w = 0.0;

	/* Break case */

	break;
      }

      /**********************************************************************/

    case SURF_PZ:
      {
	/* Put vector */

	*u = 0.0;
	*v = 0.0;
	*w = 1.0;

	/* Break case */

	break;
      }

      /***********************************************************************/	
      
    case SURF_PLANE:
      {
	/* Check number of parameters */

	if (np == 9)
	  {
	    /* Surface defined by three points, read coordinates */
	    
	    x1 = RDB[ptr + n++];
	    y1 = RDB[ptr + n++];
	    z1 = RDB[ptr + n++];
	    x2 = RDB[ptr + n++];
	    y2 = RDB[ptr + n++];
	    z2 = RDB[ptr + n++];
	    x3 = RDB[ptr + n++];
	    y3 = RDB[ptr + n++];
	    z3 = RDB[ptr + n++];
	    
	    /* Calculate coefficients */
	    
	    A = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
	    B = z2*x3 - z3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
	    C = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
	  }
	else
	  {
	    /* Surface defined by parameters, reset optional*/

	    B = 0.0;
	    C = 0.0;

	    /* Get parameters */
	    
	    A = RDB[ptr + n++];
	    
	    if (n < np)
	      B = RDB[ptr + n++];
	    
	    if (n < np)
	      C = RDB[ptr + n++];
	  }

	/* Calculate denominator */

	if ((d = sqrt(A*A + B*B + C*C)) == 0.0)
	  Die(FUNCTION_NAME, "Error in surface definition");
	
	/* Calculate components */
	
	*u = A/d;
	*v = B/d;
	*w = C/d;

	/* Break case */

	break;
      }
      
      /**********************************************************************/

    default:
      Die(FUNCTION_NAME, "Invalid surface type");
    }

  /* Check direction cosines */
  
  CheckValue(FUNCTION_NAME, "Normal vector", "",
	     (*u)*(*u) + (*v)*(*v) + (*w)*(*w) - 1.0, -1E-4, 1E-4);

  /***************************************************************************/

}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
