/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testsurface.c                                  */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2015/03/05 (JLe)                                           */
/* Version:       2.1.23                                                      */
/*                                                                           */
/* Description: Checks if point is inside a surface or not                   */
/*                                                                           */
/* Comments: - Adopted from Serpent 1.1.14                                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestSurface:"

/*****************************************************************************/

long TestSurface(long surf, double x, double y, double z, long on, long id)
{
  double d, l, r, h, R1, R2, alpha, k, A, B, C, D, E, F, G, H, J, K;
  double r2, x0, y0, z0, x1, y1, z1, x2, y2, z2, La, Lb, Lc, th, ps, ph;
  long ptr, type, np, in, n;

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

  if((type != SURF_INF) && (ptr < 0))
    Die(FUNCTION_NAME, "Surface %s has no parameters",
	GetText(surf + SURFACE_PTR_NAME));

#endif

  /* Get number of parameters */

  np = (long)RDB[surf + SURFACE_N_PARAMS];

  /***************************************************************************/

  switch(type)
    {
      /***********************************************************************/

    case SURF_CYL:
    case SURF_CYLZ:
      {
	/* NOTE: Tähän lisättiin katkaisu z-suunnassa 3.12.2009 (JLE) */

	x = x - RDB[ptr];
	y = y - RDB[ptr + 1];
	r = RDB[ptr + 2];
	
	/* Check if cylinder is cut */

	if (np == 5)
	  {
	    /* Test z-planes */
  
	    if (RDB[ptr + 3] > RDB[ptr + 4])
              Error(surf, "Last two parameters given in incorrect order");
	    else if (z < RDB[ptr + 3])
	      return NO;
	    else if (z > RDB[ptr + 4])
	      return NO;
	  }

	if (x*x + y*y > r*r)
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_CYLX:
      {
	y = y - RDB[ptr];
	z = z - RDB[ptr + 1];
	r = RDB[ptr + 2];
	
	/* Check if cylinder is cut */

	if (np == 5)
	  {
	    /* Test x-planes */

	    if (RDB[ptr + 3] > RDB[ptr + 4])
              Error(surf, "Last two parameters given in incorrect order");
	    else if (x < RDB[ptr + 3])
	      return NO;
	    else if (x > RDB[ptr + 4])
	      return NO;
	  }

	if (y*y + z*z > r*r)
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_CYLY:
      {
	x = x - RDB[ptr];
	z = z - RDB[ptr + 1];
	r = RDB[ptr + 2];
	
	/* Check if cylinder is cut */

	if (np == 5)
	  {
	    /* Test y-planes */

	    if (RDB[ptr + 3] > RDB[ptr + 4])
              Error(surf, "Last two parameters given in incorrect order");
	    else if (y < RDB[ptr + 3])
	      return NO;
	    else if (y > RDB[ptr + 4])
	      return NO;
	  }

	if (x*x + z*z > r*r)
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_PAD:
      {
	x = x - RDB[ptr];
	y = y - RDB[ptr + 1];
	R1 = RDB[ptr + 2];
	R2 = RDB[ptr + 3];

	/* Test annular regions */

	r2 = x*x + y*y;

	if (r2 < R1*R1)
	  return NO;
	else if (r2 > R2*R2)
	  return NO;
	else
	  {
	    /* Test first sector surface */

	    alpha = RDB[ptr + 4] - (long)(RDB[ptr + 4]/360)*360.0;
	    
	    if ((alpha == 90.0) || (alpha == 270.0))
	      k = -INFTY;
	    else if ((alpha == -90.0) || (alpha == -270.0))
	      k = INFTY;
	    else
	      k = -sin(PI*alpha/180.0)/cos(PI*alpha/180.0);
	    
	    if (((fabs(alpha) > 90) && (fabs(alpha) <= 270)))
	      {
		if (y < k*x)
		  return NO;
	      }
	    else 
	      {
		if (y > k*x)
		  return NO;
	      }

	    /* Test second sector surface */

	    alpha = RDB[ptr + 5] - (long)(RDB[ptr + 5]/360)*360.0;
	    if ((alpha == 90.0) || (alpha == 270.0))
	      k = -INFTY;
	    else if ((alpha == -90.0) || (alpha == -270.0))
	      k = INFTY;
	    else
	      k = -sin(PI*alpha/180.0)/cos(PI*alpha/180.0);
	    
	    if (((fabs(alpha) > 90) && (fabs(alpha) <= 270)))
	      {
		if (y > k*x)
		  return NO;
	      }
	    else 
	      {
		if (y < k*x)
		  return NO;
	      }
	  }

	return YES;
      }

      /**********************************************************************/

    case SURF_SPH:
      {
	x = x - RDB[ptr];
	y = y - RDB[ptr + 1];
	z = z - RDB[ptr + 2];
	r = RDB[ptr + 3];

	if (x*x + y*y + z*z > r*r)
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_CUBE:
      {
	x = x - RDB[ptr];
	y = y - RDB[ptr + 1];
	z = z - RDB[ptr + 2];
	r = RDB[ptr + 3];
	
        if ((x > r) || (x < -r) || (y > r) || (y < -r)|| (z > r) || (z < -r))
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_CUBOID:
      {
	if (x < RDB[ptr++])
	  return NO;
	else if (x > RDB[ptr++])
	  return NO;
	else if (y < RDB[ptr++])
	  return NO;
	else if (y > RDB[ptr++])
	  return NO;
	else if (z < RDB[ptr++])
	  return NO;
	else if (z > RDB[ptr])
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_RECT:
      {
	if (x < RDB[ptr++])
	  return NO;
	else if (x > RDB[ptr++])
	  return NO;
	else if (y < RDB[ptr++])
	  return NO;
	else if (y > RDB[ptr++])
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_PX:
      {
	x = x - RDB[ptr];

	if (x > 0)
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_PY:
      {
	y = y - RDB[ptr];

	if (y > 0)
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_PZ:
      {
	z = z - RDB[ptr];

	if (z > 0)
	  return NO;
	else
	  return YES;
      }      

      /**********************************************************************/

    case SURF_CONE:
      {
	x = x - RDB[ptr];
	y = y - RDB[ptr + 1];
	z = z - RDB[ptr + 2];
	h = RDB[ptr + 4];

	/* Check half-plane */

	if (h > 0.0)
	  {
	    /* Downward pointing cone */
	    
	    if (z > h)
	      return NO;	    
	  }
	else
	  {
	    /* Upward pointing cone */
	    
	    if (z < h)
	      return NO;	    
	  }

	/* Calculate radius */

	r = RDB[ptr + 3]*(1.0 - z/h);
	
	/* Neutron inside half-plane */

	if (x*x + y*y > r*r)
	  return NO;
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_INF:
      {
	return YES;
      }

      /**********************************************************************/

    case SURF_SQC:
      {
	/* Check definition of rounded corners */

	if (np == 3)
	  {
	    /* Sharp corners */
	    
	    x = x - RDB[ptr];
	    y = y - RDB[ptr + 1];
	    r = RDB[ptr + 2];
	    
	    if (fabs(x) > r)
	      return NO;
	    else if (fabs(y) > r)
	      return NO;
	    else
	      return YES;
	  }
	else
	  {
	    /* Rounded corners. Try outer boundaries first */

	    x = x - RDB[ptr];
	    y = y - RDB[ptr + 1];
	    r = RDB[ptr + 2];
	    
	    if (fabs(x) > r)
	      return NO;
	    else if (fabs(y) > r)
	      return NO;
	    
	    /* Inner boundaries */
	    
	    r = r - RDB[ptr + 3];
	    
	    if (fabs(x) < r)
	      return YES;
	    else if (fabs(y) < r)
	      return YES;
	    
	    /* corners */
	    
	    x = fabs(x) - r;
	    y = fabs(y) - r;
	    
	    if (x*x + y*y < RDB[ptr + 3]*RDB[ptr + 3])
	      return YES;
	    else
	      return NO;
	  }
      }

      /**********************************************************************/

    case SURF_HEXYC:
      {
	/* Check definition of rounded corners */

	if (np == 3)
	  {
	    /* Sharp corners */

	    x = fabs(x - RDB[ptr]);
	    y = fabs(y - RDB[ptr + 1]);
	    r = RDB[ptr + 2];

	    if (y > r)
	      return NO;
	    else if (y > -SQRT3*x + 2*r)
	      return NO;
	    else 
	      return YES;
	  }
	else
	  {
	    /* Rounded corners. Try outer boundaries first */

	    x = fabs(x - RDB[ptr]);
	    y = fabs(y - RDB[ptr + 1]);
	    r = RDB[ptr + 2];
	    
	    if (y > r)
	      return NO;
	    else if (y > -SQRT3*x + 2*r)
	      return NO;
	    
	    /* Inner boundaries */
	    
	    r = r - RDB[ptr + 3]*(1.0 - SIN30);
	    
	    if ((y < r) && (y > SQRT3*x - 2*r))
	      return YES;
	    else if (y < -SQRT3*x + 2*r)
	      return YES;
	    
	    /* Corners */
	    
	    x0 = (RDB[ptr + 2] - RDB[ptr + 3])/COS30;
	    r2 = RDB[ptr + 3]*RDB[ptr + 3];
	    
	    if ((x - x0)*(x - x0) + y*y < r2)
	      return YES;
	    
	    x0 = x0/2.0;
	    y0 = (RDB[ptr + 2] - RDB[ptr + 3]);
	    
	    if ((x - x0)*(x - x0) + (y - y0)*(y - y0) < r2)
	      return YES;
	    
	    /* Neutron is outside */
	    
	    return NO;
	  }
      }

      /**********************************************************************/

    case SURF_HEXXC:
      {
	/* Check definition of rounded corners */

	if (np == 3)
	  {
	    /* Sharp corners */

	    y0 = fabs(x - RDB[ptr]);
	    x = fabs(y - RDB[ptr + 1]);
	    y = y0;
	    r = RDB[ptr + 2];
	    
	    if (y > r)
	      return NO;
	    else if (y > -SQRT3*x + 2*r)
	      return NO;
	    else 
	      return YES;
	  }
	else
	  {
	    /* Rounded corners. Try outer boundaries first */

	    y0 = fabs(x - RDB[ptr]);
	    x = fabs(y - RDB[ptr + 1]);
	    y = y0;
	    r = RDB[ptr + 2];
	    
	    if (y > r)
	      return NO;
	    else if (y > -SQRT3*x + 2*r)
	      return NO;
	    
	    /* Inner boundaries */
	    
	    r = r - RDB[ptr + 3]*(1.0 - SIN30);
	    
	    if ((y < r) && (y > SQRT3*x - 2*r))
	      return YES;
	    else if (y < -SQRT3*x + 2*r)
	      return YES;
	    
	    /* Corners */
	    
	    x0 = (RDB[ptr + 2] - RDB[ptr + 3])/COS30;
	    r2 = RDB[ptr + 3]*RDB[ptr + 3];
	    
	    if ((x - x0)*(x - x0) + y*y < r2)
	      return YES;
	    
	    x0 = x0/2.0;
	    y0 = (RDB[ptr + 2] - RDB[ptr + 3]);
	    
	    if ((x - x0)*(x - x0) + (y - y0)*(y - y0) < r2)
	      return YES;
	    
	    /* Neutron is outside */
	    
	    return NO;
	  }
      }

      /**********************************************************************/

    case SURF_HEXYPRISM:
      {
	/* Check z-boundaries */

	if (RDB[ptr + 3] > RDB[ptr + 4])
	  Error(surf, "Last two parameters given in incorrect order");
	else if (z < RDB[ptr + 3])
	  return NO;
	else if (z > RDB[ptr + 4])
	  return NO;

	/* Transfer co-ordinates */
	
	x = fabs(x - RDB[ptr]);
	y = fabs(y - RDB[ptr + 1]);
	r = RDB[ptr + 2];
	
	/* Test sides */

	if (y > r)
	  return NO;
	else if (y > -SQRT3*x + 2*r)
	  return NO;
	else 
	  return YES;
      }

      /**********************************************************************/
      
    case SURF_HEXXPRISM:
      {
	/* Check z-boundaries */

	if (RDB[ptr + 3] > RDB[ptr + 4])
	  Error(surf, "Last two parameters given in incorrect order");
	else if (z < RDB[ptr + 3])
	  return NO;
	else if (z > RDB[ptr + 4])
	  return NO;

	/* Transfer co-ordinates */

	y0 = fabs(x - RDB[ptr]);
	x = fabs(y - RDB[ptr + 1]);
	y = y0;
	r = RDB[ptr + 2];
	
	/* Test sides */
	
	if (y > r)
	  return NO;
	else if (y > -SQRT3*x + 2*r)
	  return NO;
	else 
	  return YES;
	
      }
      
      /**********************************************************************/
      
    case SURF_CROSS:
      {
	x = fabs(x - RDB[ptr]);
	y = fabs(y - RDB[ptr + 1]);
	l = RDB[ptr + 2];
	d = RDB[ptr + 3];

	/* Check rounded corners (fifth co-ordinate is dummy?) */

	if (np == 5)
	  r = d;
	else
	  r = 0;

	if (x > l)
	  return NO;
	else if (y > l)
	  return NO;
	else if ((x > d) && (y > d))
	  return NO;

	if (r == 0)
	  return YES;

	/* test corners */

	l = l - r;

	if (x > l)
	  {
	    x = x - l;
	    if (x*x + y*y > r*r)
	      return NO;
	    else 
	      return YES;
	  }
	else if (y > l)
	  {
	    y = y - l;
	    if (x*x + y*y > r*r)
	      return NO;
	    else 
	      return YES;
	  }
	else
	  return YES;
      }

      /**********************************************************************/

    case SURF_SVC:
      {

	x = x - RDB[ptr];
	y = y - RDB[ptr + 1];

	d = RDB[ptr + 2];

	if ((fabs(y - x)  < SIN45*(2.96 + 2*d)) &&
	    (fabs(y + x)  < SIN45*(2.96 + 2*d)))
	  return YES;

	if (fabs(x) > 6.87)
	  return NO;

	if (fabs(y) > 6.87)
	  return NO;
       
	if (fabs(x) < d)
	  return YES;

	if (fabs(y) < d)
	  return YES;

	x = x + 4.7;

	if ((fabs(y) < 0.2 + d) &&
	    (fabs(y - x)  < SIN45*(2.5 + 2*d)) &&
	    (fabs(y + x)  < SIN45*(2.5 + 2*d)))
	  return YES;

	x = x - 2*4.7;

	if ((fabs(y) < 0.2 + d) &&
	    (fabs(y - x)  < SIN45*(2.5 + 2*d)) &&
	    (fabs(y + x)  < SIN45*(2.5 + 2*d)))
	  return YES;

	x = x + 4.7;
	y = y + 4.7;

	if ((fabs(x) < 0.2 + d) &&
	    (fabs(y - x)  < SIN45*(2.5 + 2*d)) &&
	    (fabs(y + x)  < SIN45*(2.5 + 2*d)))
	  return YES;

	y = y - 2*4.7;

	if ((fabs(x) < 0.2 + d) &&
	    (fabs(y - x)  < SIN45*(2.5 + 2*d)) &&
	    (fabs(y + x)  < SIN45*(2.5 + 2*d)))
	  return YES;

	return NO;
      }

      /**********************************************************************/

    case SURF_DODE:
      {
	/* Sharp corners */
	
	x = fabs(x - RDB[ptr]);
	y = fabs(y - RDB[ptr + 1]);
	r = RDB[ptr + 2];

	if (np == 4)
	  l = RDB[ptr + 3];
	else
	  l = r;
	if (y > l)
	  return NO;
	else if (x > r)
	  return NO;	
	else if (y > -SQRT3*x + 2*l)
	  return NO;
	else if (x > -SQRT3*y + 2*r)
	  return NO;
	else
	  return YES;
      }
	
      /**********************************************************************/
	  
    case SURF_OCTA:
      {
	/* Sharp corners */
	
	x = fabs(x - RDB[ptr]);
	y = fabs(y - RDB[ptr + 1]);
	r = RDB[ptr + 2];

	if (np == 4)
	  l = RDB[ptr + 3];
	else
	  l = r;

	if (y > r)
	  return NO;
	else if (x > r)
	  return NO;
	else if (y > -x + SQRT2*l)
	  return NO;
	else
	  return YES;
      }
      
      /**********************************************************************/
	  
    case SURF_ASTRA:
      {
	/* Sharp corners */
	
	x = fabs(x - RDB[ptr]);
	y = fabs(y - RDB[ptr + 1]);
	r = RDB[ptr + 2];
	l = RDB[ptr + 3];

	if ((y < r) && (x < r) && (y < -x + SQRT2*l))
	  return YES;
	else 
	  /* Two cylinders */
	  x0 = -0.5*r + 0.75*SQRT2*l;
	  y0 =  0.5*r + 0.25*SQRT2*l;
	  r2 = RDB[ptr + 4]*RDB[ptr + 4];
	  
	  if ((x - x0)*(x - x0) + (y - y0)*(y - y0) < r2)
	    return YES;
	    
	  x0 =  0.5*r + 0.25*SQRT2*l;
	  y0 = -0.5*r + 0.75*SQRT2*l;
	    
	  if ((x - x0)*(x - x0) + (y - y0)*(y - y0) < r2)
	    return YES;
	    
	  /* Neutron is outside */
	    
	  return NO;
      }
    
      /***********************************************************************/	
      
    case SURF_PLANE:
      {
	/* Check number of parameters */

	if (np < 5)
	  {
	    /* Parametric form, reset optional parameters */
	    
	    B = 0.0;
	    C = 0.0;
	    D = 0.0;
	    
	    /* Get parameters */
	    
	    A = RDB[ptr];
	    
	    if (np > 1)
	      B = RDB[ptr + 1];
	    
	    if (np > 2)
	      C = RDB[ptr + 2];
	    
	    if (np > 3)
	      D = RDB[ptr + 3];
	    
	    /* Test */
		
	    if (A*x + B*y + C*z < D)
	      return YES;
	    else
	      return NO;
	  }
	else if (np == 9)
	  {
	    /* Three points, calculate vectors */

	    x0 = -RDB[ptr + 3] + x;
	    y0 = -RDB[ptr + 4] + y;
	    z0 = -RDB[ptr + 5] + z;

	    x1 = -RDB[ptr] + RDB[ptr + 3];
	    y1 = -RDB[ptr + 1] + RDB[ptr + 4];
	    z1 = -RDB[ptr + 2] + RDB[ptr + 5];

	    x2 = -RDB[ptr + 3] + RDB[ptr + 6];
	    y2 = -RDB[ptr + 4] + RDB[ptr + 7];
	    z2 = -RDB[ptr + 5] + RDB[ptr + 8];

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
		  return YES; 
		else
		  return NO;
	      }
	    else
	      {
		/* Check */
		
		if (d < 0.0)
		  return YES; 
		else
		  return NO;
	      }
	  }
	else
	  Die(FUNCTION_NAME, "Invalid number of parameters");
      }

      /**********************************************************************/

    case SURF_QUADRATIC:
      {
	/* Reset optional parameters */

	B = 0.0;
	C = 0.0;
	D = 0.0;
	E = 0.0;
	F = 0.0;
	G = 0.0;
	H = 0.0;
	J = 0.0;
	K = 0.0;

	/* Get parameters */

	A = RDB[ptr];

	if (np > 1)
	  B = RDB[ptr + 1];

	if (np > 2)
	  C = RDB[ptr + 2];

	if (np > 3)
	  D = RDB[ptr + 3];

	if (np > 4)
	  E = RDB[ptr + 4];

	if (np > 5)
	  F = RDB[ptr + 5];

	if (np > 6)
	  G = RDB[ptr + 6];

	if (np > 7)
	  H = RDB[ptr + 7];

	if (np > 8)
	  J = RDB[ptr + 8];

	if (np > 9)
	  K = RDB[ptr + 9];

	/* Test */

	if (A*x*x + B*y*y + C*z*z + D*x*y + E*y*z + F*z*x + G*x + H*y + J*z + K
	    < 0.0)
	  return YES;
	else
	  return NO;
      }
      
      /**********************************************************************/
      
    case SURF_GCROSS:
      {
	x = fabs(x - RDB[ptr]);
	y = fabs(y - RDB[ptr + 1]);

	if ((x > RDB[ptr + 2]) || (y > RDB[ptr + 2]))
	  return NO;

	for (n = 0; n < np - 3; n++)
	  if ((x > RDB[ptr + n + 3]) && (y > RDB[ptr + np - n - 1]))
	    return NO;

	return YES;
      }
      
      /**********************************************************************/
      
    case SURF_PPD:
      {
	/* Coordinate transformation */

	x = x - RDB[ptr];
	y = y - RDB[ptr + 1];
	z = z - RDB[ptr + 2];

	/* Lengths */

	La = RDB[ptr + 3];
	Lb = RDB[ptr + 4];
	Lc = RDB[ptr + 5];

	/* Angles */

	th = RDB[ptr + 6]*PI/180.0;
	ps = RDB[ptr + 7]*PI/180.0;
	ph = RDB[ptr + 8]*PI/180.0;

	/* Test */

	if (x > (La + y*(sin(ps)/cos(ps)) + z*(cos(ph)*sin(th)/cos(th))))
	  return NO;
	else if (x < y*(sin(ps)/cos(ps)) + z*(cos(ph)*sin(th)/cos(th)))
	  return NO;
	else if (y > Lb*cos(ps) + z*(sin(ph)*sin(th)/cos(th)))
	  return NO;
	else if (y < z*(sin(ph)*sin(th)/cos(th)))
	  return NO;
	else if (z > Lc*cos(th))
	  return NO;
	else if (z < 0) 
	  return NO; 
	else 
	  return YES;
      }
      
      /**********************************************************************/

    case SURF_USER:
      {
	/* User-defined surface, call special routine */
	
	UserSurf(2, np, &RDB[ptr], NULL, NULL, NULL, NULL, NULL, NULL, &in,
		 NULL, x, y, z,  0.0, 0.0, 0.0);
	
	return in;
      }
      
      /**********************************************************************/

    default:
      break;
    }
  
  /***************************************************************************/
  
  /***** Invalid surface type ************************************************/
  
  Die(FUNCTION_NAME, "Invalid surface type %ld (surface %s)",
      type, GetText(surf + SURFACE_PTR_NAME));
   
  /* Avoid warnings */

  return 0;
}

/*****************************************************************************/
