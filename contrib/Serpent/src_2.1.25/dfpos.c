/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dfpos.c                                        */
/*                                                                           */
/* Created:       2013/03/06 (JLe)                                           */
/* Last modified: 2013/06/20 (JLe)                                           */
/* Version:       2.1.14                                                     */
/*                                                                           */
/* Description: Calculates surface index and normal vector for discontinuity */
/*              factors                                                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DFPos:"

/* Surfaces */

#define WEST  0
#define SOUTH 1
#define EAST  2
#define NORTH 3

/* Corners */

#define NW 0
#define NE 1
#define SE 2
#define SW 3

/*****************************************************************************/

void DFPos(long surf, double x, double y, double z, long *n1, long *n2,
	   long *m1, long *m2, long *l1, long *l2, double *u0, double *v0, 
	   double *w0)
{
  long ptr, i, j, k, type;
  double div, p;
  const double *param;
  
  /* Get pointer to parameters */

  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get parameters */

  param = &RDB[ptr];

  /* Get type */

  type = (long)RDB[surf + SURFACE_TYPE];

  switch (type)
    {
    case SURF_PX:
      {
	/* Plane perpendicular to x-axis */

	*n1 = 0;
	*n2 = 0;
	*m1 = -1;
	*m2 = -1;
	*l1 = -1;
	*l2 = -1;

	/* Surface normal */

	*u0 = 1.0;
	*v0 = 0.0;
	*w0 = 0.0;

	break;
      }
    case SURF_PY:
      {
	/* Plane perpendicular to y-axis */

	*n1 = 0;
	*n2 = 0;
	*m1 = -1;
	*m2 = -1;
	*l1 = -1;
	*l2 = -1;

	/* Surface normal */

	*u0 = 0.0;
	*v0 = 1.0;
	*w0 = 0.0;

	break;
      }
    case SURF_PZ:
      {
	/* Plane perpendicular to x-axis */

	*n1 = 0;
	*n2 = 0;
	*m1 = -1;
	*m2 = -1;
	*l1 = -1;
	*l2 = -1;

	/* Surface normal */

	*u0 = 0.0;
	*v0 = 0.0;
	*w0 = 1.0;

	break;
      }
    case SURF_PLANE:
      {
	/* General plane perpendicular to x-axis */

	*n1 = 0;
	*n2 = 0;
	*m1 = -1;
	*m2 = -1;
	*l1 = -1;
	*l2 = -1;

	/* Surface normal (ei oo tarkastettu) */

	div = sqrt(param[0]*param[0] + param[1]*param[1] + param[2]*param[2]);

	if (div > 0.0)
	  {
	    *u0 = param[0]/div;
	    *v0 = param[1]/div;
	    *w0 = param[2]/div;
	  }

	break;
      }
    case SURF_SQC:
      {
	/* Square prism */
	
	x = x - param[0];
	y = y - param[1];

	/* Get pitch */

	p = 2.0*param[2] - 5.0*EXTRAP_L;

	/* Get lattice indexes for surfaces */

	GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_SQC);
 
	/* Check indexes */

	if ((i == 0) && (j == 0))
	  {
	    /* Error */

	    Die(FUNCTION_NAME, "Error in indexes");
	  }
	else if ((i > 0) && (j == 0))
	  {
	    /* East face (OK) */

	    *n1 = EAST;
	    *n2 = WEST;

	    /* Surface normal */

	    *u0 = -1.0;
	    *v0 =  0.0;
	    *w0 =  0.0;
	  }
	else if ((i < 0) && (j == 0))
	  {
	    /* West face (OK) */
	    
	    *n1 = WEST;
	    *n2 = EAST;

	    /* Surface normal */

	    *u0 =  1.0;
	    *v0 =  0.0;
	    *w0 =  0.0;
	  }
	else if ((j > 0) && (i == 0))
	  {
	    /* North face (OK) */

	    *n1 = NORTH;
	    *n2 = SOUTH;

	    /* Surface normal */
	    
	    *u0 =  0.0;
	    *v0 = -1.0;
	    *w0 =  0.0;
	  }
	else if ((j < 0) && (i == 0))
	  {
	    /* South face (OK) */
	    
	    *n1 = SOUTH;
	    *n2 = NORTH;

	    /* Surface normal */
	    
	    *u0 =  0.0;
	    *v0 =  1.0;
	    *w0 =  0.0;
	  }
	else
	  {
	    /* This can happen because of numerics */

	    *n1 = -1;
	    *n2 = -1;
	    *m1 = -1;
	    *m2 = -1;
	    *l1 = -1;
	    *l2 = -1;

	    return;
	  }

	/* Get pitch adjusted to mid-point */

	p = 2.0*param[2]*ADF_MID_WIDTH;

	/* Get lattice indexes for surfaces */

	GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_SQC);
 
	/* Check indexes */

	if ((i == 0) && (j == 0))
	  {
	    /* Error */

	    Die(FUNCTION_NAME, "Error in indexes");
	  }
	else if ((i > 0) && (j == 0))
	  {
	    /* East face (OK) */

	    *l1 = EAST;
	    *l2 = WEST;
	  }
	else if ((i < 0) && (j == 0))
	  {
	    /* West face (OK) */
	    
	    *l1 = WEST;
	    *l2 = EAST;
	  }
	else if ((j > 0) && (i == 0))
	  {
	    /* North face (OK) */
	    
	    *l1 = NORTH;
	    *l2 = SOUTH;
	  }
	else if ((j < 0) && (i == 0))
	  {
	    /* South face (OK) */
	    
	    *l1 = SOUTH;
	    *l2 = NORTH;
	  }
	else
	  {
	    /* Not in mid-point region */

	    *l1 = -1;
	    *l2 = -1;
	  }

	/* Get pitch adjusted to corner widths */

	p = 2.0*param[2] - 5.0*EXTRAP_L;	
	p = p*(1.0 - ADF_CORN_WIDTH);

	/* Get lattice indexes for surfaces */

	GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_SQC);
 
	/* Check indexes */

	if ((i > 0) && (j > 0))
	  {
	    /* NE */
	    
	    *m1 = NE;

	    /* Get opposite corner */

	    if (*n1 == NORTH)
	      *m2 = SE;
	    else if (*n1 == EAST)
	      *m2 = NW;
	    else
	      Die(FUNCTION_NAME, "WTF1 %ld", *n1);
	  }

	else if ((i > 0) && (j < 0))
	  {
	    /* SE */
	    
	    *m1 = SE; 

	    /* Get opposite corner */

	    if (*n1 == SOUTH)
	      *m2 = NE;
	    else if (*n1 == EAST)
	      *m2 = SW;
	    else
	      Die(FUNCTION_NAME, "WTF2 %ld", *n1);
	  }
	else if ((i < 0) && (j < 0))
	  {
	    /* SW */

	    *m1 = SW; 
	    
	    /* Get opposite corner */

	    if (*n1 == SOUTH)
	      *m2 = NW;
	    else if (*n1 == WEST)
	      *m2 = SE;
	    else
	      Die(FUNCTION_NAME, "WTF3 %ld", *n1);
	  }
	else if ((i < 0) && (j > 0))
	  {
	    /* NW */

	    *m1 = NW; 

	    /* Get opposite corner */

	    if (*n1 == NORTH)
	      *m2 = SW;
	    else if (*n1 == WEST)
	      *m2 = NE;
	    else
	      Die(FUNCTION_NAME, "WTF4 %ld", *n1);
	  }
	else
	  {
	    *m1 = -1;
	    *m2 = -1;
	  }
	break;
      }
    case SURF_HEXYC:
      {
	/* Y-type hexagonal prism */
	
	x = x - param[0];
	y = y - param[1];

	/* Get pitch */

	p = 2.0*param[2] - 5.0*EXTRAP_L;

	/* Get lattice indexes for surfaces */

	GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_HEXYC);

	/* Check indexes */
	
	if ((i == 1) && (j == 0))
	  {
	    /* North */

	    *n1 = 0;
	    *n2 = 3;

	    /* Surface normal */

	    *u0 =  0.0;
	    *v0 = -1.0;
	    *w0 =  0.0;
	  }
	else if ((i == 0) && (j == 1))
	  {
	    /* North-East */

	    *n1 = 1;
	    *n2 = 4;

	    /* Surface normal */

	    *u0 = -COS30;
	    *v0 = -SIN30;
	    *w0 =  0.0;
	  }
	else if ((i == -1) && (j == 1))
	  {
	    /* South-East */

	    *n1 = 2;
	    *n2 = 5;

	    /* Surface normal */

	    *u0 = -COS30;
	    *v0 = SIN30;
	    *w0 =  0.0;
	  }
	else if ((i == -1) && (j == 0))
	  {
	    /* South */

	    *n1 = 3;
	    *n2 = 0;

	    /* Surface normal */

	    *u0 =  0.0;
	    *v0 =  1.0;
	    *w0 =  0.0;
	  }
	else if ((i ==  0) && (j == -1))
	  {
	    /* South-West */

	    *n1 = 4;
	    *n2 = 1;

	    /* Surface normal */

	    *u0 =  COS30;
	    *v0 =  SIN30;
	    *w0 =  0.0;
	  }
	else if ((i == 1) && (j == -1))
	  {
	    /* North-West */

	    *n1 = 5;
	    *n2 = 2;

	    /* Surface normal */

	    *u0 = COS30;
	    *v0 = -SIN30;
	    *w0 =  0.0;
	  }
	else
	  {
	    /* Error */

	    Die(FUNCTION_NAME, "Error in indexes");
	  }

	/* Mid-points are not calculated */

	*l1 = -1;
	*l2 = -1;

	/* Get pitch adjusted to corner widths */

	p = 2.0*param[2]*(2.0 - ADF_CORN_WIDTH/2.0)/SQRT3;
	
	/* Get lattice indexes for corners */

	GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_HEXXC);

	/* Check indexes */
	
	if ((i == 1) && (j == 0))
	  {
	    *m1 = 0;
	    *m2 = 3;
	  }
	else if ((i == 0) && (j == 1))
	  {
	    *m1 = 1;
	    *m2 = 4;
	  }
	else if ((i == -1) && (j == 1))
	  {
	    *m1 = 2;
	    *m2 = 5;
	  }
	else if ((i == -1) && (j == 0))
	  {
	    *m1 = 3;
	    *m2 = 0;
	  }
	else if ((i ==  0) && (j == -1))
	  {
	    *m1 = 4;
	    *m2 = 1;
	  }
	else if ((i == 1) && (j == -1))
	  {
	    *m1 = 5;
	    *m2 = 2;
	  }
	else
	  {
	    *m1 = -1;
	    *m2 = -1;
	  }

	break;
      }
    case SURF_HEXXC:
      {
	/* X-type hexagonal prism (Tää on otettu suoraan tosta edellisestä  */
	/* vaihtamalla hilatyypit ja suuntavektorit päittäin keskenään. Noi */
	/* kommentit suunnista ei päde.) */
	
	x = x - param[0];
	y = y - param[1];

	/* Get pitch */

	p = 2.0*param[2] - 5.0*EXTRAP_L;

	/* Get lattice indexes for surfaces */

	GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_HEXXC);

	/* Check indexes */
	
	if ((i == 1) && (j == 0))
	  {
	    /* North */

	    *n1 = 0;
	    *n2 = 3;

	    /* Surface normal */

	    *v0 =  0.0;
	    *u0 = -1.0;
	    *w0 =  0.0;
	  }
	else if ((i == 0) && (j == 1))
	  {
	    /* North-East */

	    *n1 = 1;
	    *n2 = 4;

	    /* Surface normal */

	    *v0 = -COS30;
	    *u0 = -SIN30;
	    *w0 =  0.0;
	  }
	else if ((i == -1) && (j == 1))
	  {
	    /* South-East */

	    *n1 = 2;
	    *n2 = 5;

	    /* Surface normal */

	    *v0 = -COS30;
	    *u0 = SIN30;
	    *w0 =  0.0;
	  }
	else if ((i == -1) && (j == 0))
	  {
	    /* South */

	    *n1 = 3;
	    *n2 = 0;

	    /* Surface normal */

	    *v0 =  0.0;
	    *u0 =  1.0;
	    *w0 =  0.0;
	  }
	else if ((i ==  0) && (j == -1))
	  {
	    /* South-West */

	    *n1 = 4;
	    *n2 = 1;

	    /* Surface normal */

	    *v0 =  COS30;
	    *u0 =  SIN30;
	    *w0 =  0.0;
	  }
	else if ((i == 1) && (j == -1))
	  {
	    /* North-West */

	    *n1 = 5;
	    *n2 = 2;

	    /* Surface normal */

	    *v0 = COS30;
	    *u0 = -SIN30;
	    *w0 =  0.0;
	  }
	else
	  {
	    /* Error */

	    Die(FUNCTION_NAME, "Error in indexes");
	  }

	/* Mid-points are not calculated */

	*l1 = -1;
	*l2 = -1;

	/* Get pitch adjusted to corner widths */

	p = 2.0*param[2]*(2.0 - ADF_CORN_WIDTH/2.0)/SQRT3;
	
	/* Get lattice indexes for corners */

	GetLatticeIndexes(p, p, p, x, y, z, &i, &j, &k, SURF_HEXYC);

	/* Check indexes */
	
	if ((i == 1) && (j == 0))
	  {
	    *m1 = 0;
	    *m2 = 3;
	  }
	else if ((i == 0) && (j == 1))
	  {
	    *m1 = 1;
	    *m2 = 4;
	  }
	else if ((i == -1) && (j == 1))
	  {
	    *m1 = 2;
	    *m2 = 5;
	  }
	else if ((i == -1) && (j == 0))
	  {
	    *m1 = 3;
	    *m2 = 0;
	  }
	else if ((i ==  0) && (j == -1))
	  {
	    *m1 = 4;
	    *m2 = 1;
	  }
	else if ((i == 1) && (j == -1))
	  {
	    *m1 = 5;
	    *m2 = 2;
	  }
	else
	  {
	    *m1 = -1;
	    *m2 = -1;
	  }

	break;
      }
    case SURF_CUBE:
    case SURF_CUBOID:
    case SURF_HEXXPRISM:
    case SURF_HEXYPRISM:
    default:
      {
	Error(0, "Surface %s is wrong type for df calculation",
	      GetText(surf + SURFACE_PTR_NAME));
      }
    }
}

/*****************************************************************************/
