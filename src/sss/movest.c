#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : movest.c                                       */
/*                                                                           */
/* Created:       2012/10/05 (JLe)                                           */
/* Last modified: 2015/09/13 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Moves particle forward using surface-tracking                */
/*                                                                           */
/* Comments: - WhereAmI() pitää nyt kutsua juuri ennen NearestBoundary():a,  */
/*             koska jälkimmäinen funktio tarvitsee tallennetut suunta-      */
/*             kosinit, jotka on saattaneet muuttua törmäyksessä. Tässä on   */
/*             potentiaalinen ongelmakohta, jos tuota dataa aletaan käyttää  */
/*             myös muuhun tarkoitukseen, ja se unohdetaan päivittää samalla */
/*             tavalla. --> voisi miettiä jonkun paremman ratkaisun. Voi     */
/*             myös olla että noita kutsuja WhereAmI():hin tehdään           */
/*             tarpeettoman monta.                                           */
/*                                                                           */
/*           - Toi uusi moodi ei siirrä pistettä pinnan yli vaan lähelle     */
/*             sitä, minkä jälkeen seuraava piste arvotaan DT:llä. Tuon      */
/*             pitäisi toimia paremmin STL-geometriatyypin kanssa, ja        */
/*             nopeuttaa erityisesti fuusiogeometrioissa missä on isoja      */
/*             vakuumialueita joihin DT jää helposti jumiin.                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveST:"

/*****************************************************************************/

long MoveST(long part, double totxs, double minxs, long *cell, double *xs0, 
	    double *x, double *y, double *z, double *l0, double u, double v, 
	    double w, long id)
{
  long ptr, type, cell0;
  double d, xs, l, x0, y0, z0;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
 
  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", *x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", *y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", *z, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "u", "", u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", w, -1.0, 1.0);

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];  

  /* Add to ST fraction counter */
		
  ptr = (long)RDB[RES_ST_TRACK_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, 1.0, ptr, id, 2 - type);
  
  /* Get cross section */

  if (totxs < minxs)
    xs = minxs;
  else
    xs = totxs;
  
  /* Check cross section and sample path length */
  
  CheckValue(FUNCTION_NAME, "xs1", "", xs, ZERO, INFTY);
  l = -log(RandF(id))/xs;

  /* This call must be made to get the stored direction cosines right */
  /* (may have changed in a collision or when symmetries are invoked) */
  
  cell0 = WhereAmI(*x, *y, *z, u, v, w, id);
  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell0);
  
  /* Get distance to nearest boundary */

  d = NearestBoundary(id);

  x0 = *x;
  y0 = *y;
  z0 = *z;

  /* Compare distances */

  if (l < d)
    {
      /* Remember initial coordinates */

      x0 = *x;
      y0 = *y;
      z0 = *z;

      /* Move particle to collision site */

      *x = *x + u*l;
      *y = *y + v*l;
      *z = *z + w*l;

      /* Find location */

      *cell = WhereAmI(*x, *y, *z, u, v, w, id);
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);
      
      /* Tätä tarvitaan hoitamaan geometriavirheet UMSH-geometrioissa */
      /* (Mikähän tääkin on? 22.7.2015) */

      if ((long)RDB[DATA_PTR_UMSH0] > VALID_PTR)      
	if ((*cell < VALID_PTR) || 
	    ((long)RDB[*cell + CELL_PTR_MAT] < VALID_PTR))
	  {
	    /* Set cross section and path length */
	    
	    *xs0 = -1.0;
	    *l0 = l;
	    
	    /* Check distance */

	    CheckValue(FUNCTION_NAME, "l0", "1", *l0, ZERO, INFTY);

	    /* Return virtual collision */

	    return TRACK_END_VIRT;
	  }      

      /* Check with previous (tää antoi aiheettomia varoituksia universumi- */
      /* symmetrioiden kanssa ennen kuin symmetriarajapinta lisättiin */
      /* nearestboundary.c:hen). */

      if (cell0 > VALID_PTR)
	if (*cell != cell0)
	  {
	    printf("x = %1.14E; y = %1.14E; z = %1.14E; u = %1.14E; v = %1.14E; w = %1.14E;\n", x0, y0, z0, u, v, w);
	  Warn(FUNCTION_NAME, 
	      "Surface crossing missed between [%1.2E, %1.2E, %1.2E] and [%1.2E %1.2E %1.2E] (d = %E, l = %E)",
	      x0, y0, z0, *x, *y, *z, d, l);
	  }

      /* Check cross section (could be from previous time step) */

      CheckValue(FUNCTION_NAME, "xs2", "", xs, ZERO, INFTY);

      /* Set cross section and path length */

      *xs0 = xs;
      *l0 = l;

      /* Check distance */

      CheckValue(FUNCTION_NAME, "l0", "2", *l0, ZERO, INFTY);

      /* Sample between real and virtual collision */

      if (RandF(id) < totxs/minxs)
	{
	  /* Return collision */

	  return TRACK_END_COLL;
	}
      else
	{
	  /* Return virtual */

	  return TRACK_END_VIRT;
	}
    }
  else
    {
      /* Check mode */

      if (1 != 2)
	{
	  /* Conventional surface-tracking, move particle over surface */

	  *x = *x + u*(d + EXTRAP_L);
	  *y = *y + v*(d + EXTRAP_L);
	  *z = *z + w*(d + EXTRAP_L);
	  
	  /* Find location */
	  
	  *cell = WhereAmI(*x, *y, *z, u, v, w, id);
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);
	  
	  /* Set cross section and path length */
	  
	  *xs0 = -1.0;
	  *l0 = d + EXTRAP_L;
	  
	  /* Check distance */

	  if (*l0 < 0.0)
	    printf("%1.14E %1.14E %1.14E : %1.14E %1.14E %1.14E\n", x0, y0, z0, u, v, w);

	  CheckValue(FUNCTION_NAME, "l0", "3", *l0, ZERO, INFTY);

	  /* Return surface crossing */
	  
	  return TRACK_END_SURF;
	}
      else
	{
	  /* New mode, adjust surface distance */

	  if ((d = d - 0.05) < 0.0)
	    d = 0.0;

	  /* Move particle near surface */

	  *x = *x + u*d;
	  *y = *y + v*d;
	  *z = *z + w*d;
	  
	  /* Find location (tarvitaankohan tätä muuhun kuin */
	  /* tohon testiin?) */
	  
	  *cell = WhereAmI(*x, *y, *z, u, v, w, id);
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

	  /* Check */

	  if (cell0 > VALID_PTR)
	    if (*cell != cell0)
	      Die(FUNCTION_NAME, "Possible overlap at (%1.2E, %1.2E, %1.2E)",
		  *x, *y, *z);

	  /* Set cross section and path length */
	  
	  *xs0 = -1.0;
	  *l0 = d;

	  /* Check distance */

	  CheckValue(FUNCTION_NAME, "l0", "4", *l0, 0.0, INFTY);

	  /* Enforce delta-tracking for next collision */

	  ptr = (long)RDB[DATA_DT_ENFORCE_NEXT_TRACK];
	  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
	  PutPrivateData(ptr, YES, id);
	  
	  /* Return virtual collision */
	  
	  return TRACK_END_VIRT;
	}
    }

  /* Avoid compiler warning */

  return -1;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
