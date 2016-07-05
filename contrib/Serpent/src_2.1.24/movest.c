/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : movest.c                                       */
/*                                                                           */
/* Created:       2012/10/05 (JLe)                                           */
/* Last modified: 2014/05/30 (JLe)                                           */
/* Version:       2.1.22                                                     */
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
  double d, xs, l;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

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

  /* Compare distances */

  if (l < d)
    {
      /* Move particle to collision site */

      *x = *x + u*l;
      *y = *y + v*l;
      *z = *z + w*l;

      /* Find location */

      *cell = WhereAmI(*x, *y, *z, u, v, w, id);
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);
      
      /* Tätä tarvitaan hoitamaan geometriavirheet UMSH-geometrioissa */

      if ((long)RDB[DATA_PTR_UMSH0] > VALID_PTR)      
	if ((*cell < VALID_PTR) || 
	    ((long)RDB[*cell + CELL_PTR_MAT] < VALID_PTR))
	  {
	    /* Set cross section and path length */
	    
	    *xs0 = -1.0;
	    *l0 = l;
	    
	    /* Return virtual collision */

	    return TRACK_END_VIRT;
	  }      

      /* Check with previous (tää antoi aiheettomia varoituksia universumi- */
      /* symmetrioiden kanssa ennen kuin symmetriarajapinta lisättiin */
      /* nearestboundary.c:hen). */

      if (cell0 > VALID_PTR)
	if (*cell != cell0)
	  Warn(FUNCTION_NAME, "Possible overlap at (%1.2E, %1.2E, %1.2E)",
	       *x, *y, *z);

      /* Tää ei nyt varmaan mee ihan oikein, mutta tehdään näin että saadaan */
      /* M&C 2013 -laskut menemään läpi. */

      /* NOTE: Selvitä mitä toi edellinen kommentti meinaa. Tota pointteria */
      /* ei nyt voi käyttää mihinkään syystä että ton col-flagin takia      */
      /* WhereAmI():a ei kutsuta joka kerta (pitäisi nopeuttaa unstructured */
      /* meshissä träkkäystä). Pointterin voisi sen sijaan ottaa talteen,   */
      /* ja käyttää sitä */

      /* NOTE 28.5.2014 / 2.1.21: tän voi varmaan poistaa */

#ifdef mmmmmmmmmmmmm

      if (*cell != cell0)
	{
	  /* Set cross section and path length */

	  *xs0 = -1.0;
	  *l0 = l;
 
	  /* Return virtual collision */

	  return TRACK_END_VIRT;
	}

#endif

      /* Check cross section (could be from previous time step) */

      CheckValue(FUNCTION_NAME, "xs2", "", xs, ZERO, INFTY);

      /* Set cross section and path length */

      *xs0 = xs;
      *l0 = l;

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
      /* Move particle over the surface */

      *x = *x + u*(d + EXTRAP_L);
      *y = *y + v*(d + EXTRAP_L);
      *z = *z + w*(d + EXTRAP_L);

      /* Find location */

      *cell = WhereAmI(*x, *y, *z, u, v, w, id);
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

      /* Set cross section and path length */

      *xs0 = -1.0;
      *l0 = d + EXTRAP_L;
 
      /* Return surface crossing */

      return TRACK_END_SURF;
    }

  /* Avoid compiler warning */

  return -1;
}

/*****************************************************************************/
