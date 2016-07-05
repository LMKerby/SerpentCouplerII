/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : matptr.c                                       */
/*                                                                           */
/* Created:       2012/05/10 (JLe)                                           */
/* Last modified: 2015/06/25 (Vva)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Returns material pointer for depletion zones                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MatPtr:"

/*****************************************************************************/

long MatPtr(long mat, long id)
{
  long ptr, lvl0, lvl, div, lst, idx, nax, nrad, nseg, reg, uni, i, j, k, fpe;
  double x, y, z, t, r, phi;

  /* Check null pointer */

  if (mat < VALID_PTR)
    return mat;

  /* Check divisor */

  if ((div = (long)RDB[mat + MATERIAL_PTR_DIV]) < VALID_PTR)
    return mat;

  /* Check type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_NEW)
    Die(FUNCTION_NAME, "Error");

  /***************************************************************************/

  /***** Division into zones *************************************************/

  /* Pointer to material list */

  lst = (long)RDB[div + DIV_PTR_MAT_LIST];
  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);

  /* Check division */

  if ((long)RDB[div + DIV_SEP] == YES)
    {
      /* Get level pointer */

      lvl0 = (long)RDB[DATA_PTR_LVL0];
      CheckPointer(FUNCTION_NAME, "(lvl0)", DATA_ARRAY, lvl0);
      lvl0 = LastItem(lvl0);

      /* Tää on nyt aika kamala viritys. Ongelma on se, että tuo vanha,      */
      /* versioon 2.1.14 implementoitu jako ei toimi jos jaettava materiaali */
      /* ei olekaan ylimmällä levelillä (ylempänä pointteri viimeiseen ja    */
      /* alempana siitä luuppi oikealle levelille.) Tämä on nyt korjattu     */
      /* lisäämällä tohon silmukkaan ylimääräisten levelien lukumäärä. Ei    */
      /* ole takuuta siitä toimiiko kaikissa tilanteissa, mutta ei ainakaan  */
      /* pitäisi rikkoa vanhaa menetelmää joka toimi testitapauksissa.       */
      /* (JLe 10.12.2013 / 2.1.16) */

      i = (long)RDB[DATA_GEOM_LEVELS] - (long)RDB[div + DIV_LVL_MAX] - 1;

     /* Loop to correct level */

      for (j = 0; j < (long)RDB[div + DIV_SEP_LVL] + i; j++)
	if ((lvl0 = PrevItem(lvl0)) < VALID_PTR)
	  Die(FUNCTION_NAME, "Error in level");
      
      /* Pointer to private data */

      lvl = (long)RDB[lvl0 + LVL_PTR_PRIVATE_DATA];
      CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

      /* Get index (tää luetaan nyt level-kohtaisesta muuttujasta,    */
      /* aikaisemmin käytettiin tota globaalia, joka vastaa viimeisen */
      /* tason arvoa (JLe 5.6.2013 / 2.1.15)) */

      idx = (long)GetPrivateData(lvl + LVL_PRIV_ZONE_IDX, id);

      /*
      ptr = (long)RDB[DATA_PTR_ZONE_IDX];
      idx = (long)GetPrivateData(ptr, id);
      */

      /* Find region */

      if ((lst = SeekList(lst, DIV_MAT_LIST_ZONE_IDX, idx, SORT_MODE_ASCEND))
	  < VALID_PTR)
	Die(FUNCTION_NAME, "Unable to find divide zone");
	
      /* Check index */
      
      if (idx != (long)RDB[lst + DIV_MAT_LIST_ZONE_IDX])
	Die(FUNCTION_NAME, "WTF?");
    }

  /***************************************************************************/

  /***** Division into sub-regions *******************************************/

  /* Get universe pointer */

  uni = (long)RDB[lst + DIV_MAT_LIST_PTR_UNIV];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Get pointer to sub-regions */

  reg = (long)RDB[lst + DIV_MAT_LIST_PTR_REG];
  CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

  /* Get bin sizes */

  nax = (long)RDB[div + DIV_NAX];
  nrad = (long)RDB[div + DIV_NRAD];
  nseg = (long)RDB[div + DIV_NSEG];

  /* Check for single-valued */

  if (nax*nrad*nseg == 1)
    idx = 0;
  else
    {
      /* Get local coordinates */

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_X];
      x = GetPrivateData(ptr, id);
      
      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Y];
      y = GetPrivateData(ptr, id);
      
      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Z];
      z = GetPrivateData(ptr, id);		  

      /* Get time */

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_T];
      t = GetPrivateData(ptr, id);

      /* Coordinate change to cold conditions */

      if ((fpe = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) > VALID_PTR)
	CoordExpans(fpe, &x, &y, &z, t, 1);

      /* Axial coordinate */

      if (nax == 1)
	i = 0;
      else
	{
	  /* Calculate normalized coordinate */
	  
	  z = (z - RDB[div + DIV_ZMIN])/
	    (RDB[div + DIV_ZMAX] - RDB[div + DIV_ZMIN]);

	  /* Check that point is inside region */

	  if ((z < 0.0) || (z >= 1.0))
	    Error(div, 
		  "Axial dimension doesn't cover region (geometry error?)");

	  /* Calculate index */

	  i = (long)(z*nax);	  
	}
	    
      /* Radial coordinate */

      if (nrad == 1)
	j = 0;
      else
	{
	  /* Calculate square radius */

          r = x*x + y*y;

          /* Calculate portion of volume inside this radius */

          r = (r - RDB[div + DIV_RMIN]*RDB[div + DIV_RMIN])/
            (RDB[div + DIV_RMAX]*RDB[div + DIV_RMAX] - 
	     RDB[div + DIV_RMIN]*RDB[div + DIV_RMIN]);

          /* Check that point is inside region */
	  
          if ((r < 0.0) || (r >= 1.0))
            Error(div, "Radial dimension doesn't cover region");
	  
          /* Calculate index */

          j = (long)(r*nrad);
	}

      /* Angular segment */

      if (nseg == 1)
	k = 0;
      else
	{
	  /* Calculate angle */

	  phi = PolarAngle(x, y);

	  /* Add tilt */

	  phi = phi + RDB[div + DIV_SEG0];

	  /* Normalize */

	  phi = phi/(2.0*PI);
	  phi = phi - (double)((long)phi);

	  /* Check phi */

	  CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 1.0);

	  /* Calculate index */

	  k = (long)(phi*nseg);
	}
      
      /* Calculate index */
      
      idx = i + j*nax + k*nax*nrad;
      CheckValue(FUNCTION_NAME, "idx", "", idx, 0, nax*nrad*nseg - 1);
    }
  
  /* Get pointer */
  
  mat = (long)RDB[reg + idx];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /***************************************************************************/

  /* Return pointer */

  return mat;
}

/*****************************************************************************/
