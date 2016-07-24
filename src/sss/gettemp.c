#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gettemp.c                                      */
/*                                                                           */
/* Created:       2012/01/11 (JLe)                                           */
/* Last modified: 2016/02/16 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Returns material temperature at given location               */
/*                                                                           */
/* Comments: - Multifysiikkarajapinnasta haettava lämpötila siirrettiin      */
/*             yhdessä tiheyskertoimen kanssa ifcpoint.c -aliohjelmaan       */
/*             18.3.2013 / 2.1.13.                                           */
/*                                                                           */
/*           - Noi muuttujat pitäis jotenkin pystyä kuitenkin välittämään    */
/*             tonne IFCPoint():lle tässä ja GetTemp():ssä.                  */
/*                                                                           */
/*           - Return value is zero if no temperature is given.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetTemp:"

/*****************************************************************************/

double GetTemp(long mat, long id)
{
  double T, f;

  /***************************************************************************/

  /***** Temperature from interface ******************************************/

  /* Reset values */

  f = 1.0;
  T = -1.0;

  /* Get density factor from multi-physics interface */

  IFCPoint(mat, &f, &T, id);

  /* Check value */

  if (T > 0.0)
    return T;
  
  /***************************************************************************/

  /***** Temperature from material *******************************************/
  
  /* Get Temperature */

  if ((mat > VALID_PTR) && (RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE))
    {
      /* Use minimum temperature */
      
      T = RDB[mat + MATERIAL_TMS_TMIN];

      /* Check with maximum */

      if (T != RDB[mat + MATERIAL_TMS_TMAX])
	Die(FUNCTION_NAME, "Error in temperature (%s) %E %E", 
	    GetText(mat + MATERIAL_PTR_NAME),
	    RDB[mat + MATERIAL_TMS_TMAX], T);
    }
  else
    T = 0.0;

  /* Return value */

  return T;

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
