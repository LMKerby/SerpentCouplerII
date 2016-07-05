/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : gettemp.c                                      */
/*                                                                           */
/* Created:       2012/01/11 (JLe)                                           */
/* Last modified: 2015/05/29 (JLe)                                           */
/* Version:       2.1.24                                                     */
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
  long ptr, ncol, tfb, nst, reg, uni;
  double T, f, x, y, r;

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

  /***** Temperature from feedback iteration *********************************/

  /* Check if feedback is in use */

  if ((long)RDB[DATA_USE_TFB] == YES)
    {
      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      ncol = (long)GetPrivateData(ptr, id);
      
      /* Loop over temperature feedbacks */
      
      tfb = (long)RDB[DATA_PTR_TFB0];
      while (tfb > VALID_PTR)
	{
	  /* Pointer to nest */
	  
	  nst = (long)RDB[tfb + TFB_PTR_NST];
	  CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);
	  
	  /* Get pointer to region */
	  
	  if ((reg = (long)TestValuePair(nst + NEST_PTR_COL_REG, ncol, id))
	      > VALID_PTR)
	    {
	      /* Pointer to feedback region */
	      
	      if ((reg = (long)RDB[reg + NEST_REG_PTR_TFB_REG]) > VALID_PTR)
		{
		  /* Pointer to universe */

		  uni = (long)RDB[nst + NEST_PTR_UNI];
		  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
		  
		  /* Get coordinates */
		  
		  ptr = RDB[uni + UNIVERSE_PTR_PRIVA_X];
		  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
		  x = GetPrivateData(ptr, id);

		  ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Y];
		  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
		  y = GetPrivateData(ptr, id);
		  
		  /*
		  ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Z];
		  z = GetPrivateData(ptr, id);		  
		  */

		  /* Get temperature (tähän voi sitten lyödä sen mallin */
		  /* joka laskee pinnin säteestä riippuvan lämpötilan.  */
		  /* Noi koordinaatit on suhteessa pinnin origoon.) */

                  r= sqrt(x*x + y*y);
		  
		  T = RDB[reg + TFB_REG_ITER_C0]*r*r +
		    RDB[reg + TFB_REG_ITER_C1]*log(r) + 
		    RDB[reg + TFB_REG_ITER_C2];
		  
		  /* Check value */

		  if (T < 0.0)
		    Die(FUNCTION_NAME, "Negative temperature in %s",
			GetText(mat + MATERIAL_PTR_NAME));
		  else if (T < RDB[mat + MATERIAL_TMS_TMIN])
		    {
		      /* Score failure */
		      
		      ptr = RDB[RES_TMS_FAIL_STAT];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      AddBuf1D(1.0, 1.0, ptr, id, 2);

		      /* Print warning */
		      /*
		      Warn(FUNCTION_NAME, 
			   "Temperature %E below maximum (%E) (%s) ",
			   T, RDB[mat + MATERIAL_TFB_TMIN],
			   GetText(mat + MATERIAL_PTR_NAME));
		      */
		      /* Set limiting value */

		      T = RDB[mat + MATERIAL_TMS_TMIN];
		    }
		  else if (T > RDB[mat + MATERIAL_TMS_TMAX])
		    {
		      /* Score failure */

		      ptr = RDB[RES_TMS_FAIL_STAT];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      AddBuf1D(1.0, 1.0, ptr, id, 3);
		      
		      /* Print warning */
		      /*
		      Warn(FUNCTION_NAME, 
			   "Temperature %E above maximum (%E) (%s) ",
			   T, RDB[mat + MATERIAL_TFB_TMAX],
			   GetText(mat + MATERIAL_PTR_NAME));
		      */
		      /* Set limiting value */

		      T = WDB[mat + MATERIAL_TMS_TMAX];
		    }

		  /* Return value */

		  return T;
		}
	      else
		{
		  /* Break loop */
		  
		  break;
		}
	    }

	  /* Next feedback */

	  tfb = NextItem(tfb);
	}
    }

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
