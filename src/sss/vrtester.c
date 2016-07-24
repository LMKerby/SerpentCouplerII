#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : vrtester.c                                     */
/*                                                                           */
/* Created:       2011/04/19 (JLe)                                           */
/* Last modified: 2013/05/03 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Tester routine for variance readuction techniques            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "VrTester:"

/*****************************************************************************/

long VrTester(long part, double x, double y, double z, double u, double v,
	      double w, double E, double *wgt, long id)
{
#ifdef mmmmmmmmm
  long new;
  double u0, v0, w0, f;

  /* Check active cycle */
  Die(FUNCTION_NAME, "Ei nyt");
  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return part;

  /* Lasketaan suuntavektori */

  f = sqrt(x*x + y*y);

  u0 = x/f;
  v0 = y/f;
  w0 = z/f;
  
  /* Pistetulo */

  f = u0*u + v0*v;

  /* Check value */

  CheckValue(FUNCTION_NAME, "f", "", f, -1.0, 1.0);

  /* Check */

  if (f > 0.0)
    {
      /* Oikea suunta */

      if (RandF(id) < *wgt*f)
	{
	  /* Adjust weight */

	  *wgt = *wgt/2.0;
      
	  /* Duplicate particle */
	  
	  new = DuplicateParticle(part, id);
	  
	  /* Put variables */
	  
	  WDB[new + PARTICLE_X] = x;
	  WDB[new + PARTICLE_Y] = y;
	  WDB[new + PARTICLE_Z] = z;
	  
	  WDB[new + PARTICLE_U] = u;
	  WDB[new + PARTICLE_V] = v;
	  WDB[new + PARTICLE_W] = w;
	  
	  WDB[new + PARTICLE_E] = E;
	  WDB[new + PARTICLE_WGT] = *wgt;
	  
	  /* Put new particle in que */
	  
	  ToQue(new, id);	  
	}
    }
  else
    {
      /* Väärä suunta */

      if (RandF(id) < -f/(*wgt))
	{
	  /* Russian roulette */

	  if (RandF(id) < 0.5)
	    *wgt = *wgt/0.5;
	  else
	    {
	      /* Put particle back in stack */
	      
	      ToStack(part, id);
	      
	      /* Exit subroutine */
	      
	      return -1;
	    }
	}

    }



#ifdef mmmm

  return part;

  /* Yritetään saada enemmän fotoneja pihviin BBQ-keississä */

  if (RandF(id) < -w/(*wgt))
    {
      /* Russian roulette */

      if (RandF(id) < 0.5)
	*wgt = *wgt/0.5;
      else
	{
	  /* Put particle back in stack */

	  ToStack(part);
	  
	  /* Exit subroutine */
	  
	  return -1;
	}
    }  
  else if (RandF(id) < *wgt*w)
    {
      /* Adjust weight */

      *wgt = *wgt/2.0;
      
      /* Duplicate particle */
      
      new = DuplicateParticle(part);
      
      /* Put variables */
      
      WDB[new + PARTICLE_X] = x;
      WDB[new + PARTICLE_Y] = y;
      WDB[new + PARTICLE_Z] = z;
      
      WDB[new + PARTICLE_U] = u;
      WDB[new + PARTICLE_V] = v;
      WDB[new + PARTICLE_W] = w;
      
      WDB[new + PARTICLE_E] = E;
      WDB[new + PARTICLE_WGT] = *wgt;
      
      /* Put new particle in que */
      
      ToQue(new);	  
    }

#endif

  return part;

#endif 

  return -1.0;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
