/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : weightwindow.c                                 */
/*                                                                           */
/* Created:       2012/04/13 (JLe)                                           */
/* Last modified: 2012/10/10 (JLe)                                           */
/* Version:       2.1.9                                                      */
/*                                                                           */
/* Description: Handles particle hit in weight window                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WeightWindow:"

/*****************************************************************************/

long WeightWindow(long wwp, long part, double x, double y, double z, double u,
		  double v, double w, double E, double *wgt, double t, 
		  long id)
{
  long np, n, kill, new;
  double min, max;
  return NO;
  /* Check pointer */
  /*
  CheckPointer(FUNCTION_NAME, "(wwp)", DATA_ARRAY, wwp);
  */
  /* Get weight boundaries */
  /*
  min = RDB[wwp + WWIN_MIN_WGT];
  max = RDB[wwp + WWIN_MAX_WGT];
  */

  min = 0.1;
  max = 0.5;
  /*
  min = 0.0;
  max = 10.0;

  if ((x > -15.0) && (x < 0.0))
    max = (40 - x)/(40 - -15);
  else
    max = INFTY;
  */
  min = 0.1*max;

  /* Reset kill-flag */

  kill = NO;



  /* Compare weight to boundaries */

  if (*wgt < min)
    {
      /************************************************************************/

      /***** Play russian roulette ********************************************/

      /* Russian roulette */
  
      if (RandF(id) < RDB[DATA_OPT_ROULETTE_P0])
	{
	  /* Preserve particle */

	  kill = NO;

	  /* Increase weight */

	  *wgt = *wgt/RDB[DATA_OPT_ROULETTE_P0];
	}
      else
	{
	  /* Kill particle */

	  kill = YES;

	  /* Put particle back in stack */
	  
	  ToStack(part, id);
	}

      /************************************************************************/
    }
  else if (*wgt > max)
    {
      /************************************************************************/

      /***** Split history ****************************************************/

      /* Calculate number of splits */

      np = (long)(*wgt/max);

      /* Divide weight */

      *wgt = *wgt/((double)(np + 1));

      /* Loop over splits */

      for (n = 0; n < np; n++)
	{
      	  /* Duplicate incident neutron */
	  
	  new = DuplicateParticle(part, id);

	  /* Put state variables */

	  WDB[new + PARTICLE_X] = x;
	  WDB[new + PARTICLE_Y] = y;
	  WDB[new + PARTICLE_Z] = z;

	  WDB[new + PARTICLE_U] = u;
	  WDB[new + PARTICLE_V] = v;
	  WDB[new + PARTICLE_W] = w;
	  
	  WDB[new + PARTICLE_E] = E;
	  WDB[new + PARTICLE_WGT] = *wgt;
	  WDB[new + PARTICLE_T] = t;

	  /* Put particle in que */

	  ToQue(new, id);
	}

      /* Preserve initial particle */

      kill = NO;

      /************************************************************************/
    }

  /* Return status */

  return kill;
}

/*****************************************************************************/
