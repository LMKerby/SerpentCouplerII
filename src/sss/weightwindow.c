#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : weightwindow.c                                 */
/*                                                                           */
/* Created:       2012/04/13 (JLe)                                           */
/* Last modified: 2015/10/12 (JLe)                                           */
/* Version:       2.1.25                                                     */
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

long WeightWindow(long trk, long part, double x, double y, double z, double u,
		  double v, double w, double E, double *wgt, double t, 
		  long bound, long id)
{
  long np, new;
  double min, max, p;

  /* Check if weight windows are used */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == NO)
    return trk;
  
  /* Get importance */

  if ((p = WWImportance(x, y, z, E)) == 0.0)
    return trk;

  /* Check mode */

  if (1 == 2)
    {
      /***********************************************************************/

      /***** Splitting *******************************************************/
      
      /* Check if called at boundary */

      if (bound == NO)
	return trk;

      /* Get integer part */
      
      np = (long)p;

      /* Sample extra particle */

      if (RandF(id) < p - (double)np)
	np++;

      /* Check */

      if (np < 2)
	return trk;

      /* Divide weight */

      *wgt = *wgt/((double)np);	  

      /* Weight cut-off */

      if (*wgt < 1E-24)
	{
	  /* Put particle back in stack */

	  if (part > VALID_PTR)
	    ToStack(part, id);
	  
	  /* Return weight cut-off */
	      
	  return TRACK_END_WCUT;
	}	
            
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
      
      /* Set multiplicity */

      WDB[new + PARTICLE_MULTIPLICITY] = (double)(np - 2);
	  
      /* Put particle in que */
      
      ToQue(new, id);

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/
      
      /***** Conventional weight windows *************************************/

      /* Boundaries */
      
      min = RDB[DATA_WWD_LOWER_BOUND]/p;
      max = RDB[DATA_WWD_UPPER_BOUND]/p;

      /* Compare weight to boundaries */
      
      if (*wgt < min)
	{
	  /*******************************************************************/

	  /***** Play russian roulette ***************************************/

	  /* Calculate survival probability */
	  
	  p = *wgt/(1.01*min);
	  
	  /* Russian roulette */
	  
	  if (RandF(id) < p)
	    {
	      /* Increase weight */
	      
	      *wgt = *wgt/p;
	    }
	  else
	    {
	      /* Put particle back in stack */
	      
	      if (part > VALID_PTR)
		ToStack(part, id);
	      
	      /* Return weight cut-off */
	      
	      return TRACK_END_WCUT;
	    }
	  
	  /*******************************************************************/
	}
      else if (*wgt > max)
	{
	  /*******************************************************************/

	  /***** Split history ***********************************************/

	  /* Not split if called from source routine */

	  if (part < VALID_PTR)
	    return -1;
	  
	  /* Calculate number of splits */
	  
	  np = (long)(*wgt/max);
	  
	  /* Divide weight */
	  
	  *wgt = *wgt/((double)(np + 1));
	  
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
	  
	  /* Set multiplicity */
	  
	  WDB[new + PARTICLE_MULTIPLICITY] = (double)(np - 1);
	  
	  /* Put particle in que */
	  
	  ToQue(new, id);
	  
	  /*******************************************************************/
	}

      /* Check */

      if ((*wgt < min) || (*wgt > max))
	Die(FUNCTION_NAME, "Weight is out of bounds : %E %E %E\n", min, 
	    *wgt, max);

      /***********************************************************************/
    }

  /* Exit subroutine */

  return trk;
}

/****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
