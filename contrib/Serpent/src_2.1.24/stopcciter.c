/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : stopcciter.c                                   */
/*                                                                           */
/* Created:       2015/01/19 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Updates coupled calculation stopping criterion               */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StopCCIter:"

/*****************************************************************************/

void StopCCIter()
{
  long niter, tb;

  if(RDB[DATA_SIMULATION_MODE] == (double)SIMULATION_MODE_DYN)
    {
      /* Get time bin index */

      tb = (long)RDB[DATA_DYN_TB];

      /* Get current number of iteration */

      niter = (long)RDB[DATA_SOL_REL_ITER];

      if(tb > 0)
	{
	  if (niter > 4)
	    WDB[DATA_ITERATE] = (double)NO;
	}
      /* Additional iterations for first step?*/
      else if(niter > 0)
	WDB[DATA_ITERATE] = (double)NO;

    }
}
