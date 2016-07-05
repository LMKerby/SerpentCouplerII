/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : samplenu.c                                     */
/*                                                                           */
/* Created:       2011/03/10 (JLe)                                           */
/* Last modified: 2012/11/18 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Samples number of fission neutrons                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleNu:"

/*****************************************************************************/

long SampleNu(double nubar, long id)
{
  double k;
  long nu;

  /* Get eigenvalue */

  if ((long)RDB[DATA_WIELANDT_MODE] != WIELANDT_MODE_NONE)
    {
      k = RDB[DATA_WIELANDT_KP];
      CheckValue(FUNCTION_NAME, "k", "", k, 0.001, 100.0);
    }
  else
    {
      k = RDB[DATA_CYCLE_KEFF];
      CheckValue(FUNCTION_NAME, "k", "", k, 0.001, 10.0);
    }

  /* Adjust nubar */

  nubar = nubar/k;

  /* Get integer part */

  nu = (long)nubar;

  /* Sample extra neutron */

  if (RandF(id) < nubar - (double)nu)
    nu++;

  /* Return value */

  return nu;
}

/*****************************************************************************/
