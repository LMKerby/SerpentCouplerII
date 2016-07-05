/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : speed.c                                        */
/*                                                                           */
/* Created:       2012/10/11 (JLe)                                           */
/* Last modified: 2012/10/11 (JLe)                                           */
/* Version:       2.1.9                                                      */
/*                                                                           */
/* Description: Calculates speed from energy                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Speed:"

/*****************************************************************************/

double Speed(long type, double E)
{
  double spd;
  
  /* Avoid compiler warning */

  spd = -1.0;

  /* Check type */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      spd = E/NEUTRON_E0 + 1.0;
      spd = SPD_C*sqrt(1.0 - 1.0/(spd*spd));
    }
  else if (type == PARTICLE_TYPE_GAMMA)
    spd = SPD_C;
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Return speed */

  return spd;
}

/*****************************************************************************/
