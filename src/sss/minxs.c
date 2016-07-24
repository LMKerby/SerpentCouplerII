#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : minxs.c                                        */
/*                                                                           */
/* Created:       2012/11/04 (JLe)                                           */
/* Last modified: 2014/04/04 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Determines the minimum cross section used for sampling       */
/*              path lengths                                                 */
/*                                                                           */
/* Comments: - The idea is to get better statistics for the collision        */
/*             estimator.                                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MinXS:"

/*****************************************************************************/

double MinXS(long type, double spd, long id)
{
  long ptr;
  double minxs, dl, dt;

  /* Get minimum cross section (NOTE: toi arvo) */

  if ((long)RDB[DATA_OPT_USE_DT] == NO)
    minxs = 0.2;
  else if (type == PARTICLE_TYPE_NEUTRON)
    minxs = RDB[DATA_MIN_NMACROXS];
  else
    minxs = RDB[DATA_MIN_PMACROXS];

  /* Get dinstance and time intervals */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      dl = RDB[DATA_CFE_N_MIN_L];
      dt = RDB[DATA_CFE_N_MIN_T];
    }
  else
    {
      dl = RDB[DATA_CFE_G_MIN_L];
      dt = RDB[DATA_CFE_G_MIN_T];
    }

  /* Compare to minimum */

  if (dl > 0.0)
    if (minxs < 1.0/dl)
      minxs = 1.0/dl;

  if ((dl = dt*spd) > 0.0)
    if (minxs < 1.0/dl)
      minxs = 1.0/dl;

  /* Check value */
  
  CheckValue(FUNCTION_NAME, "minxs", "", minxs, ZERO, 1E+12);

  /* Score */

  ptr = (long)RDB[RES_MIN_MACROXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(minxs, 1.0, ptr, id, 2 - type);

  /* Return value */

  return minxs;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
