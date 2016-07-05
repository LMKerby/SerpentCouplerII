/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : profile.c                                      */
/*                                                                           */
/* Created:       2014/12/24 (JLe)                                           */
/* Last modified: 2014/12/24 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Records CPU time interval for profiling                      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Profile:"

/*****************************************************************************/

#ifdef PROFILING

void Profile(double *t0, long bin, long id)
{
  long ptr;
  double t;
  
  /* Check if timer is initialized */

  if (bin < 0)
    *t0 = clock();
  else
    {
      /* Calculate interval */
  
      t = clock() - *t0;

      /* Store */
      
      ptr = (long)RDB[RES_CPU_TIME_PROFILE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(t, 1.0, ptr, id, bin);
    }
}

#endif

/*****************************************************************************/
