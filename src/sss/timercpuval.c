#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : timercpuval.c                                  */
/*                                                                           */
/* Created:       2010/12/11 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns CPU time from timer                                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "TimerCPUVal:"

/*****************************************************************************/

double TimerCPUVal(long i)
{
  double t;

  /* Check index */

  if ((i < 0) || (i > TOT_TIMERS))
    Die(FUNCTION_NAME, "Invalid timer index %d", i);

  /* Return total time if timer is stopped */

  if (timer[i].on == NO)
    return timer[i].cpu_t/((double)CLOCKS_PER_SEC);

  /* Get current CPU time */
  
  t = clock();
  
  /* Return time */

  return (timer[i].cpu_t + t - timer[i].cpu_t0)/((double)CLOCKS_PER_SEC);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
