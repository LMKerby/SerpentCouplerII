#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resetoption.c                                  */
/*                                                                           */
/* Created:       2010/11/22 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Resets option (flag) in parameter                            */
/*                                                                           */
/* Comments: - From Serpent 1.1.12                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResetOption:"

/*****************************************************************************/

void ResetOption(long ptr, long opt)
{
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Set option */

  WDB[ptr] = (double)((long)RDB[ptr] & ~opt);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
