/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : setoption.c                                    */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2011/11/11 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Sets option (flag) in parameter                              */
/*                                                                           */
/* Comments: - From Serpent 1.1.12                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SetOption:"

/*****************************************************************************/

void SetOption(long ptr, long opt)
{
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Set option */

  WDB[ptr] = (double)((long)RDB[ptr] | opt);
}

/*****************************************************************************/
