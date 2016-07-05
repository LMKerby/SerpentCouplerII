/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getlinenumber.c                                */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2013/08/01 (JLe)                                           */
/* Version:       2.1.16                                                     */
/*                                                                           */
/* Description: Line number from position                                    */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetLineNumber:"

/*****************************************************************************/

long GetLineNumber(char *file, long i0)
{
  long nl, i;

  /* Check if blocked */

  if ((long)RDB[DATA_NO_LINE_NUMBER_CALC] == YES)
    return -1;

  /* Reset number of lines */

  nl = 1;

  /* Loop to position */

  for (i = 0; i < i0; i++)
    {
      /* Check eof */

      if (file[i] == '\0')
	Die(FUNCTION_NAME, "End-of-file encountered");

      /* Check newline and add counter */

      else if (file[i] == '\n')
	nl++;
    }

  /* Return line number */

  return nl;
}

/*****************************************************************************/
