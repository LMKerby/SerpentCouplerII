#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : Error.c                                        */
/*                                                                           */
/* Created:       2010/09/23 (JLe)                                           */
/* Last modified: 2013/08/01 (JLe)                                           */
/* Version:       2.1.16                                                     */
/*                                                                           */
/* Description: Prints input error message                                   */
/*                                                                           */
/* Comments: - Converted from Serpent 1.1.12                                 */
/*           - Use this for user errors, errors in code terminate the run    */
/*             with Die()                                                    */
/*           - This and Due() are the only functions that should make        */
/*             a call to exit().                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Error:"

/*****************************************************************************/

void Error(long ptr, ...)
{
  char param[MAX_STR], fname[MAX_STR];
  long line;
  va_list argp;
  va_start (argp, ptr);

  /* Initialize variables to avoid compiler warning */

  *param = '\0';
  *fname = '\0';
  line = 0;

  /* Check pointer */

  if (ptr > VALID_PTR)
    {
      /* Put parameter name */

      if ((long)RDB[ptr + PARAM_PTR_NAME] > VALID_PTR)
	strcpy(param, GetText(ptr + PARAM_PTR_NAME));
      else
	ptr = 0;

      /* Check pointer */

      if (ptr > VALID_PTR)
	{
	  /* Put file name */
	  
	  if ((long)RDB[ptr + PARAM_PTR_FNAME] > VALID_PTR)
	    strcpy(fname, GetText(ptr + PARAM_PTR_FNAME));
	  else
	    Die(FUNCTION_NAME, "File name not set");
	  
	  /* Put line number */
      
	  line = (long)RDB[ptr + PARAM_LINE];
	}
    }
  else if (ptr < 0)
    {
      /* Use argument values */

      strcpy(param, va_arg(argp, char *));
      strcpy(fname, va_arg(argp, char *));
      line = va_arg(argp, int);
    }

  /* Print error message */

  fprintf(stdout, "\n***** %s\n\n", TimeStamp());

  if (ptr != 0)
    {
      if (line > -1)
	fprintf(stdout, "Input error in parameter \"%s\" on line %ld ", param, 
		line);
      else
	fprintf(stdout, "Input error in parameter \"%s\" ", param);

      fprintf(stdout, "in file \"%s\":\n\n", fname);
    }
  else
    fprintf(stdout, "Input error:\n\n");

  vfprintf(stdout, va_arg(argp, char *), argp);

  fprintf(stdout, "\n\n");

  /* Exit subroutine */

  exit(-1);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
