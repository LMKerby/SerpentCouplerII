#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : atof.c                                         */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2014/03/18 (JLe)                                           */
/* Version:       2.1.19                                                      */
/*                                                                           */
/* Description: Array to double conversion with type checking                */
/*                                                                           */
/* Comments: - From Serpent 1.1.8                                            */
/*           - Tätä pitäisi kutsua pelkästään GetParam():sta.                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "AtoF:"

/*****************************************************************************/

double AtoF(char *str, char *param, char *file, long line)
{
  long n;

  /* Loop over string */

  for (n = 0; n < (long)strlen(str); n++)
    {
      /* Check that entry is numerical */

      if (!isdigit(str[n]) && (str[n] != '.') && (str[n] != '-') && 
	  (str[n] != '+') && (str[n] != 'e') && (str[n] != 'E'))
	{
	  /* Check mode */
	  
	  if (line > 0)
	    {
	      /* Error in input file */
	      
	      Error(-1, param, file, line, "Invalid numerical entry \"%s\"",
		    str);
	    }
	  else
	    {
	      /* Error in string */
	      
	      Die(FUNCTION_NAME, "Invalid numerical entry \"%s\"", str);
	    }
	}

      /* Check for common shorthand notation */

      if (n > 0)
	if (str[n] == '-')
	  if ((str[n - 1] != 'e') && (str[n - 1] != 'E'))
	    Error(-1, param, file, line, "Invalid numerical entry \"%s\"", 
		  str);
    }

  /* Return number */

  return atof(str);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
