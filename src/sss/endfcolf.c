#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : endfcolf.c                                     */
/*                                                                           */
/* Created:       2010/11/20 (JLe)                                           */
/* Last modified: 2013/12/23 (JLe)                                           */
/* Version:       2.1.16                                                     */
/*                                                                           */
/* Description: Reads double from column formatted ENDF file                 */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ENDFColF:"

/*****************************************************************************/

double ENDFColF(long col, char *line)
{
  long i, n, y;
  double x;
  char str[MAX_STR], c;

  /* Check column */

  if ((col < 1) || (col > 6))
    Die(FUNCTION_NAME, "Invalid column index %d", col);

  /* Try normal transformation first */
  
  n = 11*(col - 1);

  for (i = 0; i < 11; i++)
    if (((c = line[n + i]) == 'e') || (c == 'E'))
      {
	/* Format OK. Transform to float */
	
	strcpy(str, &line[n]);
	str[11] = '\0';

	/* Return value */

	return atof(str);
      }

  /* Copy first part of string */

  n = 11*(col - 1);

  i = 0;
  while(((c = line[n++]) != '+') && (c != '-') && (i < 11))
    str[i++] = c;

  str[i] = '\0';

  /* Get float value */

  x = atof(str);

  /* Check length */

  if (i == 11)
    return x;

  /* Copy remaining part */

  strcpy(str, &line[n-1]);
  str[11 - i] = '\0';

  /* Get exponent */

  y = atoi(str);

  /* Make new string */

  sprintf(str, "%1.6fE%ld", x, y);

  /* Return value */
  
  return atof(str);
}

/*****************************************************************************/



#ifdef __cplusplus 
} 
#endif 
