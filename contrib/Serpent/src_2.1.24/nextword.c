/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : nextword.c                                     */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Reads next white-space separated word from string            */
/*                                                                           */
/* Comments: - From Serpent 1.1.0                                            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "NextWord:"

/*****************************************************************************/

long NextWord(char *from, char *to)
{
  long n, i;

  /* skip leading white spaces */

  n = 0;
  while(((from[n] == ' ') || (from[n] == '\n')) && (from[n] != EOF))
    n++;

  /* check if string is in quotations */
  
  if (from[n] == '\"')
    {
      n++;

      /* read string */

      i = 0;

      while ((from[n] != '\"') && (from[n] != '\0'))
	to[i++] = from[n++];

      n++;
    }
  else
    {
      /* read string */

      i = 0;

      while((from[n] != ' ') && (from[n] != '\n') && (from[n] != '\0'))
	to[i++] = from[n++];
    }

  /* terminate string */
  
  to[i] = '\0';

  return n;
}

/*****************************************************************************/
