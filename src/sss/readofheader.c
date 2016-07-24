#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofheader.c                                 */
/*                                                                           */
/* Created:       2013/12/27 (JLe)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Reads header data from an OpenFOAM format file               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ReadOFHeader:"

/*****************************************************************************/

void ReadOFHeader(FILE *fp, long *type, long *sz, long *dim)
{
  long n;
  char tmpstr[MAX_STR], c;
  
  /* Check file pointer */

  if (fp == NULL)
    Die(FUNCTION_NAME, "Pointer error");

  /* Reset type */

  *type = -1;

  /* Reset dimension vector */

  for (n = 0; n < 7; n++)
    dim[n] = 0;

  /* Read first string */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    Die(FUNCTION_NAME, "Unexpected EOF");

  /* Check size */

  if ((n = atol(tmpstr)) > 0)
    {
      /* Put size */

      *sz = n;

      /* Exit subroutine */

      return;
    }
  
  /* Look for key word "object" */
  
  do 
    {
      /* Check string */

      if (!strcmp(tmpstr, "object"))
	{
	  /* Read next word */
	  
	  if (fscanf(fp, "%s", tmpstr) == EOF)
	    Die(FUNCTION_NAME, "Unexpected EOF");
	  
	  /* Check and set type */
	  
	  if (!strcmp(tmpstr, "points;"))
	    *type = OF_FILE_POINTS;
	  else if (!strcmp(tmpstr, "faces;"))
	    *type = OF_FILE_FACES;
	  else if (!strcmp(tmpstr, "owner;"))
	    *type = OF_FILE_OWNER;
	  else if (!strcmp(tmpstr, "neighbour;"))
	    *type = OF_FILE_NEIGHBOUR;
	  else if ((!strcmp(tmpstr, "rho;")) || 
		   (!strcmp(tmpstr, "rhok;")))
	    *type = OF_FILE_DENSITY;
	  else if (!strcmp(tmpstr, "T;"))
	    *type = OF_FILE_TEMP;
	  
	  /* Break loop */
	  
	  break;
	}
    }
  while (fscanf(fp, "%s", tmpstr) != EOF);

  /* Loop until end of block */

  while ((n = fscanf(fp, "%s", tmpstr)) != EOF)
    if (!strcmp(tmpstr, "}"))
      break;

  /* Check */

  if (n == EOF)
    Die(FUNCTION_NAME, "Unexpected EOF");

  /* Read lines until size is found */

  while (1 != 2)
    {
      if (fgets(tmpstr, MAX_STR, fp) == NULL)
	Die(FUNCTION_NAME, "Unexpected EOF");

      /* Skip comment lines */

      if ((tmpstr[0] == '/') && (tmpstr[1] == '/'))
	continue;

      /* Check dimensions array */

      if (!strncmp(tmpstr, "dimensions", 10))
	{
	  /* Read dimensions */

	  sscanf(tmpstr, "dimensions [%ld %ld %ld %ld %ld %ld %ld]", &dim[0], 
		 &dim[1], &dim[2], &dim[3], &dim[4], &dim[5], &dim[6]);
	}
      else
	{
	  /* Convert to integer */

	  n = (long)atoi(tmpstr);

	  /* Check */

	  if (n > 0)
	    {
	      /* Put size */
	      
	      *sz = n;
	      
	      /* Loop until the first '(' */
	      
	      while ((c = fgetc(fp)) != EOF)
		if (c == '(')
		  break;
	      
	      /* Exit subroutine */
	      
	      return;
	    }
	}
    }

  /* Something wrong */

  Die(FUNCTION_NAME, "Shouldn't be here");
}

/*****************************************************************************/



#ifdef __cplusplus 
} 
#endif 
