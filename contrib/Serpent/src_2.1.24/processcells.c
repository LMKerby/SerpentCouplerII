/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcells.c                                 */
/*                                                                           */
/* Created:       2010/10/11 (JLe)                                           */
/* Last modified: 2013/12/27 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Creates surface lists, etc.                                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessCells:"

/* Väliaikainen testifunktio */

void Testaa(long, long *, long);

/*****************************************************************************/

void ProcessCells()
{
  long cell, ptr, infix[10000], surf, loc0, nspec, nmax;

  fprintf(out, "Processing cells...\n");

  /* Avoid compiler warning */
  
  nmax = -1;

  /* Loop over cells */

  cell = RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /***********************************************************************/

      /***** Process surface lists according to format ***********************/

      /* Add counter */

      WDB[DATA_N_TOT_CELLS] = RDB[DATA_N_TOT_CELLS] + 1.0;

      /* Get parameters in infix format */

      nmax = ReadInfix(cell, infix, &nspec);

      /* Check special count (unions and parenthesis) */

      if (nspec == 0)
	{
	  /* Only intersections in list */
	  
	  IntersectionList(cell, infix, nmax);
	}
      else
	{
	  /* Convert list to postfix */
	  
	  ShuntingYard(cell, infix, nmax);

	  /* Add counter */

	  WDB[DATA_N_UNION_CELLS] = RDB[DATA_N_UNION_CELLS] + 1.0;

	  /* Test notation */
	  /*
	  Testaa(cell, infix, nmax);
	  */
	}

#ifndef TEST
    
      /***********************************************************************/
      
      /***** Link surfaces ***************************************************/

      /* Allocate memory for surface list */
	      
      loc0 = ReallocMem(DATA_ARRAY, nmax + 1);
      WDB[cell + CELL_PTR_SURF_LIST] = (double)loc0;
	      
      /* Check type */
      
      if ((ptr = RDB[cell + CELL_PTR_SURF_COMP]) > VALID_PTR)
	{
	  /*******************************************************************/

	  /***** Composition list given **************************************/
      
	  /* Loop over list */

	  while ((long)RDB[ptr] != 0)
	    {
	      /* Check type */

	      if ((long)RDB[ptr] > 0)
		{
		  /* Loop over surfaces */
		  
		  surf = RDB[DATA_PTR_S0];
		  while (surf > VALID_PTR)
		    {
		      /* Compare */
		      
		      if (!strcmp(GetText(ptr),
				  GetText(surf + SURFACE_PTR_NAME)))
			{
			  /* Put pointers */
			  
			  WDB[ptr] = (double)surf;
			  WDB[loc0++] = (double)surf;
			  
			  /* Break loop */
			  
			  break;
			}
		      
		      /* Next */
		      
		      surf = NextItem(surf);
		    }
		  
		  /* Check */
		  
		  if (surf < 0)
		    Error(cell, "Surface %s is not defined", GetText(ptr));
		}
	  
	      /* Next */

	      ptr++;
	    }

	  /* Put null terminator for surface list */
	  
	  WDB[loc0] = NULLPTR;

	  /*******************************************************************/
	}
      else if ((ptr = RDB[cell + CELL_PTR_SURF_INSC]) > VALID_PTR)
	{
	  /*******************************************************************/

	  /***** Intersection list given *************************************/

	  /* Close list */
	  
	  CloseList(ptr);

	  /* Loop over items */
  
	  while (ptr > VALID_PTR)
	    {
	      /* Loop over surfaces */
		  
	      surf = RDB[DATA_PTR_S0];
	      while (surf > VALID_PTR)
		{
		  /* Compare */
		  
		  if (CompareStr(surf + SURFACE_PTR_NAME, 
				 ptr + CELL_INSC_PTR_SURF))
		    {
		      /* Put pointers */
			  
		      WDB[ptr + CELL_INSC_PTR_SURF] = (double)surf;
		      WDB[loc0++] = (double)surf;

		      /* Break loop */
		      
		      break;
		    }
		  
		  /* Next */
		      
		  surf = NextItem(surf);
		}
		  
	      /* Check */
		  
	      if (surf < 0)
		Error(cell, "Surface %s is not defined", 
		      GetText(ptr + CELL_INSC_PTR_SURF));
	    
	      /* Next */

	      ptr = NextItem(ptr);
	    }
	  
	  /* Put null terminator for surface list */
	  
	  WDB[loc0] = (double)NULLPTR;
	  
	  /*******************************************************************/
	}
      else
	Die(FUNCTION_NAME, "No surface list");

      /***********************************************************************/

#endif
      
      /* Next cell */
      
      cell = NextItem(cell);
    }

  fprintf(out, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/

void Testaa(long cell, long *infix, long ni)
{
  long n, ptr, stack[10000], ns, a, b;
	
  fprintf(out, "\ncell %s: \n", GetText(cell + CELL_PTR_NAME));
  
  /* Print infix */
  
  fprintf(out, "infix:   ");
  
  for (n = 0; n < ni; n++)
    {
      if (infix[n] == SURF_OP_OR)
	fprintf(out, " + ");
      else if (infix[n] == SURF_OP_AND)
	fprintf(out, " * ");
      else if (infix[n] == SURF_OP_NOT)
	fprintf(out, " - ");
      else if (infix[n] == SURF_OP_LEFT)
	fprintf(out, " ( ");
      else if (infix[n] == SURF_OP_RIGHT)
	fprintf(out, " ) ");
      else
	{
	  WDB[DATA_DUMMY] = (double)infix[n];
	  fprintf(out, " %s ", GetText(DATA_DUMMY));
	}
    }
	
  /* Print postfix */
      
  fprintf(out, "\n");
  
  fprintf(out, "postfix: ");
	      
  ptr = RDB[cell + CELL_PTR_SURF_COMP];	  
  while ((long)RDB[ptr] != 0)
    {
      if ((long)RDB[ptr] == SURF_OP_OR)
	fprintf(out, " + ");
      else if ((long)RDB[ptr] == SURF_OP_AND)
	fprintf(out, " * ");
      else if ((long)RDB[ptr] == SURF_OP_NOT)
	fprintf(out, " - ");
      else if ((long)RDB[ptr] == SURF_OP_LEFT)
	fprintf(out, " ( ");
      else if ((long)RDB[ptr] == SURF_OP_RIGHT)
	fprintf(out, " ) ");
      else
	fprintf(out, " %s ", GetText(ptr));
      
      /* Next */
      
      ptr++;
    }

  fprintf(out, "\n");
	      
  /* Evaluate */
	      
  ns = 0;
  
  ptr = RDB[cell + CELL_PTR_SURF_COMP];	  
  while ((long)RDB[ptr] != 0)
    {
      if ((long)RDB[ptr] == SURF_OP_OR)
	{
	  /* Pop last two values */
	  
	  a = stack[--ns];
	  b = stack[--ns];
	  
	  /* Push union to stack */
	  
	  stack[ns++] = a + b;
	}
      else if ((long)RDB[ptr] == SURF_OP_NOT)
	{
	  /* Pop last value */
	  
	  a = stack[--ns];
	  
	  /* Push complement to stack */
	  
	  stack[ns++] = -a;
	}
      else if ((long)RDB[ptr] == SURF_OP_AND)
	{
	  /* Pop last two values */
	  
	  a = stack[--ns];
	  b = stack[--ns];
	  
	  /* Push intersection to stack */
	  
	  stack[ns++] = a * b;
	}
      else
	{
	  /* Push value to stack */
	  
	  stack[ns++] = atoi(GetText(ptr));
	}
      
      ptr++;
    }
  
  fprintf(out, "result:   %ld\n", stack[0]);
}
