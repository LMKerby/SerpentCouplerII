/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdivisors.c                              */
/*                                                                           */
/* Created:       2012/05/10 (JLe)                                           */
/* Last modified: 2014/03/29 (JLe)                                           */
/* Version:       2.1.19                                                     */
/*                                                                           */
/* Description: Processes material burnup divisors                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessDivisors:"

/*****************************************************************************/

void ProcessDivisors()
{
  long div, mat;
  double tmp;

  /* Reset maximum separation level */

  WDB[DATA_MAX_DIV_SEP_LVL] = -1.0;

  /* Loop over divisors */

  div = RDB[DATA_PTR_DIV0];
  while (div > VALID_PTR)
    {     
      /* Swap axial boundaries */

      if (RDB[div + DIV_ZMIN] > RDB[div + DIV_ZMAX])
	{
	  tmp = RDB[div + DIV_ZMIN];
	  WDB[div + DIV_ZMIN] = RDB[div + DIV_ZMAX];
	  WDB[div + DIV_ZMAX] = tmp;
	}
      else if (RDB[div + DIV_ZMIN] == RDB[div + DIV_ZMAX])
	Error(div, "Zero height");

      /* Swap radial boundaries */

      if (RDB[div + DIV_RMIN] > RDB[div + DIV_RMAX])
	{
	  tmp = RDB[div + DIV_RMIN];
	  WDB[div + DIV_RMIN] = RDB[div + DIV_RMAX];
	  WDB[div + DIV_RMAX] = tmp;
	}
      else if (RDB[div + DIV_RMIN] == RDB[div + DIV_RMAX])
	Error(div, "Zero thickness");

      /* Compare levels */

      if (((long)RDB[div + DIV_SEP] == YES) &&
	  (RDB[div + DIV_SEP_LVL] > RDB[DATA_MAX_DIV_SEP_LVL])) 
	WDB[DATA_MAX_DIV_SEP_LVL] = RDB[div + DIV_SEP_LVL];

      /* Convert tilt angle to rad */

      WDB[div + DIV_SEG0] = PI*RDB[div + DIV_SEG0]/180.0;

      /* Link material pointers */

      mat = RDB[DATA_PTR_M0];
      if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, 
			     GetText(div + DIV_PTR_MAT))) > VALID_PTR)
	{
	  /* Put pointers */

	  WDB[div + DIV_PTR_MAT] = (double)mat;
	  WDB[mat + MATERIAL_PTR_DIV] = (double)div;
	}
      else
	Error(div, "Material %s is not defined", GetText(div + DIV_PTR_MAT));

      /* Next divisor */

      div = NextItem(div);
    }
}

/*****************************************************************************/
