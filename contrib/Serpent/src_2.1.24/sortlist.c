/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sortlist.c                                     */
/*                                                                           */
/* Created:       2010/09/18 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: - Sorts list items                                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SortList:"

/*****************************************************************************/

void SortList(long lst, long param, long mode)
{
  long ptr, ptp, loc0, next, count;
  double val1, val2;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);
  
  /* Check if access to private data is allowed */
  
  if ((mode == SORT_MODE_ASCEND_PRIVA) || (mode == SORT_MODE_DESCEND_PRIVA))
    if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
      Die(FUNCTION_NAME, "PRIVA array not ready for access");

  /* Go to first item */
  
  lst = FirstItem(lst);

  /* Sort loop */

  do
    {
      /* Reset count */

      count = 0;

      /* Check mode */

      if ((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_ASCEND_PRIVA))
	{
	  /***** Sort items in ascending order *******************************/

	  /* Set starting point for iteration */

	  loc0 = FirstItem(lst);

	  /* Sorting loop */
	  
	  while (loc0 > VALID_PTR)
	    {
	      /* Set pointers for current and next */

	      ptr = loc0;
	      next = NextItem(ptr);

	      /* Set starting point for next iteration */
	      
	      loc0 = next;

	      /* Loop over list */
	  
	      while (next > VALID_PTR)
		{
		  /* Get values */

		  if (mode == SORT_MODE_ASCEND)
		    {
		      /* Get direct values */

		      val1 = RDB[ptr + param];
		      val2 = RDB[next + param];
		    }
		  else
		    {
		      /* Get values from private block */

		      ptp = (long)RDB[ptr + param];
		      CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
		      val1 = SumPrivateData(ptp);

		      ptp = (long)RDB[next + param];
		      CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
		      val2 = SumPrivateData(ptp);
		    }

		  /* Compare items, move one step forward or break loop */
		  
		  if (val1 > val2)
		    {
		      MoveItemRight(ptr);
		      count++;
		    }
		  else
		    break;
		  
		  /* Next */
		  
		  next = NextItem(ptr);
		}
	    }

	  /*******************************************************************/
	}
      else if((mode == SORT_MODE_DESCEND) || (mode == SORT_MODE_DESCEND_PRIVA))
	{
	  /***** Sort items in descending order ******************************/

	  /* Set starting point for iteration */

	  loc0 = FirstItem(lst);

	  /* Sorting loop */
	  
	  while (loc0 > VALID_PTR)
	    {
	      /* Set pointers for current and next */

	      ptr = loc0;
	      next = NextItem(ptr);

	      /* Set starting point for next iteration */
	      
	      loc0 = next;

	      /* Loop over list */
	  
	      while (next > VALID_PTR)
		{
		  /* Get values */

		  if (mode == SORT_MODE_DESCEND)
		    {
		      /* Get direct values */

		      val1 = RDB[ptr + param];
		      val2 = RDB[next + param];
		    }
		  else
		    {
		      /* Get values from private block */
		      
		      ptp = (long)RDB[ptr + param];
		      CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
		      val1 = SumPrivateData(ptp);

		      ptp = (long)RDB[next + param];
		      CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
		      val2 = SumPrivateData(ptp);
		    }

		  /* Compare items, move one step forward or break loop */
		  
		  if (val1 < val2)
		    {
		      MoveItemRight(ptr);
		      count++;
		    }
		  else
		    break;
		  
		  /* Next */
		  
		  next = NextItem(ptr);
		}
	    }

	  /*******************************************************************/
	}
      else
	Die(FUNCTION_NAME, "Invalid sort mode");
    }
  while (count > 0);

#ifdef DEBUG

  /* Check order */

  ptr = FirstItem(lst);

  while (ptr > 0)
    {
      if ((next = NextItem(ptr)) > 0)
	{
	  /* Get values */
	  
	  if ((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_DESCEND))
	    {
	      /* Get direct values */
	      
	      val1 = RDB[ptr + param];
	      val2 = RDB[next + param];
	    }
	  else
	    {
	      /* Get values from private block */

	      ptp = (long)RDB[ptr + param];
	      CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
	      val1 = SumPrivateData(ptp);
	      
	      ptp = (long)RDB[next + param];
	      CheckPointer(FUNCTION_NAME, "(ptp)", PRIVA_ARRAY, ptp);
	      val2 = SumPrivateData(ptp);
	    }
	  
	  if (((mode == SORT_MODE_ASCEND) || (mode == SORT_MODE_ASCEND_PRIVA))
	      && (val1 > val2))
	    Die(FUNCTION_NAME, "Sort error");
	  else if (((mode == SORT_MODE_DESCEND) || 
		    (mode == SORT_MODE_DESCEND_PRIVA))
		   && (val1 < val2))
	    Die(FUNCTION_NAME, "Sort error");
	}
      else
	break;
      
      /* Next */

      ptr = next;
    }

#endif

  /* Set direct pointers if list is closed */

  if ((long)RDB[lst + LIST_PTR_DIRECT] > VALID_PTR)
    SetDirectPointers(lst);
}

/*****************************************************************************/
