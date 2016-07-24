#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processevents.c                                */
/*                                                                           */
/* Created:       2011/10/01 (JLe)                                           */
/* Last modified: 2014/10/02 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Removes unused events from killed particles, allocates more  */
/*              memory for bank, etc.                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessEvents:"

/*****************************************************************************/

void ProcessEvents()
{
  long ptr, part, loc0, prev, count, n, np;

  /* Avoid compiler warning without debug mode */

  count = 0;
  part = 0;

  /* Check if events are set */

  if ((long)RDB[DATA_EVENT_RECORD_FLAGS] == 0)
    return;

#ifdef DEBUG

  /***************************************************************************/

  /***** Check history counts ************************************************/

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Reset previous */

      n = 1;

      /* Loop over events */
      
      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
	{
	  /* Check */

	  if ((long)RDB[ptr + EVENT_HIS_COUNT] < n)
	    Die(FUNCTION_NAME, "Something wrong here");
	  else
	    n = (long)RDB[ptr + EVENT_HIS_COUNT];

	  /* Next */

	  ptr = NextItem(ptr);
	}

      /* Next particle */

      part = NextItem(part);
    }

  /***************************************************************************/

  /***** Confirm count *******************************************************/

  /* Reset counter */

  count = 0;

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */
  
  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Loop over events */
      
      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
	{
	  /* Check if event was already counted */
	  
	  if ((long)RDB[ptr + EVENT_HIS_COUNT] > 0)
	    {
	      /* Add to count */

	      count++;

	      /* Mark as counted */

	      WDB[ptr + EVENT_HIS_COUNT] = 
		-RDB[ptr + EVENT_HIS_COUNT];
	    }
	  else if ((long)RDB[ptr + EVENT_HIS_COUNT] == 0)
	    Die(FUNCTION_NAME, "WTF?");

	  /* Next */

	  ptr = NextItem(ptr);
	}

      /* Next particle */

      part = NextItem(part);
    }

  /* Check that all events are accounted for */

  if (count + LIFOListSize(DATA_PTR_EVENT_BANK) != 
      (long)RDB[DATA_EVENT_BANK_SZ])
    Die(FUNCTION_NAME, "Mismatch in number of events");
  
  /***************************************************************************/

  /***** Set positive values back to counters ********************************/

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Loop over events */
      
      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
	{
	  /* Check if negative and switch to positive */
	  
	  if ((long)RDB[ptr + EVENT_HIS_COUNT] < 0)
	    WDB[ptr + EVENT_HIS_COUNT] = -RDB[ptr + EVENT_HIS_COUNT];

	  /* Next */

	  ptr = NextItem(ptr);
	}

      /* Next particle */

      part = NextItem(part);
    }

#endif

  /***************************************************************************/

  /***** Remove old events ***************************************************/

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Loop over particles */

  while (part > VALID_PTR)
    {
      /* Reset count */

      n = 0;

      /* Loop over events and update history counters */
 
      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
	{
	  /* Check if fission and update counter */

	  if ((long)RDB[ptr + EVENT_TYPE] == EVENT_TYPE_FISS)
	    n++;

	  /* Check with limit */
	  
	  if (n > (long)RDB[DATA_EVENT_MAX_GEN] - 1)
	    {
	      /* Update history count */

	      WDB[ptr + EVENT_HIS_COUNT] = RDB[ptr + EVENT_HIS_COUNT] - 1.0;
	    }

	  /* Pointer to next */

	  ptr = NextItem(ptr);
	}
    
      /* Loop over events and remove unused */
 
      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
	{
	  /* Check history counter */
	
	  if ((long)RDB[ptr + EVENT_HIS_COUNT] == 0)
	    {
	      /* Remove item */
	      
	      loc0 = ptr;
	      
	      /* Pointer to next */
	      
	      ptr = NextItem(loc0);
	      
	      /* Back to bank */

	      EventToBank(loc0);

	      /* Mark as removed */

	      WDB[loc0 + EVENT_HIS_COUNT] = -2E+6;
	    }
	  else if ((long)RDB[ptr + EVENT_HIS_COUNT] > 0)
	    {
	      /* Pointer to next */
	      
	      ptr = NextItem(ptr);
	    }
	  else if ((long)RDB[ptr + EVENT_HIS_COUNT] == -2E+6)
	    break;
	  else
	    Die(FUNCTION_NAME, "WTF?");
	}

      /* Next particle */

      part = NextItem(part);
    } 

  /* Pointer to source */

  part = (long)RDB[DATA_PART_PTR_SOURCE];
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to first */

  part = FirstItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
  
  /* Get pointer to first item after dummy */

  part = NextItem(part);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);
 
  /* Loop over particles */

  while (part > VALID_PTR)
    {
       /* Loop over events and detach pointers */

      prev = -1;
      
      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
	{
	  /* Check if event was removed */

	  if ((long)RDB[ptr + EVENT_HIS_COUNT] == -2E+6)
	    {
	      /* Detach previous */

	      if (prev < VALID_PTR)
		WDB[part + PARTICLE_PTR_EVENTS] = NULLPTR;
	      else
		WDB[prev + LIFO_LIST_PTR_NEXT] = NULLPTR;

	      /* Break loop */

	      break;
	    }

	  /* Copy pointer */

	  prev = ptr;

	  /* Pointer to next */

	  ptr = NextItem(ptr);
	}

      /* Check history counts */

      ptr = (long)RDB[part + PARTICLE_PTR_EVENTS];
      while (ptr > VALID_PTR)
	{
	  if ((long)RDB[ptr + EVENT_HIS_COUNT] < 1)
	    Die(FUNCTION_NAME, "Error in history count");

	  /* Pointer to next */

	  ptr = NextItem(ptr);
	}
      
      /* Next particle */

      part = NextItem(part);
    }

  /***************************************************************************/

  /***** Allocate memory for more events if needed ***************************/

  /* Pointer to event bank */

  if ((ptr = (long)RDB[DATA_PTR_EVENT_BANK]) > VALID_PTR)
    {
      /* Check size */

      if (((double)LIFOListSize(DATA_PTR_EVENT_BANK))
	  /RDB[DATA_EVENT_BANK_SZ] < 0.2)
	{
	  /* Print */

	  Note(0, "Adjusting event bank size...");

	  /* Set new size */

	  np = (long)(0.2*RDB[DATA_EVENT_BANK_SZ]);
	  
	  /* Allow memory allocation */
	  
	  Mem(MEM_ALLOW);
	  
	  /* Allocate memory for events */
	  
	  for (n = 0; n < np; n++)
	    NewLIFOItem(DATA_PTR_EVENT_BANK, EVENT_BLOCK_SIZE);
	  
	  /* Deny memory allocation */
	  
	  Mem(MEM_DENY);
	  
	  /* Put new size */
	  
	  WDB[DATA_EVENT_BANK_SZ] = RDB[DATA_EVENT_BANK_SZ] + (double)np;
	}
    }

  /***************************************************************************/

}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
