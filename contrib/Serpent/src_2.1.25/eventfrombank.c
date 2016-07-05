/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : eventfrombank.c                                */
/*                                                                           */
/* Created:       2014/09/29 (JLe)                                           */
/* Last modified: 2014/10/04 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Retrieves a new event structure from bank                    */
/*                                                                           */
/* Comments: - Structure should be cleared when returned                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EventFromBank:"

/*****************************************************************************/

long EventFromBank(long part)
{
  long ptr, next;
  
#ifdef OPEN_MP
#pragma omp critical
#endif
  {
    /* Pointer to bank */

    if ((ptr = (long)RDB[DATA_PTR_EVENT_BANK]) < VALID_PTR)
      Die(FUNCTION_NAME, "Event bank is empty");

    /* Pointer to next */

    next = NextItem(ptr);

    /* Set pointer */

    WDB[DATA_PTR_EVENT_BANK] = (double)next;
  }

  /* Attach pointer to particle */
  
  WDB[ptr + LIFO_LIST_PTR_NEXT] = RDB[part + PARTICLE_PTR_EVENTS];
  WDB[part + PARTICLE_PTR_EVENTS] = (double)ptr;

  /* Set counter */
      
  WDB[ptr + EVENT_HIS_COUNT] = 1.0;

  /* Return pointer */
  
  return ptr;
}

/*****************************************************************************/
