#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : eventtobank.c                                  */
/*                                                                           */
/* Created:       2014/09/29 (JLe)                                           */
/* Last modified: 2014/10/04 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Returns event to bank                                        */
/*                                                                           */
/* Comments: - NOTE: Tätä pitää kutsua OpenMP barrierin takaa                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "EventToBank:"

/*****************************************************************************/

void EventToBank(long ptr)
{
  long next;

  /* Reset data */
    
  memset(&WDB[ptr + LIFO_LIST_DATA_SIZE], 0.0, 
	 (EVENT_BLOCK_SIZE - LIFO_LIST_DATA_SIZE)*sizeof(double));

  /* Pointer to next */
  
  next = (long)RDB[DATA_PTR_EVENT_BANK];

  /* Put pointers */

  WDB[DATA_PTR_EVENT_BANK] = (double)ptr;
  WDB[ptr + LIFO_LIST_PTR_NEXT] = (double)next;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
