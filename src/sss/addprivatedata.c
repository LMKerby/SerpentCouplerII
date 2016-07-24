#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addprivatedata.c                               */
/*                                                                           */
/* Created:       2011/11/11 (JLe)                                           */
/* Last modified: 2014/04/04 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Adds value in PRIVA data block                               */
/*                                                                           */
/* Comments: - Replaced by macro when not compiled in debug mode             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddPrivateData:"

/*****************************************************************************/

#ifdef DEBUG 

void AddPrivateData(long ptr, double val, long id)
{
  long sz;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);

  /* Check if access is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "PRIVA array not ready for access");

  /* Get size of data block */
  
  sz = (long)RDB[DATA_REAL_PRIVA_SIZE];
  
  /* Check pointer */

  if ((ptr < 0) || (ptr > sz - 1))
    Die(FUNCTION_NAME, "Pointer error");

  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");
  
  /* Put value */
  
  PRIVA[ptr + id*sz] = PRIVA[ptr + id*sz] + val;
}

#endif

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
