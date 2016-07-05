/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : toque.c                                        */
/*                                                                           */
/* Created:       2011/03/09 (JLe)                                           */
/* Last modified: 2014/02/22 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Puts neutron / photon in que                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ToQue:"

/*****************************************************************************/

void ToQue(long ptr, long id)
{
  /* Check id */

  if ((id < 0) || (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1))
    Die(FUNCTION_NAME, "Error in thread id");

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Add item in list */

  AddItem(OMPPtr(DATA_PART_PTR_QUE, id), ptr);

  /* Sort list to transport neutrons before photons */

  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_QUE, id)];
  SortList(ptr, PARTICLE_TYPE, SORT_MODE_ASCEND);
}

/*****************************************************************************/
