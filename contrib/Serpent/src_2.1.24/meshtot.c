/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshtot.c                                      */
/*                                                                           */
/* Created:       2011/05/14 (JLe)                                           */
/* Last modified: 2012/01/17 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description:  Returns sum of all mesh values                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshTot:"

/*****************************************************************************/

double MeshTot(long msh)
{
  long ptr;
  double tot;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_DATA)
    Die(FUNCTION_NAME, "Invalid content type");
  
  /* Get pointer to data */
  
  ptr = (long)RDB[msh + MESH_PTR_RES2];
  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

  /* Get total */

  tot = SumPrivateRes(ptr);

  /* Return value */

  return tot;
}

/*****************************************************************************/
