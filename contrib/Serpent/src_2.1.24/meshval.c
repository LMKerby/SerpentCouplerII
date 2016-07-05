/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : meshval.c                                      */
/*                                                                           */
/* Created:       2011/05/14 (JLe)                                           */
/* Last modified: 2012/01/19 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description:  Returns value in mesh by coordinates                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MeshVal:"

/*****************************************************************************/

double MeshVal(long msh, double x, double y, double z)
{
  long ptr, idx;
  double val;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Check content */

  if ((long)RDB[msh + MESH_CONTENT] != MESH_CONTENT_DATA)
    Die(FUNCTION_NAME, "Invalid content type");

  /* Get mesh index */

  if ((idx = MeshIndex(msh, x, y, z)) < 0)
    return 0.0;

  /* Get pointer to data */
  
  ptr = (long)RDB[msh + MESH_PTR_RES2];
  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

  /* Get value */

  val = SumPrivateRes(ptr + idx + 1);
  
  /* Return value */
  
  return val;
}

/*****************************************************************************/
