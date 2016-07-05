/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexgetface.c                                   */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Gets face points from a hexahedron based on face direction   */
/*                                                                           */
/* Comments:   -Positions for face on right are hex[1] hex[2] hex[6] hex[5]  */
/*                                                                           */
/*             7-----6                                                       */
/*            /|    /|                                                       */
/*           4-----5 |                                                       */
/*           | 3---|-2                                                       */
/*           |/    |/                                                        */
/*           0-----1                                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "HexGetFace:"

/*****************************************************************************/

void HexGetFace(long hex[8], long *face, long facedirection)
{
  long pos[4], i;

  /* Check direction value */

  CheckValue(FUNCTION_NAME, "face direction", "", facedirection, 0, 5);

  /* Get face positions */

  HexGetFacePos(pos, facedirection);

  /* Store corresponding points from hex to face */

  for (i = 0; i < 4; i++)
    face[i] = hex[pos[i]];

}
