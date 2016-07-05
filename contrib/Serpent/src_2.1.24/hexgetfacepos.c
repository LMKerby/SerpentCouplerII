/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexgetfacepos.c                                */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Gets face positions of a hexahedron based on face direction  */
/*                                                                           */
/* Comments:   -Positions for face on right are 1 2 6 5 etc.                 */
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

#define FUNCTION_NAME "HexGetFacePos:"

/*****************************************************************************/

void HexGetFacePos(long *pos, long facedirection)
{
  switch (facedirection)
    {

    case (BOTTOMFACE):

      pos[0] = 0;
      pos[1] = 1;
      pos[2] = 2;
      pos[3] = 3;

      break;

    case (TOPFACE):

      pos[0] = 4;
      pos[1] = 5;
      pos[2] = 6;
      pos[3] = 7;

      break;

    case (LEFTFACE):

      pos[0] = 0;
      pos[1] = 3;
      pos[2] = 7;
      pos[3] = 4;

      break;

    case (RIGHTFACE):

      pos[0] = 1;
      pos[1] = 2;
      pos[2] = 6;
      pos[3] = 5;

      break;

    case (FRONTFACE):

      pos[0] = 0;
      pos[1] = 4;
      pos[2] = 5;
      pos[3] = 1;

      break;

    case (BACKFACE):

      pos[0] = 2;
      pos[1] = 3;
      pos[2] = 7;
      pos[3] = 6;

      break;

    default:
      Die(FUNCTION_NAME, "Invalid face index %ld\n", facedirection);      

    }
}
