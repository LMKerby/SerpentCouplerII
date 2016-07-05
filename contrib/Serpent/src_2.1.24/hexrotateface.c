/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexrotateface.c                                */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Rotates a face consisting of four points                     */
/*                                                                           */
/* Comments:   -Rotation of face 0-1-2-3 by +1 gives 3-0-1-2                 */
/*             -Rotation of face 0-1-2-3 by -1 gives 1-2-3-0                 */
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

#define FUNCTION_NAME "FixHexMesh:"

/*****************************************************************************/

void HexRotateFace(long *face, long steps)
{
  long temp;
  long i;


  if (steps > 0)
    {

      /* Positive rotation */
      /* Rotate for requested number of steps */

      for (i = 0; i < steps; i++)
	{
	  temp = face[3];
	  face[3] = face[2];
	  face[2] = face[1];
	  face[1] = face[0];
	  face[0] = temp;
	}
    }
  else
    {
      /* Negative rotation */

      steps = -steps;

      /* Rotate for requested number of steps */

      for (i = 0; i < steps; i++)
	{
	  temp = face[0];
	  face[0] = face[1];
	  face[1] = face[2];
	  face[2] = face[3];
	  face[3] = temp;
	}
    }
}
