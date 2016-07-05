/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexrotatelists.c                               */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Rotates the auxilliary directional lists for hex cells       */
/*              after the hex itself has been rotated to final position      */
/*                                                                           */
/* Comments:   -Needed for fixhexmesh.c                                      */
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

#define FUNCTION_NAME "HexRotateLists:"

/*****************************************************************************/

void HexRotateLists(long hex[8], long (*initFaces)[4], long *hexfaces, long *hexsides, long (*hexnbrs)[6])
{
  long i, j;
  long tmpfaces[6], tmpsides[6], tmpnbrs[2][6], face0[4], *face1;

  /* Store all the lists to temporary lists */

  for (i = 0; i < 6; i++)
    {

      tmpfaces[i] = hexfaces[i];
      tmpsides[i] = hexsides[i];
      tmpnbrs[0][i]  = hexnbrs[0][i];
      tmpnbrs[1][i]  = hexnbrs[1][i];

    }

  /* Loop over face directions */

  for (i = 0; i < 6; i++)
    {

      /* Get next face from hex */

      HexGetFace(hex, face0, i);

      /* Loop over initial faces and find match */

      for (j = 0; j < 6; j++)
	{
	  face1 = initFaces[j];

	  /* Check if faces match */

	  if (PolySameFace(face0,face1, 4))
	    {
	      /* If faces matched store list data for this face */

	      hexfaces[i] = tmpfaces[j];
	      hexsides[i] = tmpsides[j];
	      hexnbrs[0][i] = tmpnbrs[0][j];
	      hexnbrs[1][i] = tmpnbrs[1][j];

	      break;
	    }

	}

      /* If this rotated face did not match any of the initial faces */
      /* something is wrong in fixhexmesh.c                          */

      if (j == 6)
	Die(FUNCTION_NAME, "Could not match faces");
    }

}
