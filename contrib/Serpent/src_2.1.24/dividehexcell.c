/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fixhexmesh.c                                   */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Splits hexahedral meshes to tetrahedrons                     */
/*                                                                           */
/* Comments:                                                                 */
/*             -Based on "J. Dompierre, P. Labbe, M. Vallet and R. Camarero, */
/*              How to Subdivide Pyramids, Prisms, and Hexahedra             */ 
/*              into Tetrahedra, Proceedings,                                */
/*              8th International Meshing Roundtable,                        */
/*              South Lake Tahoe, CA, U.S.A., pp.195-204, October 1999       */
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

#define FUNCTION_NAME "DivideHexCell:"

/*****************************************************************************/

void DivideHexCell(long ifc, long cgns, long hex[8], long (*diags)[8], 
		   long (*hexnbrs)[6], long hexfaces[6], long hexsides[6], 
		   long *ownrlist2, long *nbrlist2, long *sdone)
{
  long nd, nc, i, j, idx, pts[6][4], directions[6][4];
  long p0, p1, p2, p3, firsts[6], intfirsts[6][6], tmpface[4];
  long ptr, ncgns, first, pointlist;

  /* Get pointer to parents pointlist */

  pointlist = (long)RDB[ifc + IFC_PTR_POINT_LIST_PRNTS];;
  CheckPointer(FUNCTION_NAME, "(pointlist)", DATA_ARRAY, pointlist);

  /* Count number of diagonals going through 6 */

  nd = 0;
  
  for (i = 0; i < 8; i++)
    if (diags[6][i] == 1)
      nd++;

  /* Get number of subcells based on number of diagonals */

  if (nd == 0)
    nc = 5;
  else
    nc = 6;

  /* Create the new cells */

  for (i = 0; i < nc; i++)
    {
      ncgns = NewItem(ifc + IFC_PTR_TET_MSH, IFC_TET_MSH_LIST_BLOCK_SIZE);
      CheckPointer(FUNCTION_NAME, "(ncgns)", DATA_ARRAY, ncgns);

      /* Get tet mesh index */

      if ((ptr = PrevItem(ncgns)) > VALID_PTR)
	idx = (long)RDB[ptr + IFC_TET_MSH_IDX] + 1;
      else
	idx = (long)RDB[ifc + IFC_NC_PRNTS];

      /* Put index */

      WDB[ncgns + IFC_TET_MSH_IDX] = (double)idx;

      /* Put number of faces */

      WDB[ncgns + IFC_TET_MSH_NF] = 4.0;

      /* Allocate memory for face list */

      ptr = ReallocMem(DATA_ARRAY, 4.0);

      /* Put pointer to face list */

      WDB[ncgns + IFC_TET_MSH_PTR_FACES] = (double)ptr;

      /* Allocate memory for side list */

      ptr = ReallocMem(DATA_ARRAY, 4.0);

      /* Put pointer to side list */

      WDB[ncgns + IFC_TET_MSH_PTR_SIDES] = (double)ptr;

      /* Put parent pointer */
      
      WDB[ncgns + IFC_TET_MSH_PTR_PARENT] = (double)cgns;
    }

  /* Get the first new cell */

  for (i = 0; i < nc - 1; i++)
    ncgns = PrevItem(ncgns);
    
  /* Get points and directions */

  if (nd == 0)
    {
      pts[0][0] = hex[0];
      pts[0][1] = hex[1];
      pts[0][2] = hex[2];
      pts[0][3] = hex[5];

      pts[1][0] = hex[0];
      pts[1][1] = hex[2];
      pts[1][2] = hex[7];
      pts[1][3] = hex[5];

      pts[2][0] = hex[0];
      pts[2][1] = hex[2];
      pts[2][2] = hex[3];
      pts[2][3] = hex[7];

      pts[3][0] = hex[0];
      pts[3][1] = hex[5];
      pts[3][2] = hex[7];
      pts[3][3] = hex[4];

      pts[4][0] = hex[2];
      pts[4][1] = hex[7];
      pts[4][2] = hex[5];
      pts[4][3] = hex[6];

      pts[5][0] = -1;
      pts[5][1] = -1;
      pts[5][2] = -1;
      pts[5][3] = -1;

      directions[0][0] = BOTTOMFACE;
      directions[0][1] = RIGHTFACE;
      directions[0][2] = FRONTFACE;
      directions[0][3] = -2;

      directions[1][0] = -3;
      directions[1][1] = -5;
      directions[1][2] = -1;
      directions[1][3] = -4;

      directions[2][0] = BOTTOMFACE;
      directions[2][1] = BACKFACE;
      directions[2][2] = -2;
      directions[2][3] = LEFTFACE;

      directions[3][0] = -2;
      directions[3][1] = TOPFACE;
      directions[3][2] = FRONTFACE;
      directions[3][3] = LEFTFACE;

      directions[4][0] = -2;
      directions[4][1] = TOPFACE;
      directions[4][2] = BACKFACE;
      directions[4][3] = RIGHTFACE;

      directions[5][0] = -100;
      directions[5][1] = -100;
      directions[5][2] = -100;
      directions[5][3] = -100;

    }
  else if (nd == 1)
    {
      pts[0][0] = hex[0];
      pts[0][1] = hex[5];
      pts[0][2] = hex[7];
      pts[0][3] = hex[4];

      pts[1][0] = hex[0];
      pts[1][1] = hex[1];
      pts[1][2] = hex[7];
      pts[1][3] = hex[5];

      pts[2][0] = hex[1];
      pts[2][1] = hex[6];
      pts[2][2] = hex[7];
      pts[2][3] = hex[5];

      pts[3][0] = hex[0];
      pts[3][1] = hex[7];
      pts[3][2] = hex[2];
      pts[3][3] = hex[3];

      pts[4][0] = hex[0];
      pts[4][1] = hex[7];
      pts[4][2] = hex[1];
      pts[4][3] = hex[2];

      pts[5][0] = hex[1];
      pts[5][1] = hex[7];
      pts[5][2] = hex[6];
      pts[5][3] = hex[2];

      directions[0][0] = -2;
      directions[0][1] = TOPFACE;
      directions[0][2] = FRONTFACE;
      directions[0][3] = LEFTFACE;

      directions[1][0] = -5;
      directions[1][1] = -3;
      directions[1][2] = FRONTFACE;
      directions[1][3] = -1;

      directions[2][0] = -6;
      directions[2][1] = TOPFACE;
      directions[2][2] = RIGHTFACE;
      directions[2][3] = -2;

      directions[3][0] = -5;
      directions[3][1] = BACKFACE;
      directions[3][2] = LEFTFACE;
      directions[3][3] = BOTTOMFACE;

      directions[4][0] = -2;
      directions[4][1] = -6;
      directions[4][2] = -4;
      directions[4][3] = BOTTOMFACE;

      directions[5][0] = -3;
      directions[5][1] = BACKFACE;
      directions[5][2] = -5;
      directions[5][3] = RIGHTFACE;

    }
  else if (nd == 2)
    {
      pts[0][0] = hex[0];
      pts[0][1] = hex[4];
      pts[0][2] = hex[5];
      pts[0][3] = hex[6];

      pts[1][0] = hex[0];
      pts[1][1] = hex[3];
      pts[1][2] = hex[7];
      pts[1][3] = hex[6];

      pts[2][0] = hex[0];
      pts[2][1] = hex[7];
      pts[2][2] = hex[4];
      pts[2][3] = hex[6];

      pts[3][0] = hex[0];
      pts[3][1] = hex[1];
      pts[3][2] = hex[2];
      pts[3][3] = hex[5];

      pts[4][0] = hex[0];
      pts[4][1] = hex[3];
      pts[4][2] = hex[6];
      pts[4][3] = hex[2];

      pts[5][0] = hex[0];
      pts[5][1] = hex[6];
      pts[5][2] = hex[5];
      pts[5][3] = hex[2];

      directions[0][0] = FRONTFACE;
      directions[0][1] = TOPFACE;
      directions[0][2] = -3;
      directions[0][3] = -6;

      directions[1][0] = LEFTFACE;
      directions[1][1] = BACKFACE;
      directions[1][2] = -5;
      directions[1][3] = -3;

      directions[2][0] = LEFTFACE;
      directions[2][1] = TOPFACE;
      directions[2][2] = -2;
      directions[2][3] = -1;

      directions[3][0] = BOTTOMFACE;
      directions[3][1] = RIGHTFACE;
      directions[3][2] = FRONTFACE;
      directions[3][3] = -6;

      directions[4][0] = -2;
      directions[4][1] = BACKFACE;
      directions[4][2] = BOTTOMFACE;
      directions[4][3] = -6;

      directions[5][0] = -1;
      directions[5][1] = RIGHTFACE;
      directions[5][2] = -5;
      directions[5][3] = -4;
    }
  else if (nd == 3)
    {
      pts[0][0] = hex[0];
      pts[0][1] = hex[2];
      pts[0][2] = hex[3];
      pts[0][3] = hex[6];

      pts[1][0] = hex[0];
      pts[1][1] = hex[3];
      pts[1][2] = hex[7];
      pts[1][3] = hex[6];

      pts[2][0] = hex[0];
      pts[2][1] = hex[7];
      pts[2][2] = hex[4];
      pts[2][3] = hex[6];

      pts[3][0] = hex[0];
      pts[3][1] = hex[5];
      pts[3][2] = hex[6];
      pts[3][3] = hex[4];

      pts[4][0] = hex[1];
      pts[4][1] = hex[5];
      pts[4][2] = hex[6];
      pts[4][3] = hex[0];

      pts[5][0] = hex[1];
      pts[5][1] = hex[6];
      pts[5][2] = hex[2];
      pts[5][3] = hex[0];

      directions[0][0] = BOTTOMFACE;
      directions[0][1] = BACKFACE;
      directions[0][2] = -6;
      directions[0][3] = -2;

      directions[1][0] = LEFTFACE;
      directions[1][1] = BACKFACE;
      directions[1][2] = -1;
      directions[1][3] = -3;

      directions[2][0] = LEFTFACE;
      directions[2][1] = TOPFACE;
      directions[2][2] = -2;
      directions[2][3] = -4;

      directions[3][0] = -5;
      directions[3][1] = TOPFACE;
      directions[3][2] = FRONTFACE;
      directions[3][3] = -3;

      directions[4][0] = RIGHTFACE;
      directions[4][1] = -4;
      directions[4][2] = FRONTFACE;
      directions[4][3] = -6;

      directions[5][0] = RIGHTFACE;
      directions[5][1] = -1;
      directions[5][2] = -5;
      directions[5][3] = BOTTOMFACE;
    }
  else
    {
      /* This is just to avoid compiler warnings */

      Die(FUNCTION_NAME, "Shouldn't be here");

      pts[0][0] = -1;
      pts[0][1] = -1;
      pts[0][2] = -1;
      pts[0][3] = -1;

      pts[1][0] = -1;
      pts[1][1] = -1;
      pts[1][2] = -1;
      pts[1][3] = -1;

      pts[2][0] = -1;
      pts[2][1] = -1;
      pts[2][2] = -1;
      pts[2][3] = -1;

      pts[3][0] = -1;
      pts[3][1] = -1;
      pts[3][2] = -1;
      pts[3][3] = -1;

      pts[4][0] = -1;
      pts[4][1] = -1;
      pts[4][2] = -1;
      pts[4][3] = -1;

      pts[5][0] = -1;
      pts[5][1] = -1;
      pts[5][2] = -1;
      pts[5][3] = -1;

      directions[0][0] = -1;
      directions[0][1] = -1;
      directions[0][2] = -1;
      directions[0][3] = -1;

      directions[1][0] = -1;
      directions[1][1] = -1;
      directions[1][2] = -1;
      directions[1][3] = -1;

      directions[2][0] = -1;
      directions[2][1] = -1;
      directions[2][2] = -1;
      directions[2][3] = -1;

      directions[3][0] = -1;
      directions[3][1] = -1;
      directions[3][2] = -1;
      directions[3][3] = -1;

      directions[4][0] = -1;
      directions[4][1] = -1;
      directions[4][2] = -1;
      directions[4][3] = -1;

      directions[5][0] = -1;
      directions[5][1] = -1;
      directions[5][2] = -1;
      directions[5][3] = -1;
    }

  /* Reset first flags */

  for (i = 0; i < 6; i++)
    firsts[i] = YES;

  for (i = 0; i < 6; i++)
    for (j = 0; j < 6; j++)
      intfirsts[i][j] = YES;	  

  /* Create cells */

  for (i = 0; i < nc; i++)
    {
      /* Set bounding box */
    
      p0 = pts[i][0];
      p1 = pts[i][1];
      p2 = pts[i][2];
      p3 = pts[i][3];

      p0 = pointlist + 3*p0;
      p1 = pointlist + 3*p1;
      p2 = pointlist + 3*p2;
      p3 = pointlist + 3*p3;

      tmpface[0] = p0;
      tmpface[1] = p1;
      tmpface[2] = p2;
      tmpface[3] = p3;

      TetPutBoundingBox(ncgns, pointlist, tmpface);

      for (j = 0; j < 4; j++)
	{
	  /* Get points */

	  switch (j)
	    {
	    case 0:
	      p0 = pts[i][0];
	      p1 = pts[i][1];
	      p2 = pts[i][2];

	      break;

	    case 1:
	      p0 = pts[i][1];
	      p1 = pts[i][3];
	      p2 = pts[i][2];

	      break;

	    case 2:
	      p0 = pts[i][0];
	      p1 = pts[i][3];
	      p2 = pts[i][1];

	      break;

	    case 3:
	      p0 = pts[i][0];
	      p1 = pts[i][2];
	      p2 = pts[i][3];

	      break;
	    }

	  /* Get pointer to beginning of point */

	  p0 = pointlist + 3*p0;
	  p1 = pointlist + 3*p1;
	  p2 = pointlist + 3*p2;

	  /* Get first flag and unset it for future */

	  if (directions[i][j] >= 0)
	    {
	      /* Out of cell */

	      idx = directions[i][j];

	      /* out of cell */

	      first = firsts[idx];

	      /* Next one won't be first anymore */

	      firsts[idx] = NO;
	    }
	  else
	    {
	      /* Internal face from this subcell */

	      idx = -directions[i][j] - 1;

	      /* Internal face */

	      first = intfirsts[i][idx];
		      
	      /* Next one won't be first anymore */

	      intfirsts[i][idx] = NO;
	      intfirsts[idx][i] = NO;
	    }

	  HexNewTetFace(ifc, cgns, ncgns, hexnbrs, hexfaces, hexsides, 
			directions[i][j], first, p0, p1, p2, ownrlist2, 
			nbrlist2, sdone, j, i);

	}

      /* Next subcell */

      ncgns = NextItem(ncgns);
    }
}

/*****************************************************************************/
