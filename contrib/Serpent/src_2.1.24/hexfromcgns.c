/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexfromcgns.c                                  */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Takes a IFC_HEX_MSH cgns cell and creates an 8 point hex     */
/*              out of it                                                    */
/*                                                                           */
/* Comments:                                                                 */
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

#define FUNCTION_NAME "HexFromCGNS:"

/*****************************************************************************/

void HexFromCGNS(long ifc, long cgns, long *V, 
		 long *hexfaces, long *hexsides, long (*hexnbrs)[6],
		 long (*initFaces)[4], long *ownrlist2, long *nbrlist2)
{
  long loc2, loc3;
  long ptr, ptr1, pointlist, surflist, surf;
  long ownrlist, nbrlist;
  long nf, np, i, j, n; 
  long side, oppose, id;
  long topfound, sidefound;
  long p0, p1;
  long bot[4], top[4], left[4];

  /* Get pointer to parents pointlist */

  pointlist = (long)RDB[ifc + IFC_PTR_POINT_LIST_PRNTS];;
  CheckPointer(FUNCTION_NAME, "(pointlist)", DATA_ARRAY, pointlist);

  /* Get pointer to parents surfacelist */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);
 
  /* Get pointer to parents owner list */

  ownrlist = (long)RDB[ifc + IFC_PTR_OWNR_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(ownrlist)", DATA_ARRAY, ownrlist);

  /* Get pointer to parents neighbour list */

  nbrlist  = (long)RDB[ifc + IFC_PTR_NBR_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(nbrlist)", DATA_ARRAY, nbrlist);

  /* Get number of faces */

  nf = (long)RDB[cgns + IFC_TET_MSH_NF];

  /* Get pointer to face list */

  ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(PTR)", DATA_ARRAY, ptr);

  /* Get pointer to side list */

  ptr1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
  CheckPointer(FUNCTION_NAME, "(PTR1)", DATA_ARRAY, ptr1);

  /* Make lists of hex faces, sides and neighbours */

  for (j = 0; j < 6; j++)
    {
      hexfaces[j] = (long)RDB[ptr + j];
      hexsides[j] = (long)RDB[ptr1 + j];

      if (hexsides[j] == -1)
	{
	  /* Owner */
	  hexnbrs[0][j] = (long)RDB[nbrlist + hexfaces[j]];

	  /* Get second neighbour (if any) */
	  hexnbrs[1][j] = nbrlist2[hexfaces[j]];

	}
      else
	{
	  /* Neighbour */
	  hexnbrs[0][j] = (long)RDB[ownrlist + hexfaces[j]];

	  /* Get second owner (if any) */
	  hexnbrs[1][j] = ownrlist2[hexfaces[j]];

	}

    }

  /* Get index of first face */

  n = (long)RDB[ptr + 0];

  /* Get pointer to surface */

  surf = ListPtr(surflist, n);
  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get side of first face */

  side = (long)RDB[ptr1 + 0];

  /* Get number of points on the face */

  np = (long)RDB[surf + SURFACE_N_PARAMS];

  /* Get pointer to surface parameters (point list) */
	      
  loc2 = (long)RDB[surf + SURFACE_PTR_PARAMS];

  /* Loop over points */

  for (i = 0; i < np; i++)
    {
      /* Get point */
      loc3 = (long)RDB[loc2 + i];

      /* Calculate point id */

      id = (long)((loc3 - pointlist)/3);

      /* Store first face        */
      /* This will be the bottom */

      if (side == -1)
	{
	  /* Owner */

	  /* Store index (change direction) */

	  bot[3-i] = id;

	}
      else if (side == 1)
	{
	  /* Store index (do not change) */

	  bot[i] = id;

	}
      else
	Die(FUNCTION_NAME, "WTF?");

      initFaces[0][i] = id;

    }

  /* We still need to find an opposing face (top) */
  /* And one of the others (side/left)            */
  /* Reset found flags                            */

  topfound = 0;
  sidefound = 0;

  /* Loop over other faces to find the opposing side */

  for (i = 1; i < nf; i++)
    {
      /* Get index of face */

      n = (long)RDB[ptr + i];

      /* Get side of face */

      side = (long)RDB[ptr1 + i];

      /* Get pointer to surface */

      surf = ListPtr(surflist, n);
      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

      /* Get number of points on the face */

      np = (long)RDB[surf + SURFACE_N_PARAMS];

      /* Get pointer to surface parameters (point list) */
	      
      loc2 = (long)RDB[surf + SURFACE_PTR_PARAMS];

      oppose = YES;

      /* Loop over points */

      for (j = 0; j < np; j++)
	{
	  /* Get point */
	  loc3 = (long)RDB[loc2 + j];

	  id = (long)((loc3 - pointlist)/3);

	  /* Check if this is opposing side to the one we stored */
	  /* They do not share any points */

	  if (PolyPInF(id, bot, 4) == YES)
	    oppose = NO;

	  /* Store second face (Top) */

	  if (side == -1)
	    {
	      /* Owner */

	      /* Store index (Do not rotate) */

	      if (topfound == 0)
		top[j] = id;
	      else
		left[j] = id;

	    }
	  else if (side == 1)
	    {
	      /* Store index (Rotate) */

	      if (topfound == 0)
		top[3-j] = id;
	      else
		left[3-j] = id;

	    }
	  else
	    Die(FUNCTION_NAME, "WTF?");

	  /* Store face points to initFaces */

	  initFaces[i][j] = id;

	}

      /* Check if this was opposing face */

      if (oppose == YES)
	{
	  topfound = 1;
	}
      else 
	{
	  /* One of the sides */
	  /* Store sides */

	  if (topfound == 0)
	    for (j = 0; j < 4; j++)		
	      left[j] = top[j];       	

	  sidefound = 1;

	}

    }

  /* Find pair of points from side that are connected in top and bottom */

  for (i = 0; i < 3; i++)
    if (((PolyPInF(left[i], bot, 4)) && (PolyPInF(left[i+1], top, 4))) ||
	((PolyPInF(left[i], top, 4)) && (PolyPInF(left[i+1], bot, 4))))
      break;

  if (i == 3)
    Die(FUNCTION_NAME, "Bad side");

  if (PolyPInF(left[i], bot, 4))
    {
      p0 = left[i];
      p1 = left[i+1];
    }
  else
    {
      p0 = left[i+1];
      p1 = left[i];
    }

  /* Now p0 is in bot and p1 in top */
  /* They connect these two opposite faces */
  /* Let's rotate top until p0 and p1 match */

  for (i = 0; i < 4; i++)
    if (bot[i] == p0)
      break;

  for (j = 0; j < 4; j++)
    if (top[j] == p1)
      break;

  /* Top will have to be rotated i-j steps */

  HexRotateFace(top, i-j);

  /* Create the hexahedron */
      
  for (i = 0; i < 4; i++)
    {
      V[i] = bot[i];
      V[4+i] = top[i];
    }
}

/*****************************************************************************/
