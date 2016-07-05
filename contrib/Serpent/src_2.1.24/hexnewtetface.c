/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : hexnewtetface.c                                */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/06/17 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Creates a new tet-face when splitting hexahedrons            */
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

#define FUNCTION_NAME "FixHexMesh:"

/*****************************************************************************/

void HexNewTetFace(long ifc, long cgns, long ncgns, long (*hexnbrs)[6], 
		   long hexfaces[6], long hexsides[6], long direction, 
		   long first, long p0, long p1, long p2,  
		   long *ownrlist2, long *nbrlist2, 
		   long *sdone, long fdone, long cdone)
{
  long i, j;
  long nbr0, nbr1, nfidx, fidx, other, side;
  long idxthis, idxother, ofaces, osides;
  long nfacelist, loc1, loc2, ownrlist, nbrlist;
  long nownrs, nnbrs, nfaces, nsides, nsurf, ptr;

  /* Get pointer to new owners list */

  nownrs = (long)RDB[ifc + IFC_PTR_OWNR_LIST];
  CheckPointer(FUNCTION_NAME, "(nownrs)", DATA_ARRAY, nownrs);

  /* Get pointer to new neighbours list */

  nnbrs = (long)RDB[ifc + IFC_PTR_NBR_LIST];
  CheckPointer(FUNCTION_NAME, "(nnbrs)", DATA_ARRAY, nnbrs);

  /* Get pointer to parents owner list */

  ownrlist = (long)RDB[ifc + IFC_PTR_OWNR_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(ownrlist)", DATA_ARRAY, ownrlist);

  /* Get pointer to parents neighbour list */

  nbrlist  = (long)RDB[ifc + IFC_PTR_NBR_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(nbrlist)", DATA_ARRAY, nbrlist);

  /* Get pointer to new face list*/

  nfacelist = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(nfacelist)", DATA_ARRAY, nfacelist);

  /* Get pointer to face list of subcell */

  nfaces = (long)RDB[ncgns + IFC_TET_MSH_PTR_FACES];
  CheckPointer(FUNCTION_NAME, "(nfaces)", DATA_ARRAY, nfaces);

  /* Get pointer to side list */

  nsides = (long)RDB[ncgns + IFC_TET_MSH_PTR_SIDES];
  CheckPointer(FUNCTION_NAME, "(nsides)", DATA_ARRAY, nsides);

  /* Avoid compiler warning */

  other = -1;

  /* If direction is non-negative it is out of cell */

  if (direction >= 0)
    {

      /* Get side */

      side = hexsides[direction];

      /* Check if cell behind first new face has been divided */

      nbr0 = hexnbrs[0][direction];
      nbr1 = hexnbrs[1][direction];

      if((long)RDB[nbr0 + IFC_TET_MSH_IDX] != (long)RDB[nbr1 + IFC_TET_MSH_IDX])
	{
	  /* Has been divided */

	  /* Get the other cell */
	  /* We don't know which one (nbr0 or nbr1) is the corresponding cell */
	  /* on the other side of this face     */
	  /* We'll need to check the facepoints */

	  for (i = 0; i < 2; i++)
	    {
	      /* Handle both possible neighbours */

	      if (i == 0)
		other = nbr0;
	      else
		other = nbr1;

	      /* Get pointers to other's faces and sides */

	      ofaces = (long)RDB[other + IFC_TET_MSH_PTR_FACES];
	      osides = (long)RDB[other + IFC_TET_MSH_PTR_SIDES];

	      /* Loop over faces of other */

	      for (j = 0; j < 4; j++)
		{
		  
		  /* Get face index from cell */

		  nfidx = (long)RDB[ofaces + j];			      

		  /* Get face from face list */

		  loc1 = ListPtr(nfacelist, nfidx);

		  /* Get pointer to surface parameters (point list) */
	      
		  loc2 = (long)RDB[loc1 + SURFACE_PTR_PARAMS];

		  /* Check if the face points are equal to the face points */
		  /* of the face to be found (side = 1) or the same but    */
		  /* inverted (side = -1) */
		  		 
		  if (side == -1)
		    {
		      if ((((long)RDB[loc2 + 0] == p2) && 
			   ((long)RDB[loc2 + 1] == p1) && 
			   ((long)RDB[loc2 + 2] == p0)) ||
			  (((long)RDB[loc2 + 0] == p0) && 
			   ((long)RDB[loc2 + 1] == p2) && 
			   ((long)RDB[loc2 + 2] == p1)) ||
			  (((long)RDB[loc2 + 0] == p1) && 
			   ((long)RDB[loc2 + 1] == p0) && 
			   ((long)RDB[loc2 + 2] == p2)))
			break;
		      
		    }
		  else
		    {
		      if ((((long)RDB[loc2 + 0] == p0) && 
			   ((long)RDB[loc2 + 1] == p1) && 
			   ((long)RDB[loc2 + 2] == p2)) ||
			  (((long)RDB[loc2 + 0] == p2) && 
			   ((long)RDB[loc2 + 1] == p0) && 
			   ((long)RDB[loc2 + 2] == p1)) ||
			  (((long)RDB[loc2 + 0] == p1) && 
			   ((long)RDB[loc2 + 1] == p2) && 
			   ((long)RDB[loc2 + 2] == p0)))
			break;

		    }
		      
		}

	      /* Check if found */

	      if (j < 4)
		break;
	    }

	  /* Should be found */

	  if (i == 2)
	    Die(FUNCTION_NAME, "Could not find other");

	  /* New face index is now nfidx, other cell is now in other */

	  /* Store face */

	  WDB[nfaces + fdone] = (double)nfidx;

	  /* Store side */

	  WDB[nsides + fdone] = (double)side;

	  /* Update neighbors */

	  if(side == -1)
	    {      
	      WDB[nownrs + nfidx] = (double)ncgns;

	      /* Check that other is already neighbour */

	      if ((long)RDB[nnbrs  + nfidx] != other)
		Die(FUNCTION_NAME, "Other is not other");
	    }
	  else
	    {
	      /* Check that other is already owner */

	      if ((long)RDB[nownrs  + nfidx] != other)
		Die(FUNCTION_NAME, "Other is not owner");

	      /* Store this as neighbour */

	      WDB[nnbrs  + nfidx] = (double)ncgns;
	    }

	  /* Get face index */

	  fidx = hexfaces[direction];

	  /* Update owner and neighbour to the parent lists */

	  if(side == -1)
	    {
	      if (i == 0)
		{
		  /* If this corresponded to the first subcell on this face */
		  /* Update to first ownerlist */
		  WDB[ownrlist + fidx] = (double)ncgns;
		}
	      else
		{
		  /* If this is the second subcell on this face */
		  /* Update to second ownerlist */
		  ownrlist2[fidx] = (double)ncgns;
		}
	    }
	  else
	    {

	      if (i == 0)
		{
		  /* If this is the first subcell on this face */
		  /* Update to first neighbourlist */

		  WDB[nbrlist + fidx] = (double)ncgns;
		}
	      else
		{
		  /* If this is the second subcell on this face */
		  /* Update to second neighbourlist */

		  nbrlist2[fidx] = (double)ncgns;

		}
	    }

	}
      else
	{
	  /* Has not been divided*/

	  /* Get new face */

	  nsurf = (long)RDB[ifc + IFC_PTR_SURF_LIST];

	  /* Create surface and put pointer to temporary array */
	  
	  nsurf = ListPtr(nsurf, *sdone);

	  *sdone = *sdone + 1;

	  /* Get the other cell */

	  if (first)
	    other = nbr0;
	  else
	    other = nbr1;

	  /* Get pointer to surface parameters */

	  ptr = (long)RDB[nsurf + SURFACE_PTR_PARAMS];

	  /* Put the points for the parameters        */
	  /* Invert the rotation if this is neighbour */

	  if (side == -1)
	    {
	      WDB[ptr + 0] = (double)p2;
	      WDB[ptr + 1] = (double)p1;
	      WDB[ptr + 2] = (double)p0;
	    }
	  else
	    {
	      WDB[ptr + 0] = (double)p0;
	      WDB[ptr + 1] = (double)p1;
	      WDB[ptr + 2] = (double)p2;
	    }

	  /* This surface's index can be found out from list size */
  
	  nfidx = *sdone - 1;

	  /* Store face */

	  WDB[nfaces + fdone] = (double)nfidx;

	  /* Store side */

	  WDB[nsides + fdone] = (double)side;

	  /* Update neighbors */
	  if(side == -1)
	    {      
	      WDB[nownrs + nfidx] = (double)ncgns;
	      WDB[nnbrs  + nfidx] = (double)other;
	    }
	  else
	    {
	      WDB[nownrs + nfidx] = (double)other;
	      WDB[nnbrs  + nfidx] = (double)ncgns;
	    }

	  /* Put this cell as the first neighbor of the initial (undivided) face */

	  /* Get face index */

	  fidx = hexfaces[direction];

	  if(side == -1)
	    {
	      if (first)
		{
		  /* If this is the first subcell on this face */
		  /* Update to first ownerlist */
		  WDB[ownrlist + fidx] = (double)ncgns;
		}
	      else
		{
		  /* If this is the second subcell on this face */
		  /* Update to second ownerlist */
		  ownrlist2[fidx] = (double)ncgns;
		}
	    }
	  else
	    {

	      if (first)
		{
		  /* If this is the first subcell on this face */
		  /* Update to first neighbourlist */

		  WDB[nbrlist + fidx] = (double)ncgns;
		}
	      else
		{
		  /* If this is the second subcell on this face */
		  /* Update to second neighbourlist */

		  nbrlist2[fidx] = (double)ncgns;

		}
	    }

	}



    }
  else
    {

      /* This is an internal face between the new subcells */

      if (first)
	{
	  /* Create face only on second pass */

	  /* Put a placeholder                  */
	  /* If this is internal face to cell 4 */
	  /* The placeholder will be -4         */

	  WDB[nfaces + fdone] = (double)direction;

	  /* Side will be put when the face is created */
	  /* in the neighbouring subcell               */

	  return;

	}
      else
	{

	  /* This is the second pass */

	  /* Get number of this subcell (1-6) */

	  idxthis = cdone + 1;

	  /* Get number of the other subcell (1-6) */	  

	  idxother = -direction;

	  /* Get pointer to the other cell by going backwards from this one */

	  other = ncgns;
	  for (i = 0; i < idxthis - idxother; i++)
	    other = PrevItem(other);	  

	  /* Find this subcell from the face list of the other cell */
	  /* Placeholder was stored on the first pass               */

	  /* Get pointer to other's face and side lists */

	  ofaces = (long)RDB[other + IFC_TET_MSH_PTR_FACES];
	  osides = (long)RDB[other + IFC_TET_MSH_PTR_SIDES];

	  for (i = 0; i < 4; i++)
	    {
	      /* Check if i:th face is placeholder for this */
	      if ((long)RDB[ofaces + i] == -idxthis)
		break;

	    }

	  /* Should be found */

	  if (i == 4)
	    Die(FUNCTION_NAME, "Internal face placeholder not found");

	  /* Get new face surface from surface list */

	  nsurf = (long)RDB[ifc + IFC_PTR_SURF_LIST];
	  nsurf = ListPtr(nsurf, *sdone);

	  /* Increment number of used surfaces */

	  *sdone = *sdone + 1;

	  /* Get pointer to surface parameters */

	  ptr = (long)RDB[nsurf + SURFACE_PTR_PARAMS];

	  /* Put the points for the parameters        */
	  /* This is always the owner                 */

	  WDB[ptr + 0] = (double)p2;
	  WDB[ptr + 1] = (double)p1;
	  WDB[ptr + 2] = (double)p0;

	  /* Get this surfaces index */
  
	  nfidx = *sdone - 1;

	  /* Store face to this */
	  WDB[nfaces + fdone] = (double)nfidx;

	  /* Store side (this is owner) */
	  WDB[nsides + fdone] = (double)(-1);

	  /* Store face to other */
	  WDB[ofaces + i] = (double)nfidx;

	  /* Store side to other (always neighbour) */
	  WDB[osides + i] = (double)(1);

	  /* Update owner and neighbour of face */

	  WDB[nownrs + nfidx] = (double)ncgns;
	  WDB[nnbrs  + nfidx] = (double)other;

	  
	}
    }

}


