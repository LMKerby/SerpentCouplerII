#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dividepolyhedface.c                            */
/*                                                                           */
/* Created:       2014/02/24 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Creates subcells on a polyhedral cell-face by adding points  */
/*              to face centerpoints as well as to the centerpoint of        */
/*              face centerpoints ("cell centerpoint")                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DividePolyhedFace:"

/*****************************************************************************/

void DividePolyhedFace(long ifc, long cidx, long fidx, long side, long np, 
		       long *newcells, long *subcells, long nss,  long *sdone)
{
  long  p0, p1, p2, p3;
  long sparams, side2, other0;  
  long cgns, cgns1, surflist, nfidx;
  long i, j, k, l, ncgns, nsurf, prev, ptr, ptr1;
  long idx, nf, nsides, nfaces;
  long new, *tempcells;
  long surf, other, *nbrcells;
  long nownrs, nnbrs;

  /* Get pointer to parent list */

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH_PRNTS];
  CheckPointer(FUNCTION_NAME, "(cgnslist)", DATA_ARRAY, cgns);

  /* Get pointer to cell to be divided */

  cgns = ListPtr(cgns, cidx);
  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

  /* Get face surface */
  /* Get pointer to interface surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST_PRNTS];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  /* Get pointer to surface of this face */

  surf = ListPtr(surflist, fidx);
  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  /* Get pointer of surface parameters (points) of old face */

  sparams = (long)RDB[surf + SURFACE_PTR_PARAMS];

  /* Allocate memory for new cells */

  tempcells = (long *)Mem(MEM_ALLOC, np, sizeof(long));
  nbrcells = (long *)Mem(MEM_ALLOC, np, sizeof(long));

  /* Get cell on the other side of this face */

  /* This cell owns this surface */
  if(side == -1)
    {
      /* This cell owns this surface */

      ptr = (long)WDB[ifc + IFC_PTR_NBR_LIST_PRNTS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      other = (long)RDB[ptr + fidx];
    }
  else
    {
      /* This cell is this surfaces neighbour */

      ptr = (long)WDB[ifc + IFC_PTR_OWNR_LIST_PRNTS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      other = (long)RDB[ptr + fidx];
      CheckPointer(FUNCTION_NAME, "(other)", DATA_ARRAY, other);
    }

  /* No cell on the other side */

  if (other < VALID_PTR)
    {
      if(side == 1)
	Die(FUNCTION_NAME, "This is neighbour but no owner found?!?");

      new = 1;
    }
  else
    {

      /* Try to find this face in 'other' cells face list */

      /* Get number of faces for the other cell */

      nf = (long)RDB[other + IFC_TET_MSH_NF];

      /* Get pointer to others faces */

      ptr = (long)RDB[other + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Loop over 'others' faces */
      
      for(i = 0; i < nf; i++)
	{
	  
	  /* Check for negative index */
	  /* (child face) */
      
	  if ((long)RDB[ptr + i] < 0)
	    break;

	}
            
      /* if we found a negative index, 'other' has been divided */

      if (i >= nf)
	{
	  new = 1;
	}
      else
	{
	  /* The cell on the other side has been divided */

	  /* Get cells facing this face in the neighbor cell */
	  /* They will be neighbors for this face's subcells */
	  /* Store first other for checking */

	  other0 = other;

	  for(j = 0; j < np; j++)
	    {
	      /* Store cell pointer */

	      nbrcells[j] = other;

	      /* Get pointer to others faces */

	      ptr = (long)RDB[other + IFC_TET_MSH_PTR_FACES];

	      /* Get index of second face of 'other' */
	      /* second face == forward around face  */
	      /* will be negative */

	      idx = (long)RDB[ptr + 1];

	      if(idx > 0)
		Die(FUNCTION_NAME, "Positive face index even though divided");

	      /* Get pointer to others sides */

	      ptr = (long)RDB[other + IFC_TET_MSH_PTR_SIDES];

	      /* Get side of the second face */

	      side2 = (long)RDB[ptr + 1];

	      /* Get neighbour of other through second face */

	      if(side2 == -1)
		ptr = (long)RDB[ifc + IFC_PTR_NBR_LIST];
	      else
		ptr = (long)RDB[ifc + IFC_PTR_OWNR_LIST];

	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Get neighbour cell pointer */
	      /* On last iteration, this should be the original other */

	      other = (long)RDB[ptr - idx];

	    }
	  
	  /* Check that we came back to other0 */

	  if(other != other0)
	    Die(FUNCTION_NAME, "Did not conduct a loop %ld", other);

	  new = 0;
	}
    }

  /* Get pointer to new owner list */

  nownrs = (long)RDB[ifc + IFC_PTR_OWNR_LIST];

  /* Get pointer to new neighbour list */

  nnbrs = (long)RDB[ifc + IFC_PTR_NBR_LIST];

  /* Create np cells and their face and side lists */

  for(i = 0; i < np; i++)
    {
      /* Create new tet cell */

      ncgns = NewItem(ifc + IFC_PTR_TET_MSH, IFC_TET_MSH_LIST_BLOCK_SIZE);
      CheckPointer(FUNCTION_NAME, "(ncgns)", DATA_ARRAY, ncgns);

      /* Get tet mesh index */

      if((prev = PrevItem(ncgns)) > VALID_PTR)
	idx = (long)RDB[prev + IFC_TET_MSH_IDX] + 1;
      else
	idx = (long)RDB[ifc + IFC_NC_PRNTS];

      /* Put index */

      WDB[ncgns + IFC_TET_MSH_IDX] = (double)idx;

      /* Reset cell boundaries */

      WDB[ncgns + IFC_TET_MSH_XMIN] = INFTY;
      WDB[ncgns + IFC_TET_MSH_XMAX] = -INFTY;
      WDB[ncgns + IFC_TET_MSH_YMIN] = INFTY;
      WDB[ncgns + IFC_TET_MSH_YMAX] = -INFTY;
      WDB[ncgns + IFC_TET_MSH_ZMIN] = INFTY;
      WDB[ncgns + IFC_TET_MSH_ZMAX] = -INFTY;

      /* Put number of faces */

      WDB[ncgns + IFC_TET_MSH_NF] = 4.0;

      /* Allocate memory for face list */

      nfaces = ReallocMem(DATA_ARRAY, 4.0);

      /* Put pointer to face list */

      WDB[ncgns + IFC_TET_MSH_PTR_FACES] = (double)nfaces;

      /* Allocate memory for side list */

      nsides = ReallocMem(DATA_ARRAY, 4.0);

      /* Put pointer to side list */

      WDB[ncgns + IFC_TET_MSH_PTR_SIDES] = (double)nsides;
	  
      /* Store the created cell to tempcells */

      tempcells[i] = ncgns;

    }


  /* nbrcells is now an empty table or contains the new cells on the other */
  /* side of surf its length is np */

  for(i = 0; i < np; i++)
    {
      /* Create new tet cell if the initial one cannot be recycled */
      /* Do not recycle old ones anymore to allow for parent->child updates */

      /* Get i:th tempcell */

      ncgns = tempcells[i];

      /* Get pointer to face list */

      nfaces = (long)RDB[ncgns + IFC_TET_MSH_PTR_FACES];

      /* Get pointer to side list */

      nsides = (long)RDB[ncgns + IFC_TET_MSH_PTR_SIDES];

      /* Get next two points from surface perimeter */

      p0 = (long)RDB[sparams + i];

      if(i < np - 1)
	p1 = (long)RDB[sparams + i + 1];
      else
	p1 = (long)RDB[sparams + 0];

      /* Get pointer to face centerpoint list */

      ptr = (long)RDB[ifc + IFC_PTR_PRNT_FACE_CP_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Create face centerpoint pointer */

      p2 = ptr + fidx*3;
 
      /* Get pointer to cell centerpoint list */

      ptr = (long)RDB[ifc + IFC_PTR_PRNT_CELL_CP_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Create cell centerpoint pointer */

      p3 = ptr + cidx*3;

      /* Put bounding box */

      if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[p0 + 0])
	WDB[ncgns + IFC_TET_MSH_XMIN]  = RDB[p0 + 0];
      if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[p0 + 0])
	WDB[ncgns + IFC_TET_MSH_XMAX]  = RDB[p0 + 0];      
      if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[p0 + 1])
	WDB[ncgns + IFC_TET_MSH_YMIN]  = RDB[p0 + 1];
      if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[p0 + 1])
	WDB[ncgns + IFC_TET_MSH_YMAX]  = RDB[p0 + 1];      
      if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[p0 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMIN]  = RDB[p0 + 2];
      if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[p0 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMAX]  = RDB[p0 + 2];      

      if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[p1 + 0])
	WDB[ncgns + IFC_TET_MSH_XMIN]  = RDB[p1 + 0];
      if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[p1 + 0])
	WDB[ncgns + IFC_TET_MSH_XMAX]  = RDB[p1 + 0];      
      if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[p1 + 1])
	WDB[ncgns + IFC_TET_MSH_YMIN]  = RDB[p1 + 1];
      if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[p1 + 1])
	WDB[ncgns + IFC_TET_MSH_YMAX]  = RDB[p1 + 1];      
      if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[p1 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMIN]  = RDB[p1 + 2];
      if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[p1 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMAX]  = RDB[p1 + 2];      

      if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[p2 + 0])
	WDB[ncgns + IFC_TET_MSH_XMIN]  = RDB[p2 + 0];
      if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[p2 + 0])
	WDB[ncgns + IFC_TET_MSH_XMAX]  = RDB[p2 + 0];      
      if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[p2 + 1])
	WDB[ncgns + IFC_TET_MSH_YMIN]  = RDB[p2 + 1];
      if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[p2 + 1])
	WDB[ncgns + IFC_TET_MSH_YMAX]  = RDB[p2 + 1];      
      if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[p2 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMIN]  = RDB[p2 + 2];
      if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[p2 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMAX]  = RDB[p2 + 2];      

      if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[p3 + 0])
	WDB[ncgns + IFC_TET_MSH_XMIN]   = RDB[p3 + 0];
      if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[p3 + 0])
	WDB[ncgns + IFC_TET_MSH_XMAX]   = RDB[p3 + 0];      
      if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[p3 + 1])
	WDB[ncgns + IFC_TET_MSH_YMIN]   = RDB[p3 + 1];
      if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[p3 + 1])
	WDB[ncgns + IFC_TET_MSH_YMAX]   = RDB[p3 + 1];      
      if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[p3 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMIN]   = RDB[p3 + 2];
      if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[p3 + 2])
	WDB[ncgns + IFC_TET_MSH_ZMAX]   = RDB[p3 + 2];      
      /*
      printf("Bounding box of new cell (%E %E) (%E %E) (%E %E)\n",
	     RDB[ncgns + IFC_TET_MSH_XMIN], RDB[ncgns + IFC_TET_MSH_XMAX],
	     RDB[ncgns + IFC_TET_MSH_YMIN], RDB[ncgns + IFC_TET_MSH_YMAX],
	     RDB[ncgns + IFC_TET_MSH_ZMIN], RDB[ncgns + IFC_TET_MSH_ZMAX]);

      */
      /* Create three surfaces */
      /* Fourth surface will be one of these */

      for (k = 0; k < 3; k++)
	{
	  /* Put three points to the parameters */
	  ptr1 = (long)WDB[ifc + IFC_PTR_PRNT_CELL_CP_LIST];

	  if(k == 0)
	    {
	      /* out of cell side */

	      /* If the other side has not been divided */
	      /* these surfaces have to be created      */

	      if(new == 1)
		{

		  nsurf = (long)RDB[ifc + IFC_PTR_SURF_LIST];

		  /* Create surface and put pointer to temporary array */

		  nsurf = ListPtr(nsurf, *sdone);

		  *sdone = *sdone + 1;

		  /* Get pointer to surface paramters */

		  ptr = (long)RDB[nsurf + SURFACE_PTR_PARAMS];

		  /* Put the points for the parameters */

		  WDB[ptr + 0] = (double)p0;
		  WDB[ptr + 1] = (double)p1;
		  WDB[ptr + 2] = (double)p2;

		  /* Put this surface's (negative) index to this cells surface list */
		  /* First face is out of the cell */

		  /* This surface's index can be found out from list size */

		  nfidx = *sdone - 1;

		  /* Put face */

		  WDB[nfaces + 0] = -(double)nfidx;

		  /* Put side */

		  WDB[nsides + 0] = (double)side;

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

		  ptr1 = (long)RDB[ifc + IFC_PTR_PRNT_CELL_CP_LIST];

		  /* Put this cell as the neighbor of the initial (undivided) face */

		  if(side == -1)
		    ptr = (long)RDB[ifc + IFC_PTR_OWNR_LIST_PRNTS];
		  else
		    ptr = (long)RDB[ifc + IFC_PTR_NBR_LIST_PRNTS];

		  WDB[ptr + fidx] = (double)ncgns;




		}
	      else
		{
		  /* The other side of this face has already been divided */
		  /* This means that the surfaces have already been created*/

		  /* Find correct surface */

		  for(l = 0; l < np; l++)
		    {
		      /* Get pointer to list of neighbors */

		      cgns1 = nbrcells[l];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

		      /* Get pointer to face list */

		      ptr = (long)RDB[cgns1 + IFC_TET_MSH_PTR_FACES];
		      CheckPointer(FUNCTION_NAME, "(ptrfl)", DATA_ARRAY, ptr);

		      /* Get the first face (this is still negative) */

		      idx = -(long)RDB[ptr + 0];

		      /* Get side of other cell */

		      ptr = (long)RDB[cgns1 + IFC_TET_MSH_PTR_SIDES];
		      CheckPointer(FUNCTION_NAME, "(ptrfl)", DATA_ARRAY, ptr);

		      side2 = (long)RDB[ptr + 0];

		      if(idx < 0)
			Die(FUNCTION_NAME, "wrong sing of face index");

		      /* Check if it's points are 0, 1, 2 */

		      surf = (long)RDB[ifc + IFC_PTR_SURF_LIST];

		      /* Get pointer to surface */

		      surf = ListPtr(surf, idx);
		      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
		      /* Get pointer to surface parameters */

		      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];


		      if(((long)RDB[ptr + 0] == p0) &&
			 ((long)RDB[ptr + 1] == p1) &&
			 ((long)RDB[ptr + 2] == p2))
			break;			

		    }

		  /* Check if found (should be) */

		  if (l >= np)
		    Die(FUNCTION_NAME, "Opposing surface not found");

		  /* Other cell is now in cgns1 */
		  /* Face is now in surf, face index in idx */

		  WDB[nfaces + 0] = -(double)idx;

		  /* Put side */

		  WDB[nsides + 0] = -(double)side2;
		  
		  /* Put cell pointers to owner and neighbor lists */

		  if(-side2 == -1)
		    {
		      WDB[nownrs + idx] = (double)ncgns;
		      WDB[nnbrs  + idx] = (double)cgns1;
		    }
		  else
		    {
		      WDB[nownrs + idx] = (double)cgns1;
		      WDB[nnbrs  + idx] = (double)ncgns;
		    }

		}
	    }
	  
	  /* forward on this face side */
	  /* This should be the fourth face for the next cell */
	  
	  else if (k == 1)
	    {

	      nsurf = (long)RDB[ifc + IFC_PTR_SURF_LIST];
	      
	      /* Get new surface */

	      nsurf = ListPtr(nsurf, *sdone);

	      *sdone = *sdone + 1;

	      /* Get pointer to surface parameters */

	      ptr = (long)RDB[nsurf + SURFACE_PTR_PARAMS];

	      /* Put the points for the parameters */

	      WDB[ptr + 0] = (double)p1;
	      WDB[ptr + 1] = (double)p3;
	      WDB[ptr + 2] = (double)p2;

	      /* This surface's index can be found out from list size */

	      nfidx = *sdone - 1;

	      /* Put face */

	      WDB[nfaces + 1] = -(double)nfidx;

	      /* Put side */

	      WDB[nsides + 1] = (double)side;

	      /* Get face and side list for next cell */
	      
	      if(i < np - 1)
		cgns1 = tempcells[i+1];
	      else
		cgns1 = tempcells[0];

	      /* Put cell pointers to owner and neighbor lists */

	      if(side == -1)
		{
		  WDB[nownrs + nfidx] = (double)ncgns;
		  WDB[nnbrs  + nfidx] = (double)cgns1;
		}
	      else
		{
		  WDB[nownrs + nfidx] = (double)cgns1;
		  WDB[nnbrs  + nfidx] = (double)ncgns;
		}

	      /* Get pointer to face list */

	      ptr = (long)RDB[cgns1 + IFC_TET_MSH_PTR_FACES];
	      
	      /* Put face */

	      WDB[ptr + 3] = -(double)nfidx;

	      /* Get pointer to side list */

	      ptr = (long)RDB[cgns1 + IFC_TET_MSH_PTR_SIDES];

	      /* Put side */

	      WDB[ptr + 3] = -(double)side;

	    }

	  /* To the face in this cell sharing p0 and p1 with this face */

	  else if (k == 2)
	    {

	      /* Check if the neighbor cell has already been created */
	      /* If it has, it should be in subcells-list */

	      for(l = 0; l < nss; l++)
		{

		  cgns1 = subcells[l];
		  CheckPointer(FUNCTION_NAME, "(ptrss)", DATA_ARRAY, ptr);

		  /* Get pointer to face list */

		  ptr = (long)RDB[cgns1 + IFC_TET_MSH_PTR_FACES];
		  CheckPointer(FUNCTION_NAME, "(ptrfl)", DATA_ARRAY, ptr);

		  /* Get the third face (this is still negative) */

		  idx = -(long)RDB[ptr + 2];

		  /* Get side of other cell */

		  ptr = (long)RDB[cgns1 + IFC_TET_MSH_PTR_SIDES];
		  CheckPointer(FUNCTION_NAME, "(ptrfl)", DATA_ARRAY, ptr);

		  side2 = (long)RDB[ptr + 2];

		  if(idx < 0)
		    Die(FUNCTION_NAME, "wrong sing of face index");

		  /* Check if it's points are 0, 3, 1 reversed i.e. 1 3 0 */

		  surf = (long)RDB[ifc + IFC_PTR_SURF_LIST];

		  /* Get pointer to surface */

		  surf = ListPtr(surf, idx);
		  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
		  /* Get pointer to surface parameters */

		  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];

		  if(((long)RDB[ptr + 0] == p1) &&
		     ((long)RDB[ptr + 1] == p3) &&
		     ((long)RDB[ptr + 2] == p0))
		    break;
		  else if (((long)RDB[ptr + 0] == p0) &&
		     ((long)RDB[ptr + 1] == p3) &&
		     ((long)RDB[ptr + 2] == p1))
		    break;

		}

	      /* Check if found */

	      if (l >= nss)
		{
		  nsurf = (long)RDB[ifc + IFC_PTR_SURF_LIST];

		  /* Get new surface */

		  nsurf = ListPtr(nsurf, *sdone);

		  *sdone = *sdone + 1;

		  /* Get pointer to surface parameters */

		  ptr = (long)RDB[nsurf + SURFACE_PTR_PARAMS];

		  /* Put the points for the parameters */

		  WDB[ptr + 0] = (double)p0;
		  WDB[ptr + 1] = (double)p3;
		  WDB[ptr + 2] = (double)p1;

		  /* This surface's index can be found out from list size */

		  nfidx = *sdone - 1;

		  /* Put face */

		  WDB[nfaces + 2] = -(double)nfidx;

		  /* Put side */

		  WDB[nsides + 2] = (double)side;

		}
	      else
		{
		  /* Surface already created */

		  /* Other cell is now in cgns1 */
		  /* Face is now in surf, face index in idx */

		  /* Put surface to this cells surface list */

		  /* Put face */

		  WDB[nfaces + 2] = -(double)idx;

		  /* Put side */

		  WDB[nsides + 2] = -(double)side2;
		  
		  /* Put cell pointers to owner and neighbor lists */

		  if(-side2 == -1)
		    {
		      WDB[nownrs + idx] = (double)ncgns;
		      WDB[nnbrs  + idx] = (double)cgns1;
		    }
		  else
		    {
		      WDB[nownrs + idx] = (double)cgns1;
		      WDB[nnbrs  + idx] = (double)ncgns;
		    }

		}	  
	    }

	}

    }

  /* Check that new cells are a loop */


  /* Copy created cells to newcells  */
  /* so that they can be accessed in */
  /* DividePolyhedCell */

  for(i=0; i<np;i++)
    {
      newcells[i] = tempcells[i];

      /* Get next cell */

      if(i < np - 1)
	cgns1 = tempcells[i+1];
      else
	cgns1 = tempcells[0];

      /* Get first surface of cell i and third of cell i+1 */

      ptr = (long)RDB[tempcells[i] + IFC_TET_MSH_PTR_FACES];
      ptr1 = (long)RDB[cgns1 + IFC_TET_MSH_PTR_FACES];

      if((long)RDB[ptr + 1] != (long)RDB[ptr1 + 3])
	Die(FUNCTION_NAME, "1 != 3");

      /* Get surface index */

      idx = (long)RDB[ptr + 1];      

      /* Check sides */

      ptr = (long)RDB[tempcells[i] + IFC_TET_MSH_PTR_SIDES];
      ptr1 = (long)RDB[cgns1 + IFC_TET_MSH_PTR_SIDES];

      if((long)RDB[ptr + 1] + (long)RDB[ptr1 + 3] != 0)
	Die(FUNCTION_NAME, "sides dont match");

      side2 = (long)RDB[ptr + 1];

      /* Check neighbour / owner */

      if(side2 == -1)
	{
	  /*  tempcells[i] should be the owner of the face */

	  if((long)RDB[nownrs - idx] != tempcells[i])
	    Die(FUNCTION_NAME, "inside, face %ld owner wrong %ld (s.b. %ld)", idx, (long)RDB[nownrs - idx], tempcells[i]);
	  else if((long)RDB[nnbrs - idx] != cgns1)
	    Die(FUNCTION_NAME, "inside, neighbour wrong");
	}
      else
	{
	  /*  tempcells[i] should be the neighbour of the face */

	  if((long)RDB[nnbrs - idx] != tempcells[i])
	    Die(FUNCTION_NAME, "outside, neighbour wrong");
	  else if((long)RDB[nownrs - idx] != cgns1)
	    Die(FUNCTION_NAME, "outside, owner wrong");

	}


    }
  Mem(MEM_FREE, tempcells);
  Mem(MEM_FREE, nbrcells);

}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
