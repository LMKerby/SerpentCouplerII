/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkpolyhedmesh.c                             */
/*                                                                           */
/* Created:       2015/01/12 (VVa)                                           */
/* Last modified: 2016/02/01 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Does some sanity checks to a divided polyhedral mesh         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckPolyhedMesh:"

/*****************************************************************************/

void CheckPolyhedMesh(long ifc)
{
  long cgns, ownrs, nbrs, idx, loc1, surf, surflist;
  long i,j,k, n;
  long nc, nf, np, ns, fail, failnbr;
  long ptr, pt; 
  double x,y,z, ncalc;

  /* Number of faces per cell and number of points per face is already */
  /* checked in fixpolyhedmesh.c                                       */

  /* Check neighbour relations */

  /* Get owner list */

  ownrs = (long)WDB[ifc + IFC_PTR_OWNR_LIST];
  CheckPointer(FUNCTION_NAME, "(ownrs)", DATA_ARRAY, ownrs);

  /* Get neighbour list */

  nbrs =  (long)WDB[ifc + IFC_PTR_NBR_LIST];
  CheckPointer(FUNCTION_NAME, "(nbrs)", DATA_ARRAY, nbrs);

  /* Get pointer to surface list */

  surf = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

  idx = 0;

  /* Loop over surfaces */
  failnbr = 0;
  ns = 0;

  while (surf > VALID_PTR)
    {    

      /**********************/
      /* Check neighbourity */
      /**********************/

      /* Get pointer to owner cell */

      cgns = (long)RDB[ownrs + idx];
      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      /* Get pointer to face list */
      
      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Get number of faces */

      nf = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Loop over face list to find this surfaces index */

      for (i = 0; i < nf; i++)
	if ((long)RDB[loc1 + i] == idx)
	  break;

      /* Check if this surface was not in owners face list */

      if (i == nf)
	{
#ifdef DEBUG
	  Warn(FUNCTION_NAME, "Could not find face from owner");
#endif

	  /* Count to failed */
	  
	  failnbr++;

	  /* Next surface */

	  surf = NextItem(surf);
	  continue;
	}

      /* Get pointer to side list */
      
      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Check that side for owner is -1 */

      if ((long)RDB[loc1 + i] != -1)
	{
#ifdef DEBUG
	  Warn(FUNCTION_NAME, "Side for owner is not -1");
#endif

	  /* Count to failed */
	  
	  failnbr++;

	  /* Next surface */

	  surf = NextItem(surf);
	  continue;
	}
      /* Get pointer to neighbour cell */

      cgns = (long)RDB[nbrs + idx];

      /* If neighbour cell is -1 there is no neighbour */

      if(cgns == -1)
	{

	  /* Handle next surface */

	  idx++;

	  surf = NextItem(surf);

	  continue;
	}

      CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

      /* Get pointer to face list */
      
      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over face list to find this surfaces index */

      for (i = 0; i < nf; i++)
	if ((long)RDB[loc1 + i] == idx)
	  break;
    
      if (i == nf)
	{
#ifdef DEBUG
	  fprintf(out, "Surface index %ld not owned by owner and neighbour\n", 
		  idx);
	  /* Get pointer to owner cell */

	  cgns = (long)RDB[ownrs + idx];
	  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

	  /* Get pointer to face list */
      
	  loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  
	  fprintf(out, "Surfaces for owner:\n");
	  for (i = 0; i < nf; i++)
	    fprintf(out, "%ld\n", (long)RDB[loc1 + i]);
	  
	  
	  cgns = (long)RDB[nbrs + idx];
	  CheckPointer(FUNCTION_NAME, "(cgns)", DATA_ARRAY, cgns);

	  /* Get pointer to face list */
      
	  loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  
	  fprintf(out, "Surfaces for other:\n");
	  for (i = 0; i < nf; i++)
	    fprintf(out, "%ld\n", (long)RDB[loc1 + i]);
	  
#endif 	  

	  /* Count to failed */
	  
	  failnbr++;

	  /* Cycle loop */

	  surf = NextItem(surf);
	  continue;

	}

      /* Get pointer to side list */
      
      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Check that side for neighbour is 1 */

      if ((long)RDB[loc1 + i] != 1)
	{
#ifdef DEBUG
	  Warn(FUNCTION_NAME, "Side for neighbour is not 1");
#endif

	  /* Count to failed */
	  
	  failnbr++;

	  /* Next surface */

	  surf = NextItem(surf);
	  continue;
	}

      idx++;

      ns++;

      surf = NextItem(surf);

    }

  /* Loop over cells */

  nc = 0;

  fail = 0;

  /* Get pointer to surface list */

  surflist = (long)RDB[ifc + IFC_PTR_SURF_LIST];
  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

  cgns = (long)RDB[ifc + IFC_PTR_TET_MSH];

  while (cgns > VALID_PTR)
    {

      /* Get number of faces */

      nf = (long)RDB[cgns + IFC_TET_MSH_NF];

      /* Get pointer to face list */
      
      loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop over faces to check centerpoints */

      x = 0.0;
      y = 0.0;
      z = 0.0;

      ncalc = 1.0;

      for (j = 0; j < nf; j++)
	{

	  /* Get index of face */

	  n = (long)RDB[loc1 + j];

	  /* Get pointer to surface */

	  surf = ListPtr(surflist, n);
	  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

	  /* Get number of points on the face */

	  np = (long)RDB[surf + SURFACE_N_PARAMS];

	  /* Get pointer to surface parameters */

	  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Calculate face centerpoint */

	  for (k = 0; k < np; k++)
	    {

	      /* Get pointer to beginning of point */

	      pt = (long)RDB[ptr + k];

	      /* Loop over xyz and average for cell centerpoint */

	      x = x*(ncalc - 1)/ncalc + RDB[pt + 0]/ncalc;
	      y = y*(ncalc - 1)/ncalc + RDB[pt + 1]/ncalc;
	      z = z*(ncalc - 1)/ncalc + RDB[pt + 2]/ncalc;

	      ncalc++;
	    }

	}

      if (!InTetCell(ifc, cgns, x, y, z, YES, 0))
	{

#ifdef DEBUG

	  /* Get number of faces */

	  nf = (long)RDB[cgns + IFC_TET_MSH_NF];

	  /* Get pointer to face list */
      
	  loc1 = (long)WDB[cgns + IFC_TET_MSH_PTR_FACES];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  /* Loop over faces and print points */

	  ncalc = 1.0;

	  for (j = 0; j < nf; j++)
	    {

	      /* Get index of face */

	      n = (long)RDB[loc1 + j];

	      fprintf(out, "F%ld = [ ", j+1);

	      /* Get pointer to surface */

	      surf = ListPtr(surflist, n);
	      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

	      /* Get number of points on the face */

	      np = (long)RDB[surf + SURFACE_N_PARAMS];

	      /* Get pointer to surface parameters */

	      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Calculate face centerpoint */

	      for (k = 0; k < np; k++)
		{

		  /* Get pointer to beginning of point */

		  pt = (long)RDB[ptr + k];

		  fprintf(out, "%E %E %E; ", 
			  RDB[pt + 0], 
			  RDB[pt + 1], 
			  RDB[pt + 2]);

		}

	      fprintf(out, "]\n");

	    }

	  fprintf(out, "figure();\n");
	  fprintf(out, "hold on;\n");

	  /* Get pointer to side list */
      
	  ptr = (long)WDB[cgns + IFC_TET_MSH_PTR_SIDES];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  for (j = 0; j < nf; j++)
	    {

	      /* Get index of face */

	      n = (long)RDB[loc1 + j];

	      /* Get pointer to surface */

	      surf = ListPtr(surflist, n);
	      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

	      /* Get number of points on the face */

	      np = (long)RDB[surf + SURFACE_N_PARAMS];

	      fprintf(out,
		      "fill3(F%ld(:,1),F%ld(:,2),F%ld(:,3),ones(%ld,1)*%f);\n",
		      j+1, j+1, j+1, np, (1-j*1.0/(2.0*(double)nf))*RDB[ptr + j]);
	    }

	  fprintf(out, "hold off;\n");
	  fprintf(out, "title('Cell %ld')\n", 
		  (long)RDB[cgns + IFC_TET_MSH_IDX]);

	  fprintf(out, "Barycenter (%E %E %E)\n", x, y, z);
	  Warn(FUNCTION_NAME, "Barycenter not inside cell %ld", (long)RDB[cgns + IFC_TET_MSH_IDX]);
#endif
	  fail++;
	}

      /* Increment number of cells tested */

      nc++;

      /* Next cell */

      cgns = NextItem(cgns);
    }

  fprintf(out, " - Checked neighbority for %ld surfaces, %ld failed\n", 
	  ns, failnbr);

  fprintf(out, " - Checked that barycenter is inside tet for %ld cells, %ld failed\n", nc, fail);

  /* Die if there were errors */

  if (failnbr + fail > 0)
    Die(FUNCTION_NAME, "There were errors in the mesh");

  fprintf(out, "\n");
}

/*****************************************************************************/
