/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifctetmesh.c                               */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2014/10/06 (VVa)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Reads tet-mesh multi-physics interfaces                      */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCTetMesh:"

/*****************************************************************************/

void ReadIFCTetMesh(long loc0, long update)
{
  long loc1, loc2, ptr, type, np, nc, nf, n;
  long nd, i, j, k, cell, surf;
  double xmin, xmax, ymin, ymax, zmin, zmax, dmax, Tmax, Tmin;
  double **dat, x, y, z, d, T;
  char tmpstr[MAX_STR];
  FILE *fp;

  if(update == 1)
    {
      Warn(FUNCTION_NAME, "Interface updating for type %ld not yet implemented", RDB[loc0 + IFC_TYPE]);
      return;
    }

  /* Open file for reading */
      
  if ((fp = fopen(GetText(loc0 + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(loc0, "Multi-physics interface file \"%s\" does not exist",
	  GetText(loc0 + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Read material and output flag */

  if (fscanf(fp, "%s %ld", tmpstr, &n) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc0 + IFC_PTR_MAT] = (double)PutText(tmpstr);
  WDB[loc0 + IFC_CALC_OUTPUT] = (double)n;

  /* Read output file name and axial binning for output */

  if (n == YES)
    {
	  
      /* Read only file name */

      if (fscanf(fp, "%s", tmpstr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[loc0 + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

      /* Loop over previous and check duplicate file names */

      loc1 = PrevItem(loc0);
      while (loc1 > VALID_PTR)
	{
	  /* Compare file names */

	  if ((long)RDB[loc1 + IFC_PTR_OUTPUT_FNAME] > VALID_PTR)
	    if (CompareStr(loc0 + IFC_PTR_OUTPUT_FNAME, 
			   loc1 + IFC_PTR_OUTPUT_FNAME))
	      Error(loc0, 
		    "Duplicate output file name with distribution \"%s\"",
		    GetText(loc1 + IFC_PTR_INPUT_FNAME));
	      
	  /* Pointer to previous */

	  loc1 = PrevItem(loc1);
	}	  
    }

    

  /***********************************************************************/

  /***** Unstructured tetrahedral mesh *******************************/

  /* NOTE: this is more of a prototype of the tet mesh structure */

  /* Reset mesh boundaries */
	  
  xmin =  INFTY;
  xmax = -INFTY;
  ymin =  INFTY;
  ymax = -INFTY;
  zmin =  INFTY;
  zmax = -INFTY;

  /* Read number of points */

  if (fscanf(fp, "%ld", &nd) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Allocate memory for data */
	  
  dat = (double **)Mem(MEM_ALLOC, 3, sizeof(double *));
      
  for(n = 0; n < 3; n++)
    dat[n] = (double *)Mem(MEM_ALLOC, nd, sizeof(double));

  /* Read points */

  for (n = 0; n < nd; n++)
    {
      /* Read coordinates */

      if (fscanf(fp, "%lf %lf %lf", &x, &y, &z) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      /* Put data */

      dat[0][n] = x;
      dat[1][n] = y;
      dat[2][n] = z;

      /* Compare to limits */

      if (x < xmin)
	xmin = x;
      if (x> xmax)
	xmax = x;
	      
      if (y < ymin)
	ymin = y;
      if (y > ymax)
	ymax = y;
	      
      if (z < zmin)
	zmin = z;
      if (z > zmax)
	zmax = z;
    }	      

  /* Put boundaries */

  xmin = -185;
  xmax =  185;
  ymin = -185;
  ymax =  185;
  zmin = 0.0;
  zmax = 370;

  WDB[loc0 + IFC_MESH_XMIN] = xmin;
  WDB[loc0 + IFC_MESH_XMAX] = xmax;
  WDB[loc0 + IFC_MESH_YMIN] = ymin;
  WDB[loc0 + IFC_MESH_YMAX] = ymax;
  WDB[loc0 + IFC_MESH_ZMIN] = zmin;
  WDB[loc0 + IFC_MESH_ZMAX] = zmax;

  /* Read number of cells */

  if (fscanf(fp, "%ld", &nc) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Reset limiting values */

  Tmax = -INFTY;
  Tmin = INFTY;
  dmax = 0.0; 

  /* Loop over cells */

  for (n = 0; n < nc; n++)
    {
      /* Allocate memory for tet cell */
		  
      loc1 = NewItem(loc0 + IFC_PTR_TET_MSH, 
		     IFC_TET_MSH_LIST_BLOCK_SIZE);

      /* Put index */

      WDB[loc1 + IFC_TET_MSH_IDX] = (double)(n + 1);

      /* Read density, temperature, number of surfaces and output */
      /* flag. */

      if (fscanf(fp, "%lf %lf %ld %ld", &d, &T, &nf, &i) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      /* Put values */

      WDB[loc1 + IFC_TET_MSH_DF] = d;
      WDB[loc1 + IFC_TET_MSH_TMP] = T;
      WDB[loc1 + IFC_TET_MSH_STAT_IDX] = (double)i;

      /* Compare to maxima */

      if (fabs(d) > fabs(dmax))
	dmax = d;
	      
      if (T > Tmax)
	Tmax = T;

      if (T < Tmin)
	Tmin = T;

      /* Create (geometry) cell */

      cell = NewItem(loc1 + IFC_TET_MSH_PTR_CELL, CELL_BLOCK_SIZE);

      /* Reset surface list pointer */
	  
      loc2 = -1;

      /* Reset cell boundaries */

      xmin =  INFTY;
      xmax = -INFTY;
      ymin =  INFTY;
      ymax = -INFTY;
      zmin =  INFTY;
      zmax = -INFTY;

      /* Loop over facets */
	      
      for (i = 0; i < nf; i++)
	{
	  /* Create new item in cell intersection list */
	  
	  loc2 = NewItem(cell + CELL_PTR_SURF_INSC, 
			 CELL_INSC_BLOCK_SIZE);
		  		  
	  /* Put side */
	  
	  WDB[loc2 + CELL_INSC_SIDE] = -1.0;
	  
	  /* Allocate memory for counter from private array */
	  
	  WDB[loc2 + CELL_INSC_PTR_OUT_COUNT] =
	    (double)AllocPrivateData(1, PRIVA_ARRAY);

	  /* Create surface */

	  surf = NewItem(DATA_PTR_S0, SURFACE_BLOCK_SIZE);
	  WDB[loc2 + CELL_INSC_PTR_SURF] = (double)surf;

	  /* put surface type and number of parameters */

	  WDB[surf + SURFACE_TYPE] = (double)SURF_PLANE;
	  WDB[surf + SURFACE_N_PARAMS] = 9.0;
		  
	  /* Set used-flag */

	  SetOption(surf + SURFACE_OPTIONS, OPT_USED);

	  /* Allocate memory for parameters */

	  ptr = ReallocMem(DATA_ARRAY, 9);
	  WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

	  /* Read number of points */

	  if (fscanf(fp, "%ld", &np) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  /* Check type */

	  if (np < 3)
	    Die(FUNCTION_NAME, "Not enough points");

	  /* Read points */
	
	  for (j = 0; j < np; j++)
	    {
	      /* Read point index */

	      if (fscanf(fp, "%ld", &k) == EOF)
		Die(FUNCTION_NAME, "fscanf error");

	      /* Check */

	      if ((k < 1) || (k > nd))
		Die(FUNCTION_NAME, "Invalid point index %ld", k);
		      
	      /* Get point */

	      x = dat[0][k - 1];
	      y = dat[1][k - 1];
	      z = dat[2][k - 1];

	      /* Put surface parameters */

	      if (j < 3)
		{
		  WDB[ptr++] = x;
		  WDB[ptr++] = y;
		  WDB[ptr++] = z;
		}

	      /* Compare to limits */

	      if (x < xmin)
		xmin = x;
	      if (x> xmax)
		xmax = x;
		      
	      if (y < ymin)
		ymin = y;
	      if (y > ymax)
		ymax = y;
		      
	      if (z < zmin)
		zmin = z;
	      if (z > zmax)
		zmax = z;
	    }
	}
	     
      /* Check pointer and close list */
	      
      if (loc2 < VALID_PTR)
	Die(FUNCTION_NAME, "No surfaces");
      else
	CloseList(loc2);

      /* Put cell boundaries */
	      
      WDB[loc1 + IFC_TET_MSH_XMIN] = xmin;
      WDB[loc1 + IFC_TET_MSH_XMAX] = xmax;
      WDB[loc1 + IFC_TET_MSH_YMIN] = ymin;
      WDB[loc1 + IFC_TET_MSH_YMAX] = ymax;
      WDB[loc1 + IFC_TET_MSH_ZMIN] = zmin;
      WDB[loc1 + IFC_TET_MSH_ZMAX] = zmax;
    }
	  
  /* Put maximum temperature and density */
	  
  WDB[loc0 + IFC_MAX_DENSITY] = dmax;
  WDB[loc0 + IFC_MAX_TEMP] = Tmax;
  WDB[loc0 + IFC_MIN_TEMP] = Tmin;

  /* Put search mesh size (tota formaattia pitää muuttaa että tän */
  /* saa luettua) */

  WDB[loc0 + IFC_SEARCH_MESH_NX] = 50.0;
  WDB[loc0 + IFC_SEARCH_MESH_NY] = 50.0;
  WDB[loc0 + IFC_SEARCH_MESH_NZ] = 50.0;

  /* Free temporary data */
      
  for(n = 0; n < 3; n++)
    Mem(MEM_FREE, dat[n]);
	  
  Mem(MEM_FREE, dat);
  dat = NULL;

  /* Set TMS on */

  if (Tmax > 0.0)
    WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;

  /*******************************************************************/

  /* Close file */

  fclose(fp);

}

/*****************************************************************************/
