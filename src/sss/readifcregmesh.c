#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcregmesh.c                               */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2016/02/17 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Reads regular mesh based multi-physics interfaces            */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCRegMesh:"

/*****************************************************************************/

void ReadIFCRegMesh(long loc0, long update)
{
  long loc1, msh, ptr, type, np, n;
  long nx, ny, nz, i, j, k, nr;
  double zmin, zmax, dmax, Tmax, Tmin;
  double d, T, *lims;
  char tmpstr[MAX_STR];
  FILE *fp;

  if(update)
    Warn(FUNCTION_NAME, 
	 "Interface updating for type %ld not yet well tested", 
	 (long)RDB[loc0 + IFC_TYPE]);
  
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

  if (!update)
    WDB[loc0 + IFC_PTR_MAT] = (double)PutText(tmpstr);

  WDB[loc0 + IFC_CALC_OUTPUT] = (double)n;

  /* Read output file name and axial binning for output */

  if (n == YES)
    {
      /* Read file name and binning for other types */

      if (fscanf(fp, "%s %ld %lf %lf %ld", tmpstr, &nz, &zmin, 
		 &zmax, &nr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      if(!update)
	WDB[loc0 + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);
	      
      /* Put number of axial zones */
	      
      if (nz < 1)
	Error(loc0, "Error in number of axial zones");
      else
	WDB[loc0 + IFC_NZ] = (double)nz;
	      
      /* Put axial limits */

      if (zmin < zmax)
	{
	  WDB[loc0 + IFC_ZMIN] = zmin;
	  WDB[loc0 + IFC_ZMAX] = zmax;
	}
      else
	Error(loc0, "Error in axial boundaries");

      /* Put number of radial zones */
	      
      if (nr < 1)
	Error(loc0, "Error in number of radial zones");
      else
	WDB[loc0 + IFC_NR] = (double)nr;
	
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
    
  /* Reset limiting values */

  Tmax = -INFTY;
  Tmin = INFTY;
  dmax = 0.0;
	  
  /* Get mesh type */

  if (fscanf(fp, "%ld", &n) == EOF)
    Die(FUNCTION_NAME, "fscanf error");
	  
  /* Avoid compiler warning */

  lims = NULL;

  /* Check type */

  if (n == MESH_TYPE_CARTESIAN)
    {
      /* Cartesian, allocate memory for limits */

      lims = (double *)Mem(MEM_ALLOC, 6, sizeof(double));

      /* Read data */
	  
      if (fscanf(fp, "%ld %lf %lf %ld %lf %lf %ld %lf %lf",
		 &nx, &lims[0], &lims[1], &ny, &lims[2], &lims[3], 
		 &nz, &lims[4], &lims[5]) == EOF)
	Die(FUNCTION_NAME, "fscanf error");
    }
  else if ((n == MESH_TYPE_HEXX) || (n == MESH_TYPE_HEXY))
    {
      /* Cartesian, allocate memory for limits */

      lims = (double *)Mem(MEM_ALLOC, 6, sizeof(double));

      /* Read data */
	  
      if (fscanf(fp, "%lf %lf %lf %lf %lf %ld %ld %ld",
		 &lims[0], &lims[2], &lims[1], &lims[4], &lims[5],
		 &nx, &ny, &nz) == EOF)
	Die(FUNCTION_NAME, "fscanf error");
    }
  else if (n == MESH_TYPE_ORTHOGONAL)
    {
      /* Orthogonal mesh */

      if (fscanf(fp, "%ld %ld %ld", &nx, &ny, &nz) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      /* Allocate memory */

      lims = (double *)Mem(MEM_ALLOC, nx + ny + nz + 3, sizeof(double));

      /* Read data (NOTE: noiden järjestys pitäis tarkistaa) */
	      
      for (i = 0; i < nx + ny + nz + 3; i++)
	if (fscanf(fp, "%lf", &lims[i]) == EOF)
	  Die(FUNCTION_NAME, "fscanf error");
    }
  else
    Error(loc0, "Invalid mesh type %ld", n);

  /* Calculate number of values */

  np = nx*ny*nz;

  /* Check */

  if (np < 1)
    Error(loc0, "Zero mesh size");
      
  if(!update)
    {
      /* Create mesh structure */
      
      msh = CreateMesh(n, MESH_CONTENT_PTR, -1, nx, ny, nz, lims, -1);
	  
      /* Put pointer */
      
      WDB[loc0 + IFC_PTR_SEARCH_MESH] = (double)msh;
    }
  else
    {
      /* Get pointer to search mesh */

      msh = (long)RDB[loc0 + IFC_PTR_SEARCH_MESH];
    }

  /* Loop over points and read data */
      
  n = 0;

  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++)
	{
	  /* Read values */
		  
	  if (fscanf(fp, "%lf %lf", &d, &T) == EOF)
	    Die(FUNCTION_NAME, "fscanf error %ld %ld %ld", k, j, i);

	  /* Compare to limits */
		  
	  if (fabs(d) > fabs(dmax))
	    dmax = d;
		  
	  if (T > Tmax)
	    Tmax = T;

	  if (T < Tmin)
	    Tmin = T;

	  /* Get mesh pointer */
		      
	  ptr = (long)RDB[msh + MESH_PTR_PTR] + n;
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	
	  if(!update)
	    {
	      /* Allocate memory for data */
	  
	      loc1 = NewItem(ptr, IFC_REG_MSH_LIST_BLOCK_SIZE);
	    }
	  else
	    {
	      /* Get pointer to data */

	      loc1 = (long)RDB[ptr];

	    }	      

	  /* Put data */
		  
	  WDB[loc1 + IFC_REG_MSH_DF] = d;
	  WDB[loc1 + IFC_REG_MSH_TMP] = T;
		  
	  /* Update index */
		  
	  n++;
	}

  /* Put maximum density and temperature                  */
  /* For updates these are checked in processifcregmesh.c */

  WDB[loc0 + IFC_MAX_DENSITY] = dmax;
  WDB[loc0 + IFC_MAX_TEMP] = Tmax;
  WDB[loc0 + IFC_MIN_TEMP] = Tmin;

  /* Set TMS on */

  if (Tmax > 0.0)
    if(!update)
      WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;
	  
  /* Free allocated memory */

  Mem(MEM_FREE, lims);

  /* Close file */

  fclose(fp);

}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
