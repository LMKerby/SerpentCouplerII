/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcofmesh.c                                */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2016/01/31 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Reads OpenFOAM multi-physics interfaces                      */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*             - Tiheysjakauman yksiköt luetaan nyt dimensions -vektorista   */
/*               (2.12.2014 / 2.1.23 / JLe)                                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCOFMesh:"

/*****************************************************************************/

void ReadIFCOFMesh(long loc0, long update)
{
  long loc1, ptr, ptr1, type, np, nc, nf, n, p;
  long nd, i, j, k, cell, surf, fi, nmax;
  long nbhr, divide, dtype, ttype, dim[7], *matarr, nmat;
  long pointlist, solist, snlist, sslist, **surfaces, **sides, *nsurf;
  long cellpts, surflist;
  long nhexc, ntetc, npyramc, nprismc, npolyc;
  long nfacepts[3];
  double xmin, xmax, ymin, ymax, zmin, zmax, dmax, dmin, Tmax, Tmin;
  double x, y, z, d, T, d0, T0, mul;
  char tmpstr[MAX_STR], pfile[MAX_STR], ffile[MAX_STR], ofile[MAX_STR];
  char nfile[MAX_STR], rfile[MAX_STR], tfile[MAX_STR], mapfile[MAX_STR], *line;
  char matfile[MAX_STR], mname[MAX_STR], name[MAX_STR];
  FILE *fp, *fp1;

  /* Open file for reading */
      
  if ((fp = fopen(GetText(loc0 + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(loc0, "Multi-physics interface file \"%s\" does not exist",
	  GetText(loc0 + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Read material name or uni & BG uni */
  if(type != IFC_TYPE_OF_SOLID)
    {
      if (fscanf(fp, "%s", tmpstr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      if(!update)
	WDB[loc0 + IFC_PTR_MAT] = (double)PutText(tmpstr);
    }
  else
    {
      /* Read universe name (stored in readinput) */
      
      if (fscanf(fp, "%s", tmpstr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      /* Read BG universe name (stored in readinput) */

      if (fscanf(fp, "%s", tmpstr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

    }

  /* Read output flag */

  if (fscanf(fp, "%ld", &n) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Read output file name for output */

  if (n == YES)
    {
      /* Set flag */

      WDB[loc0 + IFC_CALC_OUTPUT] = (double)YES;
 	  
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
  else
    {
      /* Reset flag */

      WDB[loc0 + IFC_CALC_OUTPUT] = (double)NO;
    }
    
  /*******************************************************************/

  /***** OpenFoam mesh ***********************************************/

  /* Avoid compiler warning */

  dtype = -1;
  ttype = -1;

  /* Reset divisor flag */

  divide = NO;
  
  /* Reset number of different cell types */

  nhexc = 0;
  ntetc = 0;
  nprismc = 0;
  npolyc = 0;      
  npyramc = 0;

  /* Read nominal temperature and density */

  if (fscanf(fp, "%lf %lf", &d0, &T0) == EOF)
    Error(loc0, "Missing nominal temperature or density");

  /* Read search mesh size */

  if (fscanf(fp, "%ld", &i) == EOF)
    Die(FUNCTION_NAME, "fscanf error (8)");
  else if (i < 2)
    Error(loc0, "Error in search mesh split");
  else
    {
      if(!update)
	WDB[loc0 + IFC_SEARCH_MESH_ADA_SPLIT] = (double)i;
    }

  if (fscanf(fp, "%ld", &i) == EOF)
    Die(FUNCTION_NAME, "fscanf error (9)");
  else if ((i < 1) || (i > 10))
    Error(loc0, "Error in search mesh depth");

  /* Allocate memory for search mesh size */

  if (!update)
    {
      ptr = ReallocMem(DATA_ARRAY, i + 1);
      WDB[loc0 + IFC_SEARCH_MESH_ADA_PTR_SZ] = (double)ptr;
    }
  else
    ptr = -1;

  /* Read sizes */

  for (n = 0; n < i; n++)
    {
      if (fscanf(fp, "%ld", &j) == EOF)
	Die(FUNCTION_NAME, "fscanf error (10)");
      else if (j < 1)
	Error(loc0, "Error in search mesh size");
      else
	{
	  if(!update)
	    WDB[ptr++] = (double)j;
	}
    }
	  
  /* Put terminator */

  if(!update)
    WDB[ptr] = -1.0;

  /* Read file names (points, faces, owner and neighbour) */

  if (fscanf(fp, "%s", pfile) == EOF)
    Error(loc0, "Missing path to points file");
  else
    TestDOSFile(pfile);

  WDB[loc0 + IFC_PTR_OF_PFILE] = PutText(pfile);

  if (fscanf(fp, "%s", ffile) == EOF)
    Error(loc0, "Missing path to faces file");
  else
    TestDOSFile(ffile);

  WDB[loc0 + IFC_PTR_OF_FFILE] = PutText(ffile);

  if (fscanf(fp, "%s", ofile) == EOF)
    Error(loc0, "Missing path to owner file");
  else
    TestDOSFile(ofile);

  WDB[loc0 + IFC_PTR_OF_OFILE] = PutText(ofile);

  if (fscanf(fp, "%s", nfile) == EOF)
    Error(loc0, "Missing path to neighbour file");
  else
    TestDOSFile(nfile);

  WDB[loc0 + IFC_PTR_OF_NFILE] = PutText(nfile);

  if ((type == IFC_TYPE_OF_MAT) || (type == IFC_TYPE_OF_SOLID))
    {

      /* Read material file */

      if (fscanf(fp, "%s", matfile) == EOF)
	Error(loc0, "Missing path to material file");
      else
	TestDOSFile(matfile);

      WDB[loc0 + IFC_PTR_OF_MFILE] = PutText(matfile);

    }

  /* Read density file */

  if (fscanf(fp, "%s", rfile) == EOF)
    Error(loc0, "Missing path to density file");
  else
    TestDOSFile(rfile);

  WDB[loc0 + IFC_PTR_OF_RFILE] = PutText(rfile);

  /* Read type flag */

  if (strcmp(rfile, "-1"))
    {
      if (fscanf(fp, "%ld", &dtype) == EOF)
	Error(loc0, "Missing density distribution type flag");
      else if (dtype != 1)
	Error(loc0, "Invalid density distribution type flag");
    }

  /* Read temperature file */

  if (fscanf(fp, "%s", tfile) == EOF)
    Error(loc0, "Missing path to temperature file");
  else
    TestDOSFile(tfile);

  WDB[loc0 + IFC_PTR_OF_TFILE] = PutText(tfile);

  /* Read type flag */

  if (strcmp(tfile, "-1"))
    {
      if (fscanf(fp, "%ld", &ttype) == EOF)
	Error(loc0, "Missing temperature distribution type flag");
      else if ((ttype != 1) && (ttype != 2))
	Error(loc0, "Invalid temperature distribution type flag");
    }

  /* Read output mapping file */

  if ((long)RDB[loc0 + IFC_CALC_OUTPUT] == YES)
    {
      if (fscanf(fp, "%s", mapfile) == EOF)
	Error(loc0, "Missing path to output map file");
      else
	TestDOSFile(mapfile);
    }

  /* Close file */

  fclose(fp);

  /*******************************************************************/

  /***** Read points *************************************************/

  /* Reset mesh boundaries */
	  
  xmin =  INFTY;
  xmax = -INFTY;
  ymin =  INFTY;
  ymax = -INFTY;
  zmin =  INFTY;
  zmax = -INFTY;

  /* Open points file for reading */
  if(!update)
    {


      /* Check file format */

      TestDOSFile(pfile);

      /* Open points file for reading */
      
      if ((fp = fopen(pfile, "r")) == NULL)
	Error(loc0, "Points file \"%s\" does not exist", 
	      pfile);

      /* Read header data */

      ReadOFHeader(fp, &n, &nd, (long *)dim);

      /* Check number of points */

      CheckValue(FUNCTION_NAME, "nd", "", nd, 4, 100000000000);

      /* Store number of points */

      WDB[loc0 + IFC_NP] = (double)nd;

      /* Allocate memory for points */
      
      pointlist = ReallocMem(DATA_ARRAY, nd*3);

      /* Store pointer to first point */
      
      WDB[loc0 + IFC_PTR_POINT_LIST] = (double)pointlist;

      /* Read points */
      
      for (n = 0; n < nd; n++)
	{
	  /* Read coordinates */

	  line = ReadOFData(fp, OF_FILE_POINTS);
	      
	  if (sscanf(line, "%lf %lf %lf", &x, &y, &z) == EOF)
	    Error(loc0, "Not enough entries in points file");
	  
	  /* Convert to cm */
	  
	  x = x*100.0;
	  y = y*100.0;
	  z = z*100.0;
	  
	  /* Put data */
	  
	  WDB[pointlist++] = x;
	  WDB[pointlist++] = y;
	  WDB[pointlist++] = z;
	  
	  /* Compare to limits */
	  
	  if (x < xmin)
	    xmin = x;
	  if (x > xmax)
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
      
      /* Close file */

      fclose(fp);

      /* Print out points */

      /* Get number of points */

      nd = (long)RDB[loc0 + IFC_NP];

      /* Get pointer to first point */

      pointlist = (long)RDB[loc0 + IFC_PTR_POINT_LIST];

      /* Loop and print points*/
      /*
      for (i = 0; i < nd; i++)
	printf("Point (%f %f %f)\n", RDB[pointlist + i*3], RDB[pointlist + i*3 + 1], RDB[pointlist + i*3 + 2]);
      */
      /* Store boundaries */
 
      WDB[loc0 + IFC_MESH_XMIN] = xmin;
      WDB[loc0 + IFC_MESH_XMAX] = xmax;
      WDB[loc0 + IFC_MESH_YMIN] = ymin;
      WDB[loc0 + IFC_MESH_YMAX] = ymax;
      WDB[loc0 + IFC_MESH_ZMIN] = zmin;
      WDB[loc0 + IFC_MESH_ZMAX] = zmax;

    }

  /*******************************************************************/

  /***** Read faces **************************************************/

  /* Open faces file for reading */
  if(!update)
    {

      /* Check file format */

      TestDOSFile(ffile);

      /* Open faces file for reading */
      
      if ((fp = fopen(ffile, "r")) == NULL)
	Error(loc0, "Faces file \"%s\" does not exist", 
	      ffile);

      /* Read header data */

      ReadOFHeader(fp, &n, &nf, (long *)dim);

      /* Check number of faces */

      CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 100000000000);
   
      /* Preallocate memory for data */

      PreallocMem(nf*(SURFACE_BLOCK_SIZE + 9) + 2*nf, DATA_ARRAY);
   
      /* Allocate memory for pointers */
      
      solist = ReallocMem(DATA_ARRAY, nf);

      WDB[loc0 + IFC_PTR_OWNR_LIST] = (double)solist;

      snlist = ReallocMem(DATA_ARRAY, nf);

      WDB[loc0 + IFC_PTR_NBR_LIST] = (double)snlist;

      /* Reset surface pointer */

      surf = -1;

      /* Get pointer to first point */

      pointlist = (long)RDB[loc0 + IFC_PTR_POINT_LIST];

      /* Loop over faces */
	     
      for (n = 0; n < nf; n++)
	{

	  /* Create surface and put pointer to temporary array */
	  
	  surf = NewItem(loc0 + IFC_PTR_SURF_LIST, SURFACE_BLOCK_SIZE);

	  /* Reset pointers */

	  WDB[solist + n] = -1;
	  WDB[snlist + n] = -1;

	  /* Read entry */

	  line = ReadOFData(fp, OF_FILE_FACES);

	  /* Read number of points */
	  
	  p = NextWord(line, tmpstr);
	  line = &line[p];
	  np = (long)atoi(tmpstr);
	
	  /* Check type */
	      
	  if (np < 3)
	    Error(loc0, "Not enough points");
	  else if (np > 3)
	    divide = YES;

	  /* Allocate memory for parameters */
		  
	  ptr = ReallocMem(DATA_ARRAY, np);
	  WDB[surf + SURFACE_PTR_PARAMS] = (double)ptr;

	  /* put surface type and number of parameters */
		  
	  WDB[surf + SURFACE_TYPE] = (double)SURF_PLANE;
	  WDB[surf + SURFACE_N_PARAMS] = (double)(np);

	  /* Read points */
	  
	  for (j = 0; j < np; j++)
	    {
	      /* Read point index */

	      p = NextWord(line, tmpstr);
	      line = &line[p];
	      k = (long)atoi(tmpstr);
		  
	      /* Check */
		  
	      if ((k < 0) || (k > nd - 1))
		Die(FUNCTION_NAME, "Invalid point index %ld", k);

	      /* Store pointer to beginning of point in point list */

	      WDB[ptr++] = (double)(pointlist + k*3);

	    }
	}

      /* Close surface list */

      CloseList(surf);
	  
      /* Close file */
      
      fclose(fp);
#ifdef mmmaaa
      /* Loop over faces and print out */

      surf = (long)RDB[loc0 + IFC_PTR_SURF_LIST];
      pointlist = (long)RDB[loc0 + IFC_PTR_POINT_LIST];

      while (surf > VALID_PTR)
	{
	  printf("\nFace:\n");
	  /* Get number of points */

	  np = (long)RDB[surf + SURFACE_N_PARAMS];

	  /* Get pointer to surface parameters (point list) */
	  
	  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];

	  /* Loop over points */

	  for (i = 0; i < np; i++)
	    {
	      /* Get point */
	      loc1 = (long)RDB[ptr + i];

	      /* Print point */ 
	      printf("Point: %f %f %f\n", RDB[loc1], RDB[loc1 + 1], RDB[loc1 + 2]);
	    }

	  surf = NextItem(surf);
	}

#endif
    }

  /*******************************************************************/

  /***** Read owner and neighbour files ******************************/
  if(!update)
    {

      /* Check file formats */

      TestDOSFile(ofile);
      TestDOSFile(nfile);

      /* Reset number of cells */

      nc = -1;

      /* Repeat for owner and neighbour files */
	  
      for (fi = 1; fi < 3; fi++)
	{
	  /* Check mode */
	  
	  if (fi == 1)
	    {
	      /* Open file */

	      if ((fp = fopen(ofile, "r")) 
		  == NULL)
		Error(loc0, "Owner file \"%s\" does not exist", 
		      ofile);
	      
	      /* Read header data */
	      
	      ReadOFHeader(fp, &n, &nmax, (long *)dim);

	    }
	  else
	    {
	      /* Open file */

	      if ((fp = fopen(nfile, "r")) 
		  == NULL)
		Error(loc0, "Neighbour file \"%s\" does not exist", 
		      nfile);
	      
	      /* Read header data */
	      
	      ReadOFHeader(fp, &n, &nmax, (long *)dim);

	    }
	  
	  /* Loop over faces */
	  
	  for (n = 0; n < nmax; n++)
	    {
	      /* Read cell index from owner/neighbour file */
	      
	      if (fi == 1)
		line = ReadOFData(fp, OF_FILE_OWNER);
	      else
		line = ReadOFData(fp, OF_FILE_NEIGHBOUR);
	      
	      if (sscanf(line, "%ld", &i) == EOF)
		{
		  if (fi == 1)
		    Error(loc0, "Not enough entries in owner file");
		  else
		    Error(loc0, "Not enough entries in neighbour file");
		}
	      
	      /* Update number of cells */
	      
	      if (i + 1 > nc)
		nc = i + 1;
	      
	      /* Put pointer in array */

	      if(fi == 1)
		WDB[solist + n] = (double)i;
	      else
		WDB[snlist + n] = (double)i;
	      
	    }
	  
	  /* Close file */
	  
	  fclose(fp);
	}

      /***********************************************************************/

      /***** Create cell structures ******************************************/
      
      /* Check number of cells */

      CheckValue(FUNCTION_NAME, "nc", "", nc, 1, 100000000000);

      /* Preallocate memory for data */

      PreallocMem(nc*(IFC_TET_MSH_LIST_BLOCK_SIZE), DATA_ARRAY);

      /* Avoid compiler warning */

      loc1 = -1;
	  
      /* Create cells */
      
      for (n = 0; n < nc; n++)
	{
	  /* Allocate memory for tet cell */

	  loc1 = NewItem(loc0 + IFC_PTR_TET_MSH, IFC_TET_MSH_LIST_BLOCK_SIZE);
	  
	  /* Put index */
	  
	  WDB[loc1 + IFC_TET_MSH_IDX] = (double)n;
	  
	  /* Reset cell boundaries */
	  
	  WDB[loc1 + IFC_TET_MSH_XMIN] = INFTY;
	  WDB[loc1 + IFC_TET_MSH_XMAX] = -INFTY;
	  WDB[loc1 + IFC_TET_MSH_YMIN] = INFTY;
	  WDB[loc1 + IFC_TET_MSH_YMAX] = -INFTY;
	  WDB[loc1 + IFC_TET_MSH_ZMIN] = INFTY;
	  WDB[loc1 + IFC_TET_MSH_ZMAX] = -INFTY;
	  
	}
      
      /* Close list */
      
      CloseList(loc1);

      /* Get pointer to surface list */

      sslist = (long)RDB[loc0 + IFC_PTR_SURF_LIST];

      /* Check number of faces */

      CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 100000000000);

      /* Put number of faces */

      WDB[loc0 + IFC_NF] = (double)nf;

      /* Put number of cells */
      /* NOTE: toisen noista voinee poistaa */

      WDB[loc0 + IFC_NC] = (double)nc;

      WDB[loc0 + IFC_MSH_N_CELLS] = (double)nc;

      /* Allocate memory for face lists for cells */

      /* List of cells */

      surfaces = (long **)Mem(MEM_ALLOC, nc, sizeof(long*));
      sides = (long **)Mem(MEM_ALLOC, nc, sizeof(long*));

      /* Allocate memory for four faces initially */

      for (n = 0; n < nc; n++)
	{
	  surfaces[n] = (long *)Mem(MEM_ALLOC, 4, sizeof(long));
	  sides[n] = (long *)Mem(MEM_ALLOC, 4, sizeof(long));
	  /* Reset faces */
	  
	  for (i = 0; i < 4; i++)
	    {
	      surfaces[n][i] = -1;
	      sides[n][i] = -1;
	    }
	}

      /* Number of surfaces per cell */

      nsurf = (long *)Mem(MEM_ALLOC, nc, sizeof(long*));

      /* Loop over faces */

      for (n = 0; n < nf; n++)
	{

	  /* Repeat for owner and neighbour files */
	  
	  for (fi = 1; fi < 3; fi++)
	    {
	      /* Get cell index */
	      if(fi == 1)
		i = (long)RDB[solist + n];
	      else
		i = (long)RDB[snlist + n];
	      
	      /* Check */
	      
	      if (i < 0)
		continue;
	      else if (i > nc - 1)
		Error(loc0, "Cell index exceeds maximum");
	      
	      /* Get pointer to tet cell */
	      
	      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	      
	      loc1 = ListPtr(loc1, i);
	      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);	      

	      /* Store straight pointer of cell to owner/neighbor list */
	      
	      if(fi == 1)
		WDB[solist + n] = loc1;
	      else
		WDB[snlist + n] = loc1;

	      /* Check index */
	      
	      if ((long)RDB[loc1 + IFC_TET_MSH_IDX] != i)
		Die(FUNCTION_NAME, "Indexing error");

	      /* Get neighbour index */

	      if(fi == 1)
		{
		  /* Neighbour list contains index (not pointer) */

		  i = (long)RDB[snlist + n];

		  if (i > -1)
		    {	      
		      if (i > nc - 1)
			Error(loc0, "Cell index exceeds maximum");
		  
		      nbhr = ListPtr(loc1, i);
		      CheckPointer(FUNCTION_NAME, "(nbhr)", DATA_ARRAY, nbhr);

		      /* Check index */
		  
		      if ((long)RDB[nbhr + IFC_TET_MSH_IDX] != i)
			Die(FUNCTION_NAME, "Indexing error (%ld %ld)",
			    (long)RDB[nbhr + IFC_TET_MSH_IDX], i);

		    }
		  else
		    nbhr = -1;
		}
	      else
		{
		  /* Owner list already contains direct pointer */

		  nbhr = (long)RDB[solist + n];
		}	      

	      /* Get pointer to surface */

	      surf = (long)ListPtr(sslist, n);
	      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);

	      /* Put surface pointer to tet cell */

	      /* Get owner / neighbour */

	      if(fi == 1)
		loc1 = (long)RDB[solist + n];
	      else
		loc1 = (long)RDB[snlist + n];

	      /* Get index of owner */

	      i = (long)RDB[loc1 + IFC_TET_MSH_IDX];

	      /* Store index of face to cell's face list */

	      /* Check that cells surface list is not full */


	      if (nsurf[i] >= 4)
		{

		  /* Allocate memory for one more surface */

		  surfaces[i] = (long*)Mem(MEM_REALLOC, surfaces[i], (nsurf[i]+1)*sizeof(long));
		  sides[i] = (long*)Mem(MEM_REALLOC, sides[i], (nsurf[i]+1)*sizeof(long));

		  /* Store face index */

		  surfaces[i][nsurf[i]] = n;

		  if(fi == 1)
		    sides[i][nsurf[i]] = -1;
		  else
		    sides[i][nsurf[i]] = 1;

		  /* Increase number of stored surfaces */

		  nsurf[i]++;

		}
	      else
		{

		  /* Store face pointer */

		  surfaces[i][nsurf[i]] = n;

		  if(fi == 1)
		    sides[i][nsurf[i]] = -1;
		  else
		    sides[i][nsurf[i]] = 1;

		  /* Increase number of stored surfaces */

		  nsurf[i]++;

		}

	      /* Get cell bounding box */
	      
	      xmin = RDB[loc1 + IFC_TET_MSH_XMIN];
	      xmax = RDB[loc1 + IFC_TET_MSH_XMAX];
	      ymin = RDB[loc1 + IFC_TET_MSH_YMIN];
	      ymax = RDB[loc1 + IFC_TET_MSH_YMAX];
	      zmin = RDB[loc1 + IFC_TET_MSH_ZMIN];
	      zmax = RDB[loc1 + IFC_TET_MSH_ZMAX];
	      
	      /* Get pointer to surface parameters */
	      
	      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	      /* Get number of parameters */
	      
	      np = (long)RDB[surf + SURFACE_N_PARAMS];
	      CheckValue(FUNCTION_NAME, "np", "", np, 3, 900);
	      
	      /* Loop over values */
	      
	      for (j = 0; j < np; j++)
		{

		  /* Get point */
		  ptr1 = (long)RDB[ptr + j];

		  /* Get coordinates */
		  
		  x = RDB[ptr1++];
		  y = RDB[ptr1++];
		  z = RDB[ptr1++];
		  
		  /* Compare to bounding box */
		  
		  if (x < xmin)
		    xmin = x;
		  if (x > xmax)
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

	      /* Put cell boundaries */
	      
	      WDB[loc1 + IFC_TET_MSH_XMIN] = xmin;
	      WDB[loc1 + IFC_TET_MSH_XMAX] = xmax;
	      WDB[loc1 + IFC_TET_MSH_YMIN] = ymin;
	      WDB[loc1 + IFC_TET_MSH_YMAX] = ymax;
	      WDB[loc1 + IFC_TET_MSH_ZMIN] = zmin;
	      WDB[loc1 + IFC_TET_MSH_ZMAX] = zmax;
	    }
	}

      /* Put faces to cells */

      pointlist = (long)RDB[loc0 + IFC_PTR_POINT_LIST];

      /* Loop over cell list */

      for (n = 0; n < nc; n++)
	{
	  /* Get pointer to tet cell */

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	      
	  loc1 = ListPtr(loc1, n);
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);	      

	  /* Get number of faces */

	  nd = nsurf[n];

	  /* Put number of faces */

	  WDB[loc1 + IFC_TET_MSH_NF] = (double)nd;

	  /* Allocate memory for face list */

	  ptr = ReallocMem(DATA_ARRAY, nd);

	  /* Put pointer to face list */

	  WDB[loc1 + IFC_TET_MSH_PTR_FACES] = (double)ptr;

	  /* Allocate memory for side list */

	  ptr1 = ReallocMem(DATA_ARRAY, nd);

	  /* Put pointer to side list */

	  WDB[loc1 + IFC_TET_MSH_PTR_SIDES] = (double)ptr1;

	  /* Loop over faces to store them */

	  for (i = 0; i < nd ; i++)
	    {

	      WDB[ptr + i] = surfaces[n][i];
	      WDB[ptr1 + i] = sides[n][i];

	    }

	}

      /* Free temporary lists */

      for (n = 0; n < nc; n++)
	Mem(MEM_FREE, surfaces[n]);

      Mem(MEM_FREE, surfaces);	  

      Mem(MEM_FREE, nsurf);	  

      /***** Read materials **************************************************/
      if ((type == IFC_TYPE_OF_MAT) || (type == IFC_TYPE_OF_SOLID))
	{

	  /* Check file format */

	  TestDOSFile(matfile);

	  /* Open materials file */
      
	  if ((fp = fopen(matfile, "r")) == NULL)
	    Error(loc0, "Materials file \"%s\" does not exist", 
		  matfile);      

	  /* Read header data */

	  ReadOFHeader(fp, &n, &i, (long *)dim);

	  /* Check size */

	  if (i != nc)
	    Error(loc0, "Invalid number of entries in material file");

	  /* Create an array for unique material names */

	  matarr = (long *)Mem(MEM_ALLOC, 1, sizeof(long)); 

	  nmat = 0;

	  /* Loop over cells */

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	  while (loc1 > VALID_PTR)
	    {
	      /* Read material name */
	  
	      line = ReadOFData(fp, OF_FILE_MATERIAL);
	  
	      if (sscanf(line, "%s", mname) == EOF)
		Error(loc0, "Not enough entries in materials file");

	      /* Check for filled cell */

	      if (!strcasecmp(mname, "fill"))
		Error(loc0, "Filled cells not allowed with unstructured meshes");

	      /* Check if the material name has been stored */

	      for(i=0;i<nmat;i++)
		{
		  if(strcmp(&ASCII[matarr[i]],mname) == 0)
		    break;
		}

	      if (i>=nmat)
		{

		  /* New material */

		  /* Reallocate material name array */

		  matarr = (long*)Mem(MEM_REALLOC,matarr,  (nmat+1)*sizeof(long));

		  /* Put pointer to material name */

		  matarr[nmat] = PutText(mname);

		  nmat++;
		}

	      WDB[loc1 + IFC_TET_MSH_PTR_CELL] = (double)i;
	  
	      /* Next */
	      
	      loc1 = NextItem(loc1);
	    }

	  /* Close file */

	  fclose(fp);

	  /* Reset cell pointer */

	  cell = -1;

	  /* Create nmat cells */

	  for (i = 0; i < nmat; i++)
	    {
	      /* Create cell */

	      cell = NewItem(loc0 + IFC_PTR_GCELL_LIST, CELL_BLOCK_SIZE);

	      /* Allocate memory for cell collision counter */

	      AllocValuePair(cell + CELL_COL_COUNT);

	      /* Allocate memory for collision tet-cell */

	      ptr = AllocPrivateData(1, PRIVA_ARRAY);
	      WDB[cell + CELL_PTR_PREV_TET] = (double)ptr;

	      /* Put name for cell */
	      /* NB: GetText can be used only with ascii pointers */
	      /* in RDB[], that's why we have to first write the pointer */
	      /* from matarr[] to CELL_PTR_NAME */

	      WDB[cell + CELL_PTR_NAME] = matarr[i];
	      sprintf(name, "ifc cell %s", GetText(cell + CELL_PTR_NAME));
	      WDB[cell + CELL_PTR_NAME] = PutText(name);

	      /* Put material name */
	      /* Linked in processifctetmesh.c */

	      WDB[cell + CELL_PTR_MAT] = -matarr[i];

	    }

	  CloseList(cell);

	  /* Link geometry cells */

	  /* Loop over cells */

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	  while (loc1 > VALID_PTR)
	    {

	      /* Get material number from tet */

	      i = (long)RDB[loc1 + IFC_TET_MSH_PTR_CELL];

	      /* Get cell pointer from cell list */

	      cell = ListPtr(cell, i);

	      /* Put cell pointer to tet */

	      WDB[loc1 + IFC_TET_MSH_PTR_CELL] = (double)cell;

	      /* Next */
	      
	      loc1 = NextItem(loc1);
	    }


	  /* One of the materials is already known */
	  
	  nmat--;

	  /* Store the number of "extra" materials */
          /* This is needed for geometry plotter   */

	  WDB[DATA_OF_N_EXTRAMAT] = (double)nmat;

	  /* Free the material array*/

	  Mem(MEM_FREE,matarr);

	}
      else
	{
	  /* Put material name to cells */

	  cell = NewItem(loc0 + IFC_PTR_GCELL_LIST, CELL_BLOCK_SIZE);

	  /* Allocate memory for cell collision counter */

	  AllocValuePair(cell + CELL_COL_COUNT);

	  /* Allocate memory for collision tet-cell */

	  ptr = AllocPrivateData(1, PRIVA_ARRAY);
	  WDB[cell + CELL_PTR_PREV_TET] = (double)ptr;

	  /* Put name for cell */

	  WDB[cell + CELL_PTR_NAME] = PutText("ifc cell");

	  /* Put material name */
	  /* Linked in processifctetmesh.c */

	  WDB[cell + CELL_PTR_MAT] = -RDB[loc0 + IFC_PTR_MAT];

	  CloseList(cell);

	  /* Link geometry cells */

	  /* Loop over cells */

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

	  while (loc1 > VALID_PTR)
	    {

	      /* Put cell pointer to tet */

	      WDB[loc1 + IFC_TET_MSH_PTR_CELL] = (double)cell;

	      /* Next */
	      
	      loc1 = NextItem(loc1);
	    }


	}

    }

  /*******************************************************************/

  /***** Read density distribution ***********************************/

  /* Check if density distribution is given */

  if (strcmp(rfile, "-1"))
    {
      /* Open density file for reading */

      if ((fp1 = fopen(rfile, "r")) == NULL)
	Error(loc0, "Density file \"%s\" does not exist", rfile);

      /* Read header data */

      ReadOFHeader(fp1, &n, &i, (long *)dim);

      /* Get number of cells */

      nc = (long)RDB[loc0 + IFC_MSH_N_CELLS];

      /* Reset multiplier */

      mul = 1.0;

      /* Check units */

      if ((dim[0] == 0) && (dim[1] == 0))
	{
	  /* Relative values, multiply by nominal density */

	  mul = d0;
	}
      else if ((dim[0] == 1) && (dim[1] == -3))
	{
	  /* kg/m3, convert to g/cm3 */

	  mul = -0.001;
	}
      else
	Die(FUNCTION_NAME, "Undefined dimensions %ld/%ld", dim[0], dim[1]);

      /* Check type */
      /*
	if (n != OF_FILE_DENSITY)
	Note(loc0, "Possibly invalid format for density file");
      */
      /* Check size */

      if (i != nc)
	Error(loc0, "Invalid number of entries in density file");

      /* Reset maximum and minimum density */
	      
      dmax = 0.0;
      dmin = INFTY;
	      
      /* Loop over cells */

      /* Check if mesh has been divided and original cells */
      /* have been moved to parent cells */
	      
      if((loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PRNTS]) < VALID_PTR)
	loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

      while (loc1 > VALID_PTR)
	{
	  /* Read density */

	  line = ReadOFData(fp1, OF_FILE_DENSITY);

	  if (sscanf(line, "%lf", &d) == EOF)
	    Die(FUNCTION_NAME, "Not enough entries in density file");

	  /* Convert to g/cm3 */
		  
	  d = d*mul;

	  /* Compare to maximum and minimum */

	  if (fabs(d) > fabs(dmax))
	    dmax = d;
	  if (fabs(d) < fabs(dmin))
	    dmin = d;
		  
	  /* Put value */
		  
	  WDB[loc1 + IFC_TET_MSH_DF] = d;
		  
	  /* Next cell */
		  
	  loc1 = NextItem(loc1);
	}
	      
      /* Close file */
	      
      fclose(fp1);
	      
      /* Put maximum and minimum density            */
      /* Majorant is checked in processifctetmesh.c */

      WDB[loc0 + IFC_MAX_DENSITY] = dmax;
      /* WDB[loc0 + IFC_MIN_DENSITY] = dmin;*/

    }
  else
    {
      /* Put nominal values */

      WDB[loc0 + IFC_MAX_DENSITY] = d0;
      WDB[loc0 + IFC_MIN_DENSITY] = d0;

    }

  /*******************************************************************/

  /***** Read temperature distribution *******************************/

  /* Check if temperature distribution is given */

  if (strcmp(tfile, "-1"))
    {
      /* Open temperature file for reading */

      if ((fp1 = fopen(tfile, "r")) == NULL)
	Error(loc0, "Temperature file \"%s\" does not exist", tfile);

      /* Read header data */

      ReadOFHeader(fp1, &n, &i, (long *)dim);

      /* Reset multiplier */

      mul = 1.0;

      /* Get number of cells */

      nc = (long)RDB[loc0 + IFC_MSH_N_CELLS];

      /* Check units */

      if (dim[3] == 1)
	{
	  /* K */

	  mul = 1.0;
	}
      else if(dim[3] == 0)
	{
	  /* Header might be missing altogether (K) */

	  mul = 1.0;

	}
      else
	Die(FUNCTION_NAME, "Undefined dimensions %ld", dim[3]);

      /* Check type */
      /*
	if (n != OF_FILE_TEMP)
	Note(loc0, "Possibly invalid format for temperature file");
      */
      /* Check size */

      if (i != nc)
	Error(loc0, "Invalid number of entries in temperature file");
	      
      /* Reset maximum and minimum temperature */
	      
      Tmax = -1.0;
      Tmin = INFTY;
	      
      /* Loop over cells */

      /* Check if mesh has been divided and original cells */
      /* have been moved to parent cells */
	      
      if((loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PRNTS]) < VALID_PTR)
	loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	      
      while (loc1 > VALID_PTR)
	{
	  /* Read temperature */
		  
	  line = ReadOFData(fp1, OF_FILE_TEMP);

	  if (sscanf(line, "%lf", &T) == EOF)
	    Die(FUNCTION_NAME, "Not enough entries in temp. file");
		  
	  /* Check type flag */

	  if (ttype == 2)
	    T = T + T0;
	  else if (ttype != 1)
	    Die(FUNCTION_NAME, "Invalid temperature type flag");

	  /* Compare to maximum and minimum */
		  
	  if (T > Tmax)
	    Tmax = T;
	  if ((T > 0.0) && (T < Tmin))
	    Tmin = T;
		  
	  /* Put value */
		  
	  WDB[loc1 + IFC_TET_MSH_TMP] = T;

	  /* Next cell */
		  
	  loc1 = NextItem(loc1);
	}

      /* Close file */
	      
      fclose(fp1);
	      
      /* Put maximum and minimum temperature */

      if(!update)
	{
	  WDB[loc0 + IFC_MAX_TEMP] = Tmax;
	  WDB[loc0 + IFC_MIN_TEMP] = Tmin;
	      
	  /* Set TMS on */
	      
	  WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;
	}

    }
  else
    {
      /* Put nominal values */

      WDB[loc0 + IFC_MAX_TEMP] = T0;
      WDB[loc0 + IFC_MIN_TEMP] = T0;
    }

  /*******************************************************************/

  /***** Read output mapping file ************************************/

  if(!update)
    {	 

      /* Get number of cells */

      nc = (long)RDB[loc0 + IFC_MSH_N_CELLS];
 
      /* Check output flag */

      if ((long)RDB[loc0 + IFC_CALC_OUTPUT] == YES)
	{
	  /* Check if file is given */
	      
	  if (strcmp(mapfile, "-1"))
	    {
	      /* Open file for reading */

	      if ((fp1 = fopen(mapfile, "r")) == NULL)
		Error(loc0, "Output map file \"%s\" does not exist", 
		      mapfile);

	      /* Read header data */
		  
	      ReadOFHeader(fp1, &n, &i, (long *)dim);
		  
	      /* Check size */
		  
	      if (i != nc)
		Error(loc0, 
		      "Invalid number of entries in output map file");
		  
	      /* Loop over cells */
		  
	      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	      while (loc1 > VALID_PTR)
		{
		  /* Read index */
		      
		  line = ReadOFData(fp1, OF_FILE_MAP);
		      
		  if (sscanf(line, "%ld", &n) == EOF)
		    Die(FUNCTION_NAME, 
			"Not enough entries in output map file");
		  
		  /* Check value */
		      
		  if ((n < 1) || (n > nc))
		    Error(loc0, "Invalid entry %ld in output map file", n);
		  
		  /* Put value */
		      
		  WDB[loc1 + IFC_TET_MSH_STAT_IDX] = (double)(n - 1);
		      
		  /* Next cell */
		      
		  loc1 = NextItem(loc1);
		}
		  
	      /* Close file */
		  
	      fclose(fp1);
	    }
	  else
	    {
	      /* Map each cell to itself */
		  
	      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	      while (loc1 > VALID_PTR)
		{
		  /* Copy index */
		      
		  WDB[loc1 + IFC_TET_MSH_STAT_IDX] = 
		    RDB[loc1 + IFC_TET_MSH_IDX];
		      
		  /* Next cell */
		      
		  loc1 = NextItem(loc1);
		}
	    }
	}

    	  
      /*******************************************************************/

      /***** Error check and finalize ************************************/

      /* Get pointer to surface list */

      sslist = (long)RDB[loc0 + IFC_PTR_SURF_LIST];

      /* Loop over cells */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
      while (loc1 > VALID_PTR)
	{
	  
	  /* Check boundaries */
	  
	  if (RDB[loc1 + IFC_TET_MSH_XMIN] < RDB[loc0 + IFC_MESH_XMIN])
	    Die(FUNCTION_NAME, "Error in boundaries");
	  if (RDB[loc1 + IFC_TET_MSH_XMAX] > RDB[loc0 + IFC_MESH_XMAX])
	    Die(FUNCTION_NAME, "Error in boundaries");
	  if (RDB[loc1 + IFC_TET_MSH_YMIN] < RDB[loc0 + IFC_MESH_YMIN])
	    Die(FUNCTION_NAME, "Error in boundaries");
	  if (RDB[loc1 + IFC_TET_MSH_YMAX] > RDB[loc0 + IFC_MESH_YMAX])
	    Die(FUNCTION_NAME, "Error in boundaries");
	  if (RDB[loc1 + IFC_TET_MSH_ZMIN] < RDB[loc0 + IFC_MESH_ZMIN])
	    Die(FUNCTION_NAME, "Error in boundaries");
	  if (RDB[loc1 + IFC_TET_MSH_ZMAX] > RDB[loc0 + IFC_MESH_ZMAX])
	    Die(FUNCTION_NAME, "Error in boundaries");

	  /* Get number of faces */

	  nd = (long)RDB[loc1 + IFC_TET_MSH_NF];

	  /* Reset number of different face types */
	  
	  /* 3 point faces */
	  nfacepts[0] = 0;
	  /* 4 point faces */
	  nfacepts[1] = 0;
	  /* other faces */
	  nfacepts[2] = 0;

	  /* Get pointer to geometry cell */

	  cell = (long)RDB[loc1 + IFC_TET_MSH_PTR_CELL];

	  /* Get pointer to material */

	  ptr = (long)WDB[loc1 + IFC_TET_MSH_PTR_FACES];

	  ptr1 = (long)WDB[loc1 + IFC_TET_MSH_PTR_SIDES];

	  /* Loop over faces */

	  for (i = 0; i < nd ; i++)
	    {

	      /* Get surface index */

	      n = (long)RDB[ptr + i];

	      surf = ListPtr(sslist, n);

	      np = (long)RDB[surf + SURFACE_N_PARAMS];

	      /* Store number of points */

	      if (np == 3)
		nfacepts[0]++;
	      else if (np == 4)
		nfacepts[1]++;
	      else
		nfacepts[2]++;

	    }

	  /* Check if polyhedral cell */

	  /* Less than four or more than six faces */
	  if ((nd < 4) || (nd > 6))
	    npolyc++;
	  /* Contains faces that have less than three or more than */
	  /* four points */
	  else if (nfacepts[2] > 0)
	    npolyc++;
	  /* Four sides with three points in each */
	  else if ((nd == 4) && (nfacepts[0] == 4))
	    ntetc++;
	  /* Five sides with three points in 4, four points in 1 */
	  else if ((nd == 5) && ((nfacepts[0] == 4) && (nfacepts[1] == 1)))
	    npyramc++;
	  /* Five sides with three points in 2, four points in 3 */
	  else if ((nd == 5) && ((nfacepts[0] == 2) && (nfacepts[1] == 3)))
	    nprismc++;
	  /* Six sides with four points in each */
	  else if ((nd == 6) && (nfacepts[1] == 6))
	    nhexc++;
	  else
	    npolyc++;

	  /* Next */
	  
	  loc1 = NextItem(loc1);
	}

      /* Read batches */

      ReadOFBatches(loc0);

      /* Print out number of different cells */
      fprintf(out, "Composition of mesh:\n");
      fprintf(out, "Tetrahedrons:      %ld\n", ntetc);
      fprintf(out, "Pyramids:          %ld\n", npyramc);
      fprintf(out, "Prisms:            %ld\n", nprismc);
      fprintf(out, "Hexahedrons:       %ld\n", nhexc);
      fprintf(out, "Other polyhedrons: %ld\n", npolyc);

      /*****************************************************/
      /* Calculate cell centerpoints */

      /* Get number of cells */

      nc = (long)RDB[loc0 + IFC_NC];

      /* Allocate memory for cell centerpoints */

      cellpts = ReallocMem(DATA_ARRAY, nc*3);

      /* Store cell centerpoint list */
      
      WDB[loc0 + IFC_PTR_PRNT_CELL_CP_LIST] = (double)cellpts;

      /* Get pointer to tet list */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);  

      /* Get pointer to surface list */

      surflist = (long)RDB[loc0 + IFC_PTR_SURF_LIST];
      CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

      /* Calculate cell centerpoints */
      
      CalculateTetCenter(loc1, surflist, cellpts, -1);

      /*****************************************************/

      /* Check initial mesh */

      fprintf(out, "\nChecking initial mesh:\n\n");
      CheckPolyhedMesh(loc0);

      /* Check cell divisor flag and convert to tet mesh */

      if (divide == YES)
	{

	  /* Move tet-list to parents and reset tet-list */

	  WDB[loc0 + IFC_PTR_TET_MSH_PRNTS] = (double)loc1;
	  WDB[loc0 + IFC_PTR_TET_MSH] = (double)NULLPTR;

	  /* Move surface list to parents and reset main list */

	  loc1 = (long)RDB[loc0 + IFC_PTR_SURF_LIST];

	  WDB[loc0 + IFC_PTR_SURF_LIST_PRNTS] = (double)loc1;
	  WDB[loc0 + IFC_PTR_SURF_LIST] = (double)NULLPTR;

	  /* Move owner list to parents and reset main list */

	  loc1 = (long)RDB[loc0 + IFC_PTR_OWNR_LIST];

	  WDB[loc0 + IFC_PTR_OWNR_LIST_PRNTS] = (double)loc1;
	  WDB[loc0 + IFC_PTR_OWNR_LIST] = (double)NULLPTR;

	  /* Move neighbour list to parents and reset main list */

	  loc1 = (long)RDB[loc0 + IFC_PTR_NBR_LIST];

	  WDB[loc0 + IFC_PTR_NBR_LIST_PRNTS] = (double)loc1;
	  WDB[loc0 + IFC_PTR_NBR_LIST] = (double)NULLPTR;

	  /* Move point list to parents and reset main list */

	  loc1 = (long)RDB[loc0 + IFC_PTR_POINT_LIST];

	  WDB[loc0 + IFC_PTR_POINT_LIST_PRNTS] = (double)loc1;
	  WDB[loc0 + IFC_PTR_POINT_LIST] = (double)NULLPTR;

	  /* Move cell centerpoint list to parents and reset main list */

	  loc1 = (long)RDB[loc0 + IFC_PTR_CELL_CP_LIST];

	  WDB[loc0 + IFC_PTR_PRNT_CELL_CP_LIST] = (double)loc1;
	  WDB[loc0 + IFC_PTR_CELL_CP_LIST] = (double)NULLPTR;

	  /* Get total number of faces */

	  WDB[loc0 + IFC_NF_PRNTS] = RDB[loc0 + IFC_NF];
	  WDB[loc0 + IFC_NF] = 0;

	  /* Get total number of cells */

	  WDB[loc0 + IFC_NC_PRNTS] = RDB[loc0 + IFC_NC];
	  WDB[loc0 + IFC_NC] = 0;

	  if (npolyc == 0)
	    {
	      /******* Calculate number of new faces ********************/
	      /* Each hex cell will have 6*2 external faces and 6 internal faces        */
	      /* Each prism cell will have 3*2+2 external faces and 2(?) internal faces */
	      /* Each pyramid cell will have 6 external faces and 1 internal face       */
	      /* Each tet cell will have 4 external faces                               */
	      /* NB; Some hex cells only have 4 internal faces                          */

	      nf = nhexc*(6*2 + 6) + nprismc*(3*2 + 2 + 2) + npyramc*(6+1) +ntetc*4;

	      fprintf(out, "\nDividing hexahedral mesh to tetrahedrons.\n");
	      FixHexMesh(loc0, nf);
	    }
	  else
	    {
	      fprintf(out, "\nDividing polyhedral mesh to tetrahedrons.\n");
	      FixPolyhedMesh(loc0);
	    }
	  
	  fprintf(out, "Checking divided mesh.\n");
	  CheckPolyhedMesh(loc0);

	  /* Store number of new cells */

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  WDB[loc0 + IFC_NC] = (double)ListSize(loc1);

	}
      else
	{

	  /* Calculate tet-cell centerpoints */

	  /* Store number of new cells */

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  /* Get number of cells */

	  nc = ListSize(loc1);

	  /* Put number of cells */

	  WDB[loc0 + IFC_NC] = (double)nc;

	  /* Allocate memory for cell centerpoints */

	  cellpts = ReallocMem(DATA_ARRAY, nc*3);

	  /* Store cell centerpoint list */

	  WDB[loc0 + IFC_PTR_CELL_CP_LIST] = (double)cellpts;

	  /* Get pointer to tet list */

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  /* Get pointer to surface list */

	  surflist = (long)RDB[loc0 + IFC_PTR_SURF_LIST];
	  CheckPointer(FUNCTION_NAME, "(surflist)", DATA_ARRAY, surflist);

	  /* Calculate cell centerpoints for new cells */

	  CalculateTetCenter(loc1, surflist, cellpts, -1);

	}

    }
  
  /* Switch type to regular tet mesh */

  WDB[loc0 + IFC_TYPE] = (double)IFC_TYPE_TET_MESH;

  /* Copy temperature and density data to children when updating */

  /* Loop over cells */

  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
  while (loc1 > VALID_PTR)
    {
      /* Get pointer to parent */

      ptr = RDB[loc1 + IFC_TET_MSH_PTR_PARENT];

      if(ptr > VALID_PTR)
	{
	  /* Copy density from parent */

	  WDB[loc1 + IFC_TET_MSH_DF] = RDB[ptr + IFC_TET_MSH_DF];

	  /* Copy temperature from parent */

	  WDB[loc1 + IFC_TET_MSH_TMP] = RDB[ptr + IFC_TET_MSH_TMP];

	  /* Copy statistics bin from parent */

	  WDB[loc1 + IFC_TET_MSH_STAT_IDX] = RDB[ptr + IFC_TET_MSH_STAT_IDX];
	}     
      /* Next */

      loc1 = NextItem(loc1);
    }
  
    
  fprintf(out, "\n");
  /*******************************************************************/

}

/*****************************************************************************/
