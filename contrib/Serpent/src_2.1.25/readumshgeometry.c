/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readumshgeometry.c                             */
/*                                                                           */
/* Created:       2013/11/23 (JLe)                                           */
/* Last modified: 2015/09/15 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Reads unstructured mesh based geometry                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadUMSHGeometry:"

/*****************************************************************************/

void ReadUMSHGeometry()
{  
  long umsh, loc0, loc1, ptr, ptr1, nd, nc, nf, np, n, i, j, k, nmax, surf;
  long cell, nbhr, fi, *matarr, p, dim[7], nmat, cellpts, surflist;
  long solist, snlist, sslist, pointlist, **surfaces, **sides, *nsurf;
  long divide;
  long nhexc, ntetc, npyramc, nprismc, npolyc;
  long nfacepts[3];
  double x, y, z, xmin, xmax, ymin, ymax, zmin, zmax;
  char mname[MAX_STR], tmpstr[MAX_STR], *line;
  FILE *fp;
#ifdef printoutumsh
  long mat, loc2, loc3;
#endif

  /* Check pointer */

  if ((umsh = (long)RDB[DATA_PTR_UMSH0]) < VALID_PTR)
    return;
  
  fprintf(out, "Reading unstructured mesh based geometries...\n");
  
  /* Loop over definitions */
  
  while (umsh > VALID_PTR)
    {

      if((long)RDB[umsh + UMSH_PTR_FNAME] > 0)
	{
	  /* Will be read in readifcofmesh.c */
	  
	  umsh = NextItem(umsh);

	  continue;
	}

      /* Pointer to interface structure */

      loc0 = (long)RDB[umsh + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Reset divide flag */

      divide = NO;

      /* Reset number of different cell types */

      nhexc = 0;
      ntetc = 0;
      nprismc = 0;
      npolyc = 0;      
      npyramc = 0;

      /***********************************************************************/

      /***** Read points *****************************************************/

      /* Reset mesh boundaries */
	  
      xmin =  INFTY;
      xmax = -INFTY;
      ymin =  INFTY;
      ymax = -INFTY;
      zmin =  INFTY;
      zmax = -INFTY;

      /* Check file format */

      TestDOSFile(GetText(umsh + UMSH_PTR_POINTS_FNAME));

      /* Open points file for reading */
      
      if ((fp = fopen(GetText(umsh + UMSH_PTR_POINTS_FNAME), "r")) == NULL)
	Error(umsh, "Points file \"%s\" does not exist", 
	      GetText(umsh + UMSH_PTR_POINTS_FNAME));

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
	    Error(umsh, "Not enough entries in points file");
	  
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

      /***********************************************************************/

      /***** Read faces ******************************************************/
      
      /* Check file format */

      TestDOSFile(GetText(umsh + UMSH_PTR_FACES_FNAME));

      /* Open faces file for reading */
      
      if ((fp = fopen(GetText(umsh + UMSH_PTR_FACES_FNAME), "r")) == NULL)
	Error(umsh, "Faces file \"%s\" does not exist", 
	      GetText(umsh + UMSH_PTR_FACES_FNAME));

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
	    Error(umsh, "Not enough points");
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

      /* Loop over faces and print out */
#ifdef printoutumsh
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
      /************************************************************************/

      /***** Read owner and neighbour files **********************************/

      /* Check file formats */

      TestDOSFile(GetText(umsh + UMSH_PTR_OWNER_FNAME));
      TestDOSFile(GetText(umsh + UMSH_PTR_NEIGHBOUR_FNAME));

      /* Reset number of cells */

      nc = -1;

      /* Repeat for owner and neighbour files */
	  
      for (fi = 1; fi < 3; fi++)
	{
	  /* Check mode */
	  
	  if (fi == 1)
	    {
	      /* Open file */

	      if ((fp = fopen(GetText(umsh + UMSH_PTR_OWNER_FNAME), "r")) 
		  == NULL)
		Error(umsh, "Owner file \"%s\" does not exist", 
		      GetText(umsh + UMSH_PTR_OWNER_FNAME));
	      
	      /* Read header data */
	      
	      ReadOFHeader(fp, &n, &nmax, (long *)dim);

	    }
	  else
	    {
	      /* Open file */

	      if ((fp = fopen(GetText(umsh + UMSH_PTR_NEIGHBOUR_FNAME), "r")) 
		  == NULL)
		Error(umsh, "Neighbour file \"%s\" does not exist", 
		      GetText(umsh + UMSH_PTR_NEIGHBOUR_FNAME));
	      
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
		    Error(umsh, "Not enough entries in owner file");
		  else
		    Error(umsh, "Not enough entries in neighbour file");
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

      /***********************************************************************/

      /***** Put geometry data to cells **************************************/

      /* Get pointer to surface list */

      sslist = (long)RDB[loc0 + IFC_PTR_SURF_LIST];

      /* Check number of faces */

      CheckValue(FUNCTION_NAME, "nf", "", nf, 4, 100000000000);

      /* Put number of faces */

      WDB[loc0 + IFC_NF] = (double)nf;

      /* Put number of cells */

      WDB[loc0 + IFC_NC] = (double)nc;

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
		Error(umsh, "Cell index exceeds maximum");
	      
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
			Error(umsh, "Cell index exceeds maximum");
		  
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

	      /* First, check that cells surface list is not full */

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

      /***********************************************************************/

      /***** Read materials **************************************************/

      /* Check file format */

      TestDOSFile(GetText(umsh + UMSH_PTR_MATERIALS_FNAME));

      /* Open materials file */
      
      if ((fp = fopen(GetText(umsh + UMSH_PTR_MATERIALS_FNAME), "r")) == NULL)
	Error(umsh, "Materials file \"%s\" does not exist", 
	      GetText(umsh + UMSH_PTR_MATERIALS_FNAME));      

      /* Read header data */

      ReadOFHeader(fp, &n, &i, (long *)dim);

      /* Check size */

      if (i != nc)
	Error(umsh, "Invalid number of entries in material file");

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
	  
	  for (i = 0; i < nmat; i++)
	    if(!strcmp(&ASCII[matarr[i]],mname))
	      break;

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

	  WDB[cell + CELL_PTR_NAME] = PutText("ifc cell");

	  /* Put material name */
	  /* Linked in processumshgeometry.c */

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

      /* Free the material array*/

      Mem(MEM_FREE,matarr);

      /***********************************************************************/
      
      /***** Error check and finalize ****************************************/

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

 #ifdef printoutumsh

	  /* Get pointer to geometry cell */

	  cell = (long)RDB[loc1 + IFC_TET_MSH_PTR_CELL];

	  /* Get pointer to material */
	  /* NO: Material names are replaced with material pointers only at processumshgeom... */

	  mat = (long)RDB[cell + CELL_PTR_MAT];

	  printf("****\nCell %ld, material %s, %ld faces\n", 
		 (long)RDB[loc1 + IFC_TET_MSH_IDX],
		 &ASCII[-mat], nd);

#endif

	  ptr = (long)WDB[loc1 + IFC_TET_MSH_PTR_FACES];

	  ptr1 = (long)WDB[loc1 + IFC_TET_MSH_PTR_SIDES];

	  /* Loop over faces */

	  for (i = 0; i < nd ; i++)
	    {
	      /*
	      printf("Face %ld, surf %ld, side %ld\n", i, (long)RDB[ptr + i], (long)RDB[ptr1 + i]);
	      */
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

 #ifdef printoutumsh

	      printf("%ld points: ",np);

	      /* Get pointer to surface parameters (point list) */
	      
	      loc2 = (long)RDB[surf + SURFACE_PTR_PARAMS];

	      /* Loop over points */

	      for (j = 0; j < np; j++)
		{
		  /* Get point */
		  loc3 = (long)RDB[loc2 + j];

		  /* Print point */ 
		  printf("(%f %f %f) ", RDB[loc3], RDB[loc3 + 1], RDB[loc3 + 2]);
		}

	      printf("\n");

#endif
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

      /* Allocate memory for previous cell pointer */

      ptr = AllocPrivateData(1, PRIVA_ARRAY);
      WDB[loc0 + IFC_PTR_PREV_CELL] = (double)ptr;

      /* Allocate memory for previous collison */

      AllocValuePair(loc0 + IFC_PTR_PREV_COL_CELL);

      /* Store size */

      WDB[umsh + UMSH_N_ORIG_CELLS] = RDB[loc0 + IFC_NC];
      WDB[umsh + UMSH_N_CELLS] = RDB[loc0 + IFC_NC];

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

	  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

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

	  /* Put initial number of cells */

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

	      fprintf(out, "\nDividing simple mesh to tetrahedrons.\n");
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
	  WDB[umsh + UMSH_N_CELLS] = (double)ListSize(loc1);

	  /* Copy statistics bins to children */

	  /* Loop over cells */

	  while (loc1 > VALID_PTR)
	    {
	      /* Get pointer to parent */

	      ptr = RDB[loc1 + IFC_TET_MSH_PTR_PARENT];

	      if(ptr > VALID_PTR)
		{
		  /* Copy statistics bin from parent */

		  WDB[loc1 + IFC_TET_MSH_STAT_IDX] = RDB[ptr + IFC_TET_MSH_STAT_IDX];
		}     
	      /* Next */

	      loc1 = NextItem(loc1);
	    }

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

      /***********************************************************************/

      /* Next geometry */

      umsh = NextItem(umsh);

      /***********************************************************************/
    }

  /* Done */

  fprintf(out, "OK.\n\n");

}

/*****************************************************************************/
