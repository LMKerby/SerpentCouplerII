/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifctetmesh.c                            */
/*                                                                           */
/* Created:       2015/01/15 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes tetra mesh based multi-physics interfaces          */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*           - Split from processinterface.c                                 */
/*           - TODO: Loop over tet cells to set materials and densities not  */
/*             working                                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCTetMesh:"

/*****************************************************************************/

void ProcessIFCTetMesh(long loc0, long update)
{
  long loc1, loc2, mat0, mat, ptr, found, msh, np;
  long cell, count, nc, ncc, cell0, override;
  char *mname;
  double xmin, xmax, ymin, ymax, zmin, zmax, dmax, Tmax, Tmin;
  double last, frac, matdens, d, T;
  double lims[6];

  /***********************************************************************/

  mat = -1;

  matdens = 0.0;

  override = 0;

  /***** Link materials ***************************************************/

  /* Tet mesh with a single material */

  if(RDB[loc0 + IFC_PTR_OF_MFILE] < VALID_PTR) 
    {

      if(!update)
	{
	  /* Reset found flag */

	  found = NO;
	  
	  /* Reset pointer */
	      
	  mat = -1;
	      
	  /* Loop over materials and find match */
	  
	  mat0 = (long)RDB[DATA_PTR_M0];
	  while (mat0 > VALID_PTR)
	    {
	      /* Compare name */
	      
	      if (CompareStr(mat0 + MATERIAL_PTR_NAME, loc0 + IFC_PTR_MAT))
		mat = mat0;
	      
	      /* Check if material was divided for burnup calculation */
	      
	      if ((ptr = (long)RDB[mat0 + MATERIAL_DIV_PTR_PARENT]) 
		  > VALID_PTR)
		if (CompareStr(ptr + MATERIAL_PTR_NAME, loc0 + IFC_PTR_MAT))
		  mat = mat0;
	      
	      /* Check material */
	      
	      if (mat > VALID_PTR)
		{
		  /* Set found flag */
		  
		  found = YES;
		}
	      
	      /* Next material */
	      
	      mat0 = NextItem(mat0);
	    }
      
	  /* Check that material was found */
      
	  if (found == NO)
	    Error(loc0, "Material %s in distribution file not defined", 
		  GetText(loc0 + IFC_PTR_MAT));

	  /* Link interface and put flag */
		  
	  WDB[mat + MATERIAL_PTR_IFC] = (double)loc0;
	  WDB[mat + MATERIAL_USE_IFC] = (double)YES;		  

	  WDB[loc0 + IFC_PTR_MAT] = (double)mat;
	}

    }


  /* Loop over tet cells to link materials */
  
  /* If cell has been divided, only do this to parents */

  if((loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PRNTS]) < VALID_PTR)
    loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

  /* Reset maximum and minimum temperature */
	      
  Tmax = -1.0;
  Tmin = INFTY;

  /* Reset previous material */

  mat0 = -1;

  while (loc1 > VALID_PTR)
    {
      /* Get geometry cell */

      cell = (long)RDB[loc1 + IFC_TET_MSH_PTR_CELL];
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

      /* Get material from cell or link material to cell */

      if(!update)
	{
	  if ((RDB[cell + CELL_PTR_MAT]) < 0)
	    {

	      /* Not yet linked */

	      if (RDB[loc0 + IFC_PTR_OF_MFILE] < VALID_PTR) 
		{
		  /* Tet mesh with a single material */
		  
		  /* Link material from interface */
		  WDB[cell + CELL_PTR_MAT] = RDB[loc0 + IFC_PTR_MAT];

		}
	      else
		{
		  /* Tet mesh with multiple materials */

		  /* Remove minus sign */

		  WDB[cell + CELL_PTR_MAT] = -RDB[cell + CELL_PTR_MAT];

		  /* Get material name */

		  mname = GetText(cell + CELL_PTR_MAT);

		  /* Check previous */

		  mat = -1;

		  if (mat0 > VALID_PTR)
		    if (!strcmp(mname, GetText(mat0 + MATERIAL_PTR_NAME)))
		      mat = mat0;

		  /* Find match */
	  
		  if (mat < VALID_PTR)
		    {
		      mat = (long)RDB[DATA_PTR_M0];
		      while (mat > VALID_PTR)
			{
			  /* Compare name */
			  
			  if (!strcmp(mname, GetText(mat + MATERIAL_PTR_NAME)))
			    break;
			  
			  /* Next material */
			  
			  mat = NextItem(mat);
			}
		    }
		  
		  /* Check match */
	  
		  if (mat < VALID_PTR)
		    Error(loc0, "Material %s is not defined", mname);

		  /* Set pointer */

		  WDB[cell + CELL_PTR_MAT] = mat;
		  mat0 = mat;

		}
	    }
	  else
	    {
	      /* Already linked (OF_SOLID) or linked */
	      /* during this run*/

	      mat = RDB[cell + CELL_PTR_MAT];

	    }
	}
      else
	mat = RDB[cell + CELL_PTR_MAT];

      /* Get tet-cell density */

      d = WDB[loc1 + IFC_TET_MSH_DF];

      /* Get maximum density */

      if ((dmax = RDB[loc0 + IFC_MAX_DENSITY]) == 0.0)
	if (strcmp(GetText(loc0 + IFC_PTR_OF_RFILE), "-1"))
	  Die(FUNCTION_NAME, "Zero maximum density");

      /* Check that density format corresponds with material definition */
      if (!update)
	{

	  /* Reset override flag */

	  override = 0;
	  
	  /* If no density file is given, we don't have to check this */

	  if (strcmp(GetText(loc0 + IFC_PTR_OF_RFILE), "-1"))
	    if (RDB[mat + MATERIAL_ADENS]*d < 0)
	      override = 1;
	}   

      /* Check material maximum density */

      if(update)
	{
	  /* When the interface is updated, the MATERIAL_ADENS */
	  /* Already contains the atomic density */

	  if (dmax < 0)
	    /* Mass densities given*/
	    matdens = -RDB[mat + MATERIAL_MDENS];	  
	  else if(dmax > 0)
	    /* Atomic densities given */
	    matdens = RDB[mat + MATERIAL_ADENS];
	  else
	    {
	      matdens = 0;
	      Error(loc0, "No non-zero densities in distribution"); 
	    }
	}
      else
	{
	  /* If this is the first time this interface is processed we will use */
	  /* MATERIAL_ADENS for the density regardless of IFC units */
	  /* They will be converted to atomic density later         */

	  matdens = RDB[mat + MATERIAL_ADENS];
	}

      /* Only increase the material atomic density before calculating majorants */
      /* (not on updates) */

      /* Only change it upwards, this way the normal density card can be used as */
      /* a maximum density, when setting up a coupled calculation */

      if (strcmp(GetText(loc0 + IFC_PTR_OF_RFILE), "-1"))
	{

	  if ((fabs(matdens) < fabs(d)) || (override == 1))
	    {
	      /* IFC density greater than material base density */

	      if(!update)
		{
		  WDB[mat + MATERIAL_ADENS] = dmax;
		  matdens = RDB[mat + MATERIAL_ADENS];
		}
	    }
	  else
	    {
	      /* Material base density greater than IFC density */

	      dmax = matdens;

	    }
	}

      /* Get tet-cell temperature */

      T = WDB[loc1 + IFC_TET_MSH_TMP];

      /* Compare to maximum and minimum */
		  
      if (T > Tmax)
	Tmax = T;
      if ((T > 0.0) && (T < Tmin))
	Tmin = T;

      /* Compare to material maximum and minimum */

      if (Tmax > RDB[mat + MATERIAL_TMS_TMAX])
	{
	  if(!update)
	    WDB[mat + MATERIAL_TMS_TMAX] = Tmax;
	  else
	    Die(FUNCTION_NAME,
		"Material temperature above TMS majorant for material %s",
		GetText(mat + MATERIAL_PTR_NAME));
	}

      /* Put minimum temperature */
		      
      if (Tmin < RDB[mat + MATERIAL_TMS_TMIN])
	{
	  if(!update)
	    WDB[mat + MATERIAL_TMS_TMIN] = Tmin;
	  else
	    Die(FUNCTION_NAME,
		"Material temperature below TMS minorant for material %s",
		GetText(mat + MATERIAL_PTR_NAME));

	}
		      
      if(!update)
	{

	  /* Set on-the-fly Doppler-broadening mode */

	  if (RDB[mat + MATERIAL_TMS_TMAX] !=
	      RDB[mat + MATERIAL_TMS_TMIN] ) 
	    WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

	  /* Link interface and put flag */
		  
	  WDB[mat + MATERIAL_PTR_IFC] = (double)loc0;
	  WDB[mat + MATERIAL_USE_IFC] = (double)YES;		  

	}

      loc1 = NextItem(loc1);
    }

  /* Copy cell pointers from parent tets to children */

  if(!update)
    {
      /* Get pointer to first tet cell */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

      /* Loop over tet cells */

      while(loc1 > VALID_PTR)
	{

	  /* Get pointer to parent */
	  
	  if ((loc2 = (long)RDB[loc1 + IFC_TET_MSH_PTR_PARENT]) < VALID_PTR)
	    {
	      /* Pointer to next */

	      loc1 = NextItem(loc1);

	      /* Cycle loop */

	      continue;
	    }
	  
	  /* Get cell pointer from parent */

	  cell0 = (long)RDB[loc2 + IFC_TET_MSH_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell0)", DATA_ARRAY, cell0);

	  /* Copy just the pointer to cell */

	  WDB[loc1 + IFC_TET_MSH_PTR_CELL] = (double)cell0;

	  /* Next tet */

	  loc1 = NextItem(loc1);
	}
    }


  /***** Unstructured tetrahedral mesh ********************************/

  /* Get limits */

  xmin = RDB[loc0 + IFC_MESH_XMIN];
  xmax = RDB[loc0 + IFC_MESH_XMAX];
  ymin = RDB[loc0 + IFC_MESH_YMIN];
  ymax = RDB[loc0 + IFC_MESH_YMAX];
  zmin = RDB[loc0 + IFC_MESH_ZMIN];
  zmax = RDB[loc0 + IFC_MESH_ZMAX];

  fprintf(out, "\nOpenFOAM mesh based interface \"%s\":\n\n",
	  GetText(loc0 + IFC_PTR_INPUT_FNAME));

  fprintf(out, " - Dimensions: x = [%1.1f, %1.1f]; ", 
	  xmin, xmax);
  fprintf(out, "y = [%1.1f, %1.1f]; ", 
	  ymin, ymax);
  fprintf(out, "z = [%1.1f, %1.1f]\n",
	  zmin, zmax);

  /* Get number of cells */

  nc = (long)RDB[loc0 + IFC_NC_PRNTS];

  /* Get number of child cells */

  ncc = (long)RDB[loc0 + IFC_NC];

  /* Print number of cells */
  if(nc != 0)
    {
      fprintf(out, " - Initial number of cells: %ld\n", nc);

      /* Print number of child cells */

      fprintf(out, " - Divided into a number of child cells: %ld\n", ncc);

    }
  else
    fprintf(out, " - Number of cells: %ld\n", ncc);

  /* Get maximum density */

  if ((dmax = RDB[loc0 + IFC_MAX_DENSITY]) == 0.0)
    Die(FUNCTION_NAME, "Zero maximum density");

  /* Print IFC-density limits */

  if (strcmp(GetText(loc0 + IFC_PTR_OF_RFILE), "-1"))
    {
      fprintf(out, " - Maximum density of interface: %E\n",
	      RDB[loc0 + IFC_MAX_DENSITY]);
    }
  else
    {
      if (RDB[loc0 + IFC_PTR_OF_MFILE] < VALID_PTR)
	fprintf(out, " - Using mat-card material density: %E\n",
		matdens);
      else
	fprintf(out, " - Using mat-card material densities.\n");
    }

  /* Print maximum density and TMS limits of linked material */

  if(RDB[loc0 + IFC_PTR_OF_MFILE] < VALID_PTR)
    {

      mat = (long)RDB[loc0 + IFC_PTR_MAT];

      fprintf(out, " - Majorant density of material %s: %E\n", 
	     GetText(mat + MATERIAL_PTR_NAME), matdens);

      /* Print material temperature TMS limits */

      if (RDB[mat + MATERIAL_TMS_MODE] == (double)NO)
	fprintf(out, " - No TMS treatment for material %s\n", 
		GetText(mat + MATERIAL_PTR_NAME));
      else
	fprintf(out, " - Material %s TMS limits: %.2f K and %.2f K\n", 
		GetText(mat + MATERIAL_PTR_NAME), 
		RDB[mat + MATERIAL_TMS_TMIN], RDB[mat + MATERIAL_TMS_TMAX]);
    }
  else
    {
      fprintf(out, " - Possibly multiple materials in interface\n");
    }

  /* Create search mesh */

  if(!update)
    {

      /* Check boundaries */

      if ((xmin >= xmax) || (ymin >= ymax) || (zmin >= zmax))
	Error(loc0, "Structure is not 3D");

      /* Adjust boundaries */
      
      xmin = xmin - 1E-6;
      xmax = xmax + 1E-6;
      ymin = ymin - 1E-6;
      ymax = ymax + 1E-6;
      zmin = zmin - 1E-6;
      zmax = zmax + 1E-6;

      /* Put mesh variables */

      lims[0] = xmin;
      lims[1] = xmax;
      lims[2] = ymin;
      lims[3] = ymax;
      lims[4] = zmin;
      lims[5] = zmax;

      /* Read mesh split criterion */

      np = (long)RDB[loc0 + IFC_SEARCH_MESH_ADA_SPLIT];
      CheckValue(FUNCTION_NAME, "np", "", np, 1, INFTY);

      /* Get pointer to size vector */

      ptr = (long)RDB[loc0 + IFC_SEARCH_MESH_ADA_PTR_SZ];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Create mesh structure */
	  
      msh = CreateMesh(MESH_TYPE_ADAPTIVE, MESH_CONTENT_PTR, 
		       MESH_CONTENT_DATA_TET, np, 0, 0, lims, ptr);
	  
      /* Put pointer */
      
      WDB[loc0 + IFC_PTR_SEARCH_MESH] = (double)msh;

      /* Print progress */

      fprintf(out, "\nCreating search mesh for interface %s:\n\n",
	      GetText(loc0 + IFC_PTR_INPUT_FNAME));

      last =  0.0;
      count = 0; 
	  
      /* Loop over tet cells to create search mesh */

      loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

      while (loc1 > VALID_PTR)
	{
	  /* Preallocate memory for counters */

	  if ((long)RDB[DATA_REAL_PRIVA_SIZE] - 
	      (long)RDB[DATA_ALLOC_PRIVA_SIZE] < 100)
	    PreallocMem(100, PRIVA_ARRAY);

	  /* Calculate fraction and print progress */
	      
	  frac = (double)(count++)/((double)ListSize(loc1));

	  if (frac - last > 0.10)
	    {
	      fprintf(out, " %3.0f%% complete\n", 100.0*frac);
	      last = frac;
	    }

	  /* Get limits */

	  xmin = RDB[loc1 + IFC_TET_MSH_XMIN];
	  xmax = RDB[loc1 + IFC_TET_MSH_XMAX];
	  ymin = RDB[loc1 + IFC_TET_MSH_YMIN];
	  ymax = RDB[loc1 + IFC_TET_MSH_YMAX];
	  zmin = RDB[loc1 + IFC_TET_MSH_ZMIN];
	  zmax = RDB[loc1 + IFC_TET_MSH_ZMAX];

	  /* Add to search mesh */

	  AddSearchMesh(msh, loc1, xmin, xmax, ymin, ymax, zmin, zmax);

	  /* Next */

	  loc1 = NextItem(loc1);
	}

    }

  fprintf(out, " %3.0f%% complete\n\n", 100.0);

  fprintf(out, "Processing densities for interface %s...\n",
	  GetText(loc0 + IFC_PTR_INPUT_FNAME));


  /* Loop over tet cells to set densities */

  loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
  while (loc1 > VALID_PTR)
    {
	    
      if(RDB[loc0 + IFC_PTR_OF_MFILE] > VALID_PTR)
	{
	  /* Get geometry cell */

	  cell = (long)RDB[loc1 + IFC_TET_MSH_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
	  /* Get material */

	  mat = RDB[cell + CELL_PTR_MAT];
	  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

	  /* Set material densities */
	  /* Only increase the material atomic density before calculating majorants */
	  /* (not on updates) */

	  /* Only change it upwards, this way the normal density card can be used as */
	  /* a maximum density, when setting up a coupled calculation */

	  if(update)
	    {
	      /* Atomic densities given in interface */
	      
	      if(RDB[mat + MATERIAL_MAX_ADENS] > 0)
		/* Atomic densities given */
		matdens = RDB[mat + MATERIAL_ADENS];
	      else
		/* Mass densities given */
		matdens = -RDB[mat + MATERIAL_MDENS];

	    }
	  else
	    {
	      /* If this is the first time this interface is processed we will use */
	      /* MATERIAL_ADENS for the density regardless of IFC units */
	      /* They will be converted to atomic density later         */

	      matdens = RDB[mat + MATERIAL_ADENS];
	    }

	    
	  if(fabs(matdens) <= fabs(RDB[mat + MATERIAL_MAX_ADENS]))
	    {
	      if(!update)
		WDB[mat + MATERIAL_ADENS] =  RDB[mat + MATERIAL_MAX_ADENS];       		  
		
	    }
	  else
	    {
	      /* Store majorant density */
	      if(!update)
		WDB[mat + MATERIAL_MAX_ADENS] = matdens;
	      else
		Die(FUNCTION_NAME, "Material density has increased over the majorant?");
	    }
	  
	  /* Store density factor */

	  if(fabs(WDB[mat + MATERIAL_MAX_ADENS]) > ZERO)
	    {
	      /* Calculate factor of material density relative to majorant density */

	      WDB[loc1 + IFC_TET_MSH_DF] = RDB[loc1 + IFC_TET_MSH_DF]/RDB[mat + MATERIAL_MAX_ADENS];

	      /* If density file was not used set density factor to 1.0 */

	      if( !strcmp(GetText(loc0 + IFC_PTR_OF_RFILE), "-1"))
		WDB[loc1 + IFC_TET_MSH_DF] = 1.0;

	      /* Check that the density factor is below unity*/

	      if(WDB[loc1 + IFC_TET_MSH_DF] > 1.0)
		Warn(FUNCTION_NAME, "Larger than unity density factor for interface %s: %E",
		     GetText(loc0 + IFC_PTR_INPUT_FNAME), WDB[loc1 + IFC_TET_MSH_DF]);

	    }
	  else
	    /* Zero material density */
	    WDB[loc1 + IFC_TET_MSH_DF] = 1.0;

	}
      else
	{

	  /* Convert densities to density factors */

	  if (RDB[loc1 + IFC_TET_MSH_DF]*dmax < 0.0)
	    Error(loc0, "Inconsistent densities given in distribution %E %E",
		  RDB[loc1 + IFC_TET_MSH_DF],dmax);
	  else
	    {
	      /* Store density factor */

	      WDB[loc1 + IFC_TET_MSH_DF] = RDB[loc1 + IFC_TET_MSH_DF]/dmax;

	      /* If density file was not used set density factor to 1.0 */

	      if( !strcmp(GetText(loc0 + IFC_PTR_OF_RFILE), "-1"))
		WDB[loc1 + IFC_TET_MSH_DF] = 1.0;

	      /* Check that the density factor is below unity*/

	      if(WDB[loc1 + IFC_TET_MSH_DF] > 1.0)
		Die(FUNCTION_NAME, "Larger than unity density factor for interface %s: %E",
		     GetText(loc0 + IFC_PTR_INPUT_FNAME), WDB[loc1 + IFC_TET_MSH_DF]);

	      if(WDB[loc1 + IFC_TET_MSH_DF] < 0)
		Die(FUNCTION_NAME, "Negative density factor for interface %s: %E",
		     GetText(loc0 + IFC_PTR_INPUT_FNAME), WDB[loc1 + IFC_TET_MSH_DF]);


	    }

	}

      /* Next */

      loc1 = NextItem(loc1);

    }

	  

  /***************************************************************************/

}

/*****************************************************************************/
