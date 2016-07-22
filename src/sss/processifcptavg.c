/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcptavg.c                              */
/*                                                                           */
/* Created:       2015/01/23 (VVa)                                           */
/* Last modified: 2016/04/04 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Processes point average multi-physics interfaces             */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCPtAvg:"

/*****************************************************************************/

void ProcessIFCPtAvg(long loc0, long update)
{
  long loc1, mat0, mat, ptr, found, msh, dim, nx, ny, nz, override;
  double x, y, z, rad, xmin, xmax, ymin, ymax, zmin, zmax, dmax;
  double lims[6], matdens;

  /***********************************************************************/

  /***** Link material ***************************************************/

  mat = -1;

  /* Check material pointer (not set in all types) */

  if ((long)RDB[loc0 + IFC_PTR_MAT] < VALID_PTR)
    Die(FUNCTION_NAME, "material pointer not set");
  else if(!update)
    {
      /* Link material and set TMS limits */

      /* Reset found flag */

      found = NO;
	  
      /* Loop over materials and find match */
	  
      mat0 = (long)RDB[DATA_PTR_M0];
      while (mat0 > VALID_PTR)
	{
	  /* Reset pointer */
	      
	  mat = -1;
	      
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
	      /* Link interface and put flag */
	      
	      WDB[mat + MATERIAL_PTR_IFC] = (double)loc0;
	      WDB[mat + MATERIAL_USE_IFC] = (double)YES;		  

	      WDB[loc0 + IFC_PTR_MAT] = (double)mat;
 
	      /* Check if temperature is given */
		  
	      if (RDB[loc0 + IFC_MAX_TEMP] > 0.0)
		{
		  /* Put maximum temperature */
		      
		  if (RDB[loc0 + IFC_MAX_TEMP] > 
		      RDB[mat + MATERIAL_TMS_TMAX])
		    {
		      if(!update)
			{
			  WDB[mat + MATERIAL_TMS_TMAX] = 
			    RDB[loc0 + IFC_MAX_TEMP];

			}
		      else
			Die(FUNCTION_NAME,
			    "Material temperature above TMS majorant for material %s",
			    GetText(mat + MATERIAL_PTR_NAME));
		    }
		  /* Put minimum temperature */
		      
		  if (RDB[loc0 + IFC_MIN_TEMP] < 
		      RDB[mat + MATERIAL_TMS_TMIN])
		    {
		      if(!update)
			{
			  WDB[mat + MATERIAL_TMS_TMIN] = 
			    RDB[loc0 + IFC_MIN_TEMP];

			}
		      else
			Die(FUNCTION_NAME,
			    "Material temperature below TMS minorant for material %s",
			    GetText(mat + MATERIAL_PTR_NAME));

		    }

		  if (!update)
		    {	      

		      /* Set on-the-fly Doppler-broadening mode */
		      /* TMS is set on if the TMS-limits differ */
		      /* TMS is always set on for mixtures */

		      if ((RDB[mat + MATERIAL_TMS_TMIN] <
			   RDB[mat + MATERIAL_TMS_TMAX]) ||
			  ((long)RDB[mat + MATERIAL_PTR_MIX] > VALID_PTR))		      
			{
			  /* Check that material doppler temperature is not already */
			  /* set by tmp-card */

			  if (RDB[mat + MATERIAL_DOPPLER_TEMP] >= 0)
			    Error(mat, "Material temperature set by tmp-card but a temperature distribution is also given by interface %s", GetText(loc0 + IFC_PTR_INPUT_FNAME));

			  /* Set TMS mode on */

			  WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

			  /* Set Doppler-preprocessor off */

			  WDB[mat + MATERIAL_DOPPLER_TEMP] = -1.0;

			}
		      else
			{

			  /* Check that material doppler temperature is not already */
			  /* set by tmp-card */

			  if (RDB[mat + MATERIAL_DOPPLER_TEMP] >= 0)
			    Error(mat, "Material temperature set by tmp-card but a temperature distribution is also given by interface %s", GetText(loc0 + IFC_PTR_INPUT_FNAME));

			  WDB[mat + MATERIAL_DOPPLER_TEMP] = -RDB[loc0 + IFC_MIN_TEMP];
			}
		    }
		}


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

    }
  else
    {
      /* Check TMS limits */

      mat = (long)RDB[loc0 + IFC_PTR_MAT];

      /* Check if temperature is given */
		  
      if (RDB[loc0 + IFC_MAX_TEMP] > 0.0)
	{
	  /* Check maximum temperature */
		      
	  if (RDB[loc0 + IFC_MAX_TEMP] > 
	      RDB[mat + MATERIAL_TMS_TMAX])
		Die(FUNCTION_NAME,
		    "Material temperature above TMS majorant for material %s",
		    GetText(mat + MATERIAL_PTR_NAME));

	  /* Checkminimum temperature */
		      
	  if (RDB[loc0 + IFC_MIN_TEMP] < 
	      RDB[mat + MATERIAL_TMS_TMIN])
	    Die(FUNCTION_NAME,
		"Material temperature below TMS minorant for material %s",
		GetText(mat + MATERIAL_PTR_NAME));	   

	}

    }

  /* Set material density to maximum density */

  mat = (long)RDB[loc0 + IFC_PTR_MAT];

  /* Get interface maximum density */

  dmax = RDB[loc0 + IFC_MAX_DENSITY];

  /* Get material density */

  override = 0;

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
      fprintf(out, "Material %s, matdens %E, dmax %E\n",GetText(mat + MATERIAL_PTR_NAME), matdens, dmax);

    }
  else
    {
      /* If this is the first time this interface is processed we will use */
      /* MATERIAL_ADENS for the density regardless of IFC units */
      /* They will be converted to atomic density later         */

      matdens = RDB[mat + MATERIAL_ADENS];

      /* If base material density is given in different units than the */
      /* IFC-density, we cannot use base density as a majorant */
      /* We will just override the base density with the interface density */

      if (matdens*dmax < 0)
	override = 1;
    }

  /* Only increase the material atomic density before calculating majorants */
  /* (not on updates) */

  /* Only change it upwards, this way the normal density card can be used as */
  /* a maximum density, when setting up a coupled calculation */

  if ((fabs(matdens) < fabs(dmax)) || (override == 1))
    {
      /* IFC density greater than material base density */

      if (!update)
	{
	  WDB[mat + MATERIAL_ADENS] = dmax;
	  matdens = RDB[mat + MATERIAL_ADENS];
	}
      else
	Die(FUNCTION_NAME, "Material density exceeds material majorant for %s",
	    GetText(mat + MATERIAL_PTR_NAME));
    }
  else
    {
      /* Material base density greater than IFC density */

      dmax = matdens;

    }

  fprintf(out, "\nPoint average based interface \"%s\":\n\n",
	  GetText(loc0 + IFC_PTR_INPUT_FNAME));  

  /* Print out the maximum density of the interface */

  fprintf(out, " - Maximum density of interface: %E\n",
	  RDB[loc0 + IFC_MAX_DENSITY]);

  /* Print out the density used for calculating the majorant */

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

  /***********************************************************************/

  /*******************************************************************/

  /***** Average of points *******************************************/

  if(!update)
    {
      /* Get limits */

      xmin = RDB[loc0 + IFC_MESH_XMIN];
      xmax = RDB[loc0 + IFC_MESH_XMAX];
      ymin = RDB[loc0 + IFC_MESH_YMIN];
      ymax = RDB[loc0 + IFC_MESH_YMAX];
      zmin = RDB[loc0 + IFC_MESH_ZMIN];
      zmax = RDB[loc0 + IFC_MESH_ZMAX];

      /* Get exclusion radius and dimension */

      rad = RDB[loc0 + IFC_EXCL_RAD];
      dim = (long)RDB[loc0 + IFC_DIM];

      /* Calculate search mesh size */

      if ((nx = (long)(xmax - xmin)/(2.0*rad)) < 1)
	nx = 1;

      if ((ny = (long)(ymax - ymin)/(2.0*rad)) < 1)
	ny = 1;

      if ((nz = (long)(zmax - zmin)/(2.0*rad)) < 1)
	nz = 1;
	  
      /* Adjust boundaries */
      
      xmin = xmin - 1E-6;
      xmax = xmax + 1E-6;
      ymin = ymin - 1E-6;
      ymax = ymax + 1E-6;
      zmin = zmin - 1E-6;
      zmax = zmax + 1E-6;
	  
      /* 1D and 2D distributions */
      
      if (dim == 1)
	{
	  xmin = -INFTY;
	  xmax = INFTY;
	  nx = 1;
	      
	  ymin = -INFTY;
	  ymax = INFTY;
	  ny = 1;
	}
      else if (dim == 2)
	{
	  zmin = -INFTY;
	  zmax = INFTY;
	  nz = 1;
	}

      /* Put mesh variables */

      lims[0] = xmin;
      lims[1] = xmax;
      lims[2] = ymin;
      lims[3] = ymax;
      lims[4] = zmin;
      lims[5] = zmax;
	  
      /* Create mesh structure */
      
      msh = CreateMesh(MESH_TYPE_CARTESIAN, MESH_CONTENT_PTR, -1,
		       50, 50, 50, lims, -1);
	  
      /* Put pointer */
      
      WDB[loc0 + IFC_PTR_SEARCH_MESH] = (double)msh;

      /* Loop over points */
	  
      loc1 = (long)RDB[loc0 + IFC_PTR_POINTS];
      while (loc1 > VALID_PTR)
	{
	  /* Get coordinates */
	      
	  x = RDB[loc1 + IFC_PT_X];
	  y = RDB[loc1 + IFC_PT_Y];
	  z = RDB[loc1 + IFC_PT_Z];
	      
	  /* Add point in search mesh */
	      
	  AddSearchMesh(msh, loc1, x - rad, x + rad, y - rad, y + rad, 
			z - rad, z + rad);	 
     
	  /* Next point */

	  loc1 = NextItem(loc1);
	}	  


    }

  /* Loop over points to set density factors */

  loc1 = (long)RDB[loc0 + IFC_PTR_POINTS];

  while (loc1 > VALID_PTR)
    {

      /* Convert densities to density factors */

      if (RDB[loc1 + IFC_PT_DF]*dmax < 0.0)
	Error(loc0, "Inconsistent densities given in distribution");
      else
	WDB[loc1 + IFC_PT_DF] = RDB[loc1 + IFC_PT_DF]/dmax;

      loc1 = NextItem(loc1);
    }	  

	  
  /********************************************************************/
}

/*****************************************************************************/
