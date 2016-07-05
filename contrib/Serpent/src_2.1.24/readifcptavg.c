/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcptavg.c                                 */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Reads point average multi-physics interfaces                 */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCPtAvg:"

/*****************************************************************************/

void ReadIFCPtAvg(long loc0, long update)
{
  long loc1;
  long n, nz, nr, np, dim, type;
  double xmin, xmax, ymin, ymax, zmin, zmax, Tmin, Tmax, dmax, x, y, z, d, T;
  double rad, ex;
  char tmpstr[MAX_STR];
  FILE *fp;

  if(update)
    Warn(FUNCTION_NAME, "Interface updating for type %ld not thoroughly tested", (long)RDB[loc0 + IFC_TYPE]);

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

  if(!update)
    {
      WDB[loc0 + IFC_PTR_MAT] = (double)PutText(tmpstr);
    }

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

    

  /***********************************************************************/

  /***** Average of points *******************************************/
	  
  /* Reset limiting values */

  Tmax = -INFTY;
  Tmin = INFTY;
  dmax = 0.0;

  /* Reset mesh boundaries */
      
  xmin =  INFTY;
  xmax = -INFTY;
  ymin =  INFTY;
  ymax = -INFTY;
  zmin =  INFTY;
  zmax = -INFTY;

  /* Read dimension, exclusion radius and exponent (TODO: Tohon */
  /* ne muut tyypit: eksponentiaali, ym.) */

  if (fscanf(fp, "%ld %lf %lf", &dim, &rad, &ex) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Put dimension */

  if ((dim < 1) || (dim > 3))
    Error(loc0, "Error in dimension");
  else
    WDB[loc0 + IFC_DIM] = (double)dim;

  /* Exclusion radius */
	  
  if (rad > 0.0)
    WDB[loc0 + IFC_EXCL_RAD] = rad;
  else
    Error(loc0, "Invalid exclusion radius entered");

  /* Exponent */
	  
  WDB[loc0 + IFC_EXP] = ex;

  /* Read number of points */
	  
  if (fscanf(fp, "%ld", &np) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Check number of points */
	  
  if (np > 0)
    WDB[loc0 + IFC_NP] = (double)np;
  else
    Error(loc0, "Invalid number of points entered");
    

  /* Get pointer to points */

  if(update)
    loc1 = (long)RDB[loc0 + IFC_PTR_POINTS];    
  else
    loc1 = -1;

  /* Loop over points and read data */
      
  for (n = 0; n < np; n++)
    {
      /* Avoid compiler warning */
	      
      x = 0.0;
      y = 0.0;
      z = 0.0;
	      
      /* Read values */
	      
      if (dim == 3)
	{
	  if (fscanf(fp, "%lf %lf %lf %lf %lf", 
		     &x, &y, &z, &d, &T) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");
	}
      else if (dim == 2)
	{
	  if (fscanf(fp, "%lf %lf %lf %lf", &x, &y, &d, &T) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");		   
	}
      else if (dim == 1)
	{
	  if (fscanf(fp, "%lf %lf %lf", &z, &d, &T) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");
	}
      else
	Die(FUNCTION_NAME, "error in dimension");
	  
      /* Allocate memory for point */
      if(!update)
	loc1 = NewItem(loc0 + IFC_PTR_POINTS, IFC_PT_LIST_BLOCK_SIZE);
	      
      /* Put data */

      WDB[loc1 + IFC_PT_DF] = d;
      WDB[loc1 + IFC_PT_TMP] = T;
      WDB[loc1 + IFC_PT_X] = x;
      WDB[loc1 + IFC_PT_Y] = y;
      WDB[loc1 + IFC_PT_Z] = z;

      /* Compare to limits */

      if (x - rad < xmin)
	xmin = x - rad;
      if (x + rad> xmax)
	xmax = x + rad;
	      
      if (y - rad < ymin)
	ymin = y - rad;
      if (y + rad > ymax)
	ymax = y + rad;
	      
      if (z - rad < zmin)
	zmin = z - rad;
      if (z + rad > zmax)
	zmax = z + rad;

      if (fabs(d) > fabs(dmax))
	dmax = d;
	      
      if (T > Tmax)
	Tmax = T;

      if (T < Tmin)
	Tmin = T;

      /* Next point */
      if(update)
	loc1 = NextItem(loc1);
    }

  /* Put maximum density and temperature                 */
  /* For updates, these are checked in processifcptavg.c */

  WDB[loc0 + IFC_MAX_DENSITY] = dmax;
  WDB[loc0 + IFC_MAX_TEMP] = Tmax;
  WDB[loc0 + IFC_MIN_TEMP] = Tmin;

  if(!update)
    {
      /* Put boundaries */

      WDB[loc0 + IFC_MESH_XMIN] = xmin;
      WDB[loc0 + IFC_MESH_XMAX] = xmax;
      WDB[loc0 + IFC_MESH_YMIN] = ymin;
      WDB[loc0 + IFC_MESH_YMAX] = ymax;
      WDB[loc0 + IFC_MESH_ZMIN] = zmin;
      WDB[loc0 + IFC_MESH_ZMAX] = zmax;

    }
  /* Set TMS on */

  if (Tmax > 0.0)
    if(!update)
      WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;

  /*******************************************************************/

  /* Close file */

  fclose(fp);

}

/*****************************************************************************/
