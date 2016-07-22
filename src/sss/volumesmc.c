/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : volumesmc.c                                    */
/*                                                                           */
/* Created:       2010/11/10 (JLe)                                           */
/* Last modified: 2015/05/08 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Calculates volumes by Monte Carlo simulation                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "VolumesMC"

/*****************************************************************************/

void VolumesMC()
{
  long cell, mat, mat0, ptr, loc0, nt, idx, m, nmax, nb, id, dim;
  unsigned long seed;
  double xmin, xmax, ymin, ymax, zmin, zmax, x, y, z, u, v, w, vol, f, T;
  double max, err, tmax, emax, t, val, est, diff, df;
  char tmpstr[MAX_STR], fname[MAX_STR];
  FILE *fp;

  /* Cut-offs */

  nmax = (long)RDB[DATA_VOLUME_MC_NMAX];
  tmax = RDB[DATA_VOLUME_MC_TMAX];
  emax = RDB[DATA_VOLUME_MC_EMAX];

  /* Check if calculation is required */

  if ((nmax < 0) && (tmax < 0.0) && (emax < 0.0))
    return;

  fprintf(out, "Calculating material volumes by Monte Carlo...\n");

  /***************************************************************************/

  /***** Allocate memory for statistics **************************************/

  /* Reset stat pointer */

  loc0 = -1;

  /* Loop over materials */

  mat = RDB[DATA_PTR_M0];		
  while (mat > VALID_PTR)
    {
      /* Allocate memory for volume and density */

      ptr = NewStat("VOLUME", 1, 1);
      WDB[mat + MATERIAL_PTR_MC_VOLUME] = (double)ptr;      

      /* Remember pointer to first */

      if (loc0 < 0)
	loc0 = ptr;

      ptr = NewStat("DENSITY", 1, 1);
      WDB[mat + MATERIAL_PTR_MC_DENSITY] = (double)ptr;      

      /* Next */
      
      mat = NextItem(mat);
    }

  /* Expand PRIVA, BUF and RES2 arrays for OpenMP parallel calculation */

  ExpandPrivateArrays();

  /***************************************************************************/

  /***** Monte carlo volume calculation **************************************/

  /* Get geometry boundaries */

  xmin = RDB[DATA_GEOM_MINX];
  xmax = RDB[DATA_GEOM_MAXX];
  ymin = RDB[DATA_GEOM_MINY];
  ymax = RDB[DATA_GEOM_MAXY];
  zmin = RDB[DATA_GEOM_MINZ];
  zmax = RDB[DATA_GEOM_MAXZ];

  /* Calculate total volume */

  vol = (xmax - xmin)*(ymax - ymin);
  
  /* Check dimension */
  
  if ((dim = (long)RDB[DATA_GEOM_DIM]) == 3)
    {
      /* Add axial dimension to volume */
      
      vol = vol*(zmax - zmin);
    }

  /* Start timer */

  StartTimer(TIMER_VOLUME_CALC);

  /* Number of samples per batch */

  if (nmax > 0)
    nt = (long)(nmax/100.0);
  else
    nt = 100000;

  /* Reset counters */

  nb = 0;

  /* Loop over batches */

  while (1 != 2)
    {
      /* Start parallel timer */

      StartTimer(TIMER_OMP_PARA);

#ifdef OPEN_MP
#pragma omp parallel private (m, idx, seed, ptr, x, y, z, u, v, w, cell, mat, id, f, T)
#endif
      {
	/* Loop over points */

#ifdef OPEN_MP
#pragma omp for      
#endif
	for (m = 0; m < nt; m++)
	  {
	    /* Get OpenMP thread num */

	    id = OMP_THREAD_NUM;

	    /* Calculate index */
	    
	    idx = nb*nt + m + 1;

	    /* Init random number sequence */
      
	    seed = ReInitRNG(idx);
      	    SEED[id*RNG_SZ] = seed;
      	    
	    /* Sample point */

	    x = RandF(id)*(xmax - xmin) + xmin;
	    y = RandF(id)*(ymax - ymin) + ymin;

	    if (dim == 3)
	      z = RandF(id)*(zmax - zmin) + zmin;
	    else
	      z = 0.0;

	    /* Sample direction (this is necessary for STL geometries) */

	    IsotropicDirection(&u, &v, &w, id);
	    
	    /* Find position */
	    
	    if ((cell = WhereAmI(x, y, z, u, v, w, id)) < 0)
	      Error(0, "Geometry error at %E %E %E", x, y, z);
	    
	    /* Check if cell has material */
	    
	    if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
	      {
		/* Get material pointer */

		mat = MatPtr(mat, id);
		
		/* Reset density and temperature */

		f = 1.0;
		T = 0.0;

		/* Get point from interface */
		
		IFCPoint(mat, &f, &T, id);

		/* Check for undefined density */

		if (f < 0.0)
		  continue;

		/* Score point estimators */
		
		ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		AddBuf1D(vol/((double)nt), 1.0, ptr, id, 0);
		
		ptr = (long)RDB[mat + MATERIAL_PTR_MC_DENSITY];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		AddBuf1D(f*vol/((double)nt), 1.0, ptr, id, 0);

		/* Check if material is divided */

		if ((mat = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) 
		    > VALID_PTR)
		  {
		    /* Add to stat */

		    ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
		    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		    AddBuf1D(vol/((double)nt), 1.0, ptr, id, 0);

		    ptr = (long)RDB[mat + MATERIAL_PTR_MC_DENSITY];
		    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		    AddBuf1D(f*vol/((double)nt), 1.0, ptr, id, 0);
		  }		  
	      }
	  }
      }
    
      /* Stop parallel timer */

      StopTimer(TIMER_OMP_PARA);

      /* Update batch number */

      nb++;

      /* Reduce scoring buffer */

      ReduceBuffer();

      /* Maximum error */

      max = 0.0;

      /* Loop over statistics */
      
      ptr = loc0;      
      while (ptr > VALID_PTR)
	{
	  /* Collect data */

	  val = BufVal(ptr, 0);
	  AddStat(val, ptr, 0);

	  /* Get relative error */

	  err = RelErr(ptr, 0);

	  /* Compare maximum error */

	  if (err > max)
	    max = err;	     
	  
	  /* Next */
	  
	  ptr = NextItem(ptr);
	}

      /* Clear buffer */
	  
      ClearBuf();
	
      /* Get time */

      t = TimerVal(TIMER_VOLUME_CALC);

      if (nb == 1)
	fprintf(out, "\nEstimated calculation time: %s\n",
		TimeStr(t*nmax/nt));

      /* Check batch cut-off */
      
      if ((nmax > 0) && (nb*nt >= nmax))
	{
	  fprintf(out, "Realized calculation time:  %s\n\n",
		  TimeStr(t));
	  break;
	}
      
      /* Check time cut-off */

      if ((tmax > 0.0) && (t > tmax))
	break;
      
      /* Check error cut-off */
      
      if ((emax > 0.0) && (max > 0.0) && (max < emax))
	break;
    }

  /* Stop timer */

  StopTimer(TIMER_VOLUME_CALC);

  /***************************************************************************/

  /***** Print volumes *******************************************************/

  /* Check dimensions */

  if ((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES)
    {
      if (dim == 3)
	fprintf(out, "Volumes (in cm3) :\n\n");
      else
	fprintf(out, "Volumes (2D problem, the values are in cm2) :\n\n");
    }

  /* Check mode */

  if ((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES)
    {
      /* Set file name */
      
      sprintf(fname, "%s.mvol", GetText(DATA_PTR_INPUT_FNAME));	      
      
      /* Open file */
      
      if ((fp = fopen(fname, "w")) == NULL)
	Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Print comments and card name */

      fprintf(fp, "%% --- Material volumes:\n\n");
      fprintf(fp, "%% Produced %s by MC volume calculation routine by\n",
	      GetText(DATA_PTR_DATE));
      fprintf(fp, "%% sampling %ld random points in the geometry.\n\n", nb*nt);
      fprintf(fp, "set mvol\n\n");
    }
  else
    fp = NULL;

  /* Put volumes */

  mat = RDB[DATA_PTR_M0]; 		
  while (mat > VALID_PTR)
    {
      /* Get density */

      ptr = (long)RDB[mat + MATERIAL_PTR_MC_DENSITY];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      df = Mean(ptr, 0);

      /* Pointer to volume */

      ptr = (long)RDB[mat + MATERIAL_PTR_MC_VOLUME];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Given or calculated volume */

      vol = RDB[mat + MATERIAL_VOLUME];

      /* Estimated volume, error and relative difference */

      est = Mean(ptr, 0);
      err = RelErr(ptr, 0);
      diff = est/vol - 1.0;

      /* Pointer to parent */

      mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT];
      
      /* Get number of zones */

      m = 0;
      if (mat0 > VALID_PTR) 
	m = (long)RDB[mat0 + MATERIAL_DIV_N_TOT_ZONES];

      /* Check mode */

      if ((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES)
	{
	  /* Print if no or more than 1 zones */

	  if (m != 1)
	    {
	      /* Get printed name */

	      sprintf(tmpstr, "Material %s", GetText(mat + MATERIAL_PTR_NAME));
	  
	      /* Print */
	      
	      if ((vol > ZERO) && (vol < INFTY))
		{
		  if (diff > 1.0)
		    fprintf(out, "%-22s : %1.4E %1.4E (%7.5f) :   > 1.0 ", 
			    tmpstr, vol, est, err);
		  else if (diff < -1.0)
		    fprintf(out, "%-22s : %1.4E %1.4E (%7.5f) :  < -1.0 ", 
			    tmpstr, vol, est, err);
		  else
		    fprintf(out, "%-22s : %1.4E %1.4E (%7.5f) : %8.5f ", 
			    tmpstr, vol, est, err, diff);
		}
	      else if (vol == 0.0)
		fprintf(out, "%-22s : %1.4E %1.4E (%7.5f) :      N/A ", 
			tmpstr, vol, est, err);
	      else
		fprintf(out, 
			"%-22s :        N/A %1.4E (%7.5f) :      N/A ", 
			tmpstr, est, err);
	    
	      if (fabs(diff) > 1.96*err)
		fprintf(out, "* ");
	      else
		fprintf(out, "  ");

	      if (est > 0.0)
		fprintf(out, "(%5.1f %% den.)\n", 100*df/est);
	      else
		fprintf(out, "\n");
	    }

	  /* Print volume to file */
	  
	  if ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT)
	    {
	      if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
		fprintf(fp, "%-10s ", GetText(mat0 + MATERIAL_PTR_NAME));
	      else
		fprintf(fp, "%-10s ", GetText(mat + MATERIAL_PTR_NAME));
	      
	      fprintf(fp, "%6ld %1.5E %% (%5.3f)\n", 
		      (long)RDB[mat + MATERIAL_DIV_ZONE_IDX], Mean(ptr, 0),
		      RelErr(ptr, 0));
	    }
	}
      else
	{
	  /* Put volume if not given */

	  if (RDB[mat + MATERIAL_VOLUME_GIVEN] < 0.0)
	    WDB[mat + MATERIAL_VOLUME] = Mean(ptr, 0);
	}

      /* Next */
      
      mat = NextItem(mat);
    }

  /* Close file */

  if (fp != NULL)
    fclose(fp);

  /***************************************************************************/

  if ((long)RDB[DATA_VOLUME_CALCULATION_MODE] == YES)
    fprintf(out, "\nVolumes written in file \"%s\"\n", fname);

  /* Exit subroutine */

  fprintf(out, "\nOK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
