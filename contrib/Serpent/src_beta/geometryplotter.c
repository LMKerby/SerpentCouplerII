/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : geometryplotter.c                              */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2015/10/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Plots geometry                                               */
/*                                                                           */
/* Comments: - Lisää optio jolla saa valita rajapintojen vs. materiaali-     */
/*             rajojen korostuksen.                                          */
/*                                                                           */
/*           - Tätä kutsutaan monesta eri paikkaa ja eri tapausten käsittely */
/*             on kamala sekamelska. Jaa koko roska jossain vaiheessa perus- */
/*             osaan joka plottaa pelkän geometrian ja erillisiin aliohjel-  */
/*             miin jotka hoitaa muut hommat.                                */
/*                                                                           */
/*           - Siirrettiin iso pätkä paletteja valmistelevaa koodia          */
/*             silmukoiden sisään 9.10.2015. Kamala sekasotku, ihme jos      */
/*             vielä toimii.                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifndef NO_GFX_MODE
#include <gd.h>
#endif

#define FUNCTION_NAME "GeometryPlotter:"

#define COLOURS 256

/* Local subroutines */

#ifndef NO_GFX_MODE

void FlashPalette(gdImagePtr im, long *r, long *g, long *b, long type);

#endif

/*****************************************************************************/

void GeometryPlotter(long ini)
{
#ifndef NO_GFX_MODE

  long n, m, gpl, xp, yp, cell, mat, nerr, old, mode, id, ptr, i, nc0;
  long R[COLOURS], G[COLOURS], B[COLOURS], count, nplot, ext;
  long ncol, nifc, nf, nftot, fiss, wwp;
  long **matrix1, **matrix2;
  double xmin, xmax, ymin, ymax, zmin, zmax, tmp, x, y, z, u, v, w, d, pixw;
  double dummy, f, T, **spt, fmean, fmin, fmax;
  gdImagePtr im;
  FILE *fp;
  long palette[COLOURS];
  unsigned long seed;
  char fname[MAX_STR];

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Pointer to geometry plot */

  if ((gpl = RDB[DATA_PTR_GPL0]) < 1)
    {
      /* Check stop mode */

      if ((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_GEOM)
	exit(-1);
      
      /* Return */

      return;
    }

  /* Check track plotter mode */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    {
      /* Check if done */

      if ((long)RDB[DATA_SIMULATION_COMPLETED] == NO)
	return;
      else if ((long)RDB[DATA_TRACK_PLOT_ANIM] == YES)
	fprintf(out, "Plotting particle tracks in geometry:\n\n");
    }
  else if (ini == YES)
    fprintf(out, "Plotting geometry:\n\n");
  else if ((long)RDB[DATA_SOURCE_PT_ANIM] == NO)
    return;

  /**************************************************************************/

  /***** Initialization *****************************************************/

  /* Plot weight window boundaries */

  if ((long)RDB[DATA_USE_WEIGHT_WINDOWS] == YES)
    wwp = YES;
  else
    wwp = NO;

  wwp = NO;
  
  /* Set plotter mode */
  
  WDB[DATA_PLOTTER_MODE] = (double)YES;

  /* Reset exit flag */

  ext = NO;

  /* Reset dummy weight used for albedo boundary conditions */

  dummy = -1.0;

  /***************************************************************************/

  /***** Loop over track plotter frames **************************************/

  /* Nuber of frames */

  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
    nftot = (long)RDB[DATA_TRACK_PLOT_FRAMES];
  else
    nftot = 1;

  /* Reset fission flag for track plotter */

  fiss = NO;

  /* Loop over frames */

  for (nf = 0; nf < nftot; nf++)
    {
      /***********************************************************************/
      
      /***** Plot geometry ***************************************************/
      
      /* Print progress */

      if (nftot > 200)
	{
	  if (!(nf % 20))
	    fprintf(out, " %3.0f%% complete\n", 
		    100.0*((double)nf)/((double)nftot));     
	}
      else if (nftot > 100)
	{
	  if (!(nf % 10))
	    fprintf(out, " %3.0f%% complete\n", 
		    100.0*((double)nf)/((double)nftot));     
	}
      else if (nftot > 50)
	{
	  if (!(nf % 5))
	    fprintf(out, " %3.0f%% complete\n", 
		    100.0*((double)nf)/((double)nftot));     
	}
      else if (nftot > 1)
	fprintf(out, " %3.0f%% complete\n", 
		100.0*((double)nf)/((double)nftot));     

      /* Get pointer to geometry plots */

      gpl = RDB[DATA_PTR_GPL0];
      CheckPointer(FUNCTION_NAME, "(gpl)", DATA_ARRAY, gpl);

      /* Reset counter and get number of plots */

      count = 0;
      nplot = ListSize(gpl);

      /* Loop over plots */
      
      gpl = FirstItem(gpl);
      while (gpl > VALID_PTR)
	{
	  /* Skip importance mesh plots if no weight windows */

	  if ((long)RDB[gpl + GPL_IMP_SCALE] > 0)
	    if ((long)RDB[DATA_PTR_WWD0] < VALID_PTR)
	      {
		/* Pointer to next */

		gpl = NextItem(gpl);

		/* Cycle loop */

		continue;
	      }

	  /* Loop over OpenMP threads */

	  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	    {
	      /* Init random number sequence */
	      
	      seed = ReInitRNG(id);
	      SEED[id*RNG_SZ] = seed;
	    }

	  /* Reduce number of colors if IFC or source point animation */
	  /* is in use */

	  if ((ptr = (long)RDB[DATA_PTR_IFC0]) > VALID_PTR)
	    {
	      ncol = (long)(((double)COLOURS)/2.0);
	      nifc = (long)(((double)ncol - 10)/((double)ListSize(ptr)));
	    }
	  else if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
	    {
	      ncol = (long)(((double)COLOURS)/2.0);
	      nifc = 0;
	    }
	  else if (((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0) && 
		   ((long)RDB[DATA_TRACK_PLOT_ANIM] == YES))
	    {
	      ncol = COLOURS - TRACK_PLOT_NCOL;
	      nifc = 0;
	    }
	  else
	    {
	      ncol = COLOURS;
	      nifc = 0;
	    }
	  
	  /* Random colors */
	  
	  for(n = 5; n < COLOURS; n++)
	    {
	      R[n] = (long)(255.0*(RandF(0) + 1.0)/2);
	      G[n] = (long)(255.0*(RandF(0) + 1.0)/2);
	      B[n] = (long)(255.0*(RandF(0) + 1.0)/2);
	    }
	  
	  /* Void colour */
	  
	  R[0] = 0;
	  G[0] = 0;
	  B[0] = 0;
	  
	  /* No cell error */
	  
	  R[1] = 0;
	  G[1] = 255;
	  B[1] = 0;
	  
	  /* Multiple cell error */
	  
	  R[2] = 255;
	  G[2] = 0;
	  B[2] = 0;
	  
	  /* Pointer error */
	  
	  R[3] = 255;
	  G[3] = 255;
	  B[3] = 0;
	  
	  /* Underined density factor */
	  
	  R[4] = 255;
	  G[4] = 0;
	  B[4] = 255;
	  
	  /*******************************************************************/

	  /***** Set fixed colors ********************************************/
	  
	  /* User-defined material colours */
	  
	  n = 5;
	  i = NO;
	  
	  mat = RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Check pointer to interface */
	      
	      if ((long)RDB[mat + MATERIAL_PTR_IFC] < VALID_PTR)
		{
		  /* Check if colour is defined */
		  
		  if ((m = (long)RDB[mat + MATERIAL_RGB]) > 0)
		    {
		      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) 
			  > VALID_PTR)
			WDB[mat + MATERIAL_COLOUR_IDX] = 
			  RDB[ptr + MATERIAL_COLOUR_IDX];
		      else if (n < ncol)
			{
			  /* Check index */
			  
			  if ((n < 0) || (n > COLOURS - 1))
			    Die(FUNCTION_NAME, "Index error");
			  
			  R[n] = (long)(m/1000000.0);
			  G[n] = (long)((m - 1000000.0*R[n])/1000.0);
			  B[n] = (long)(m - 1000000.0*R[n] - 1000.0*G[n]);
			  
			  /* Set index */
			  
			  WDB[mat + MATERIAL_COLOUR_IDX] = (double)n;
			  
			  /* Next colour */
			  
			  n++;
			}
		      else
			{
			  /* Use random index */
			  
			  WDB[mat + MATERIAL_COLOUR_IDX] = 
			    (double)((long)(RandF(0)*(ncol - 5))) + 5.0;
			  
			  /* Flag for warning */
			  
			  i = YES;
			}
		    }
		}
	      
	      /* Next material */
	      
	      mat = NextItem(mat);
	    }
	  
	  /* Check if all are used */
	  
	  if (i == YES)
	    Note(0, 
		 "Number of pre-assigned material colors exceeds max (252)");
    
	  /* Put remaining colour indexes */
	  
	  m = n;
	  
	  mat = RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Check pointer to interface */
	      
	      if ((long)RDB[mat + MATERIAL_PTR_IFC] < VALID_PTR)
		{
		  /* Cycle if all are used */
		  
		  if (m > ncol - 1)
		    m = n;
		  
		  /* One more check */
		  
		  if (m > ncol - 1)
		    m = ncol - 1;
		  
		  /* Check if index is not already given */
		  
		  if ((long)RDB[mat + MATERIAL_COLOUR_IDX] == 0)
		    WDB[mat + MATERIAL_COLOUR_IDX] = (double)m++;
		}
	      
	      /* Next material */
	      
	      mat = NextItem(mat);
	    }

	  /* Set total number of colors */

	  nc0 = n;

	  /*******************************************************************/

	  /***** Set IFC colors **********************************************/

	  /* User-defined material colours */
	  
	  n = ncol;
	  
	  mat = RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Check pointer to interface */

	      if ((long)RDB[mat + MATERIAL_PTR_IFC] > VALID_PTR)
		{
		  /* Number of colors */

		  nc0 = ncol;

		  /* Check if colour is defined */
		  
		  if ((m = (long)RDB[mat + MATERIAL_RGB]) > 0)
		    {
		      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > 
			  VALID_PTR)
			WDB[mat + MATERIAL_COLOUR_IDX] = 
			  RDB[ptr + MATERIAL_COLOUR_IDX];
		      else if (n < COLOURS)
			{
			  /* Check index */
			  
			  if ((n < 0) || (n > COLOURS - 1))
			    Die(FUNCTION_NAME, "Index error");
			  
			  R[n] = (long)(m/1000000.0);
			  G[n] = (long)((m - 1000000.0*R[n])/1000.0);
			  B[n] = (long)(m - 1000000.0*R[n] - 1000.0*G[n]);
			  
			  /* Set index */
			  
			  WDB[mat + MATERIAL_COLOUR_IDX] = (double)n;
			  
			  /* Next colour */
			  
			  n = n + nifc;
			}
		      else
			WDB[mat + MATERIAL_COLOUR_IDX] = 0;
		    }
		}
	      
	      /* Next material */
	      
	      mat = NextItem(mat);
	    }
	  
	  /* Put remaining colour indexes */
	  
	  m = n;
	  
	  mat = RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Check pointer to interface */
	      
	      if ((long)RDB[mat + MATERIAL_PTR_IFC] > VALID_PTR)
		{
		  /* Cycle if all are used */
		  
		  if (m > COLOURS - 1)
		    m = n;
		  
		  /* One more check */
		  
		  if (m > COLOURS - 1)
		    m = COLOURS - 1;
		  
		  /* Check if index is not already given */
		  
		  if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > 
		      VALID_PTR)
		    WDB[mat + MATERIAL_COLOUR_IDX] = 
		      RDB[ptr + MATERIAL_COLOUR_IDX];
		  else if ((long)RDB[mat + MATERIAL_COLOUR_IDX] == 0)
		    {
		      WDB[mat + MATERIAL_COLOUR_IDX] = (double)m;
		      m = m + nifc;
		    }
		}
	      
	      /* Next material */
	      
	      mat = NextItem(mat);
	    }
	  
	  /* Create shades */
	  
	  mat = RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Check pointer to interface */
	      
	      if ((long)RDB[mat + MATERIAL_PTR_IFC] > VALID_PTR)
		{
		  /* Get index to color */
		  
		  if ((n = (long)RDB[mat + MATERIAL_COLOUR_IDX]) < ncol)
		    Die(FUNCTION_NAME, "Indexing error");
		  
		  /* Loop over shades */
		  
		  for (m = 1; m < nifc; m++)
		    {
		      /* Check index */
		      
		      if ((n < 0) || (n > COLOURS - 1) || 
			  (n + nifc - m < 0) || (n + nifc - m > COLOURS - 1))
			Die(FUNCTION_NAME, "Index error");
		      
		      /* Put colours */
		      
		      R[n + nifc - m] = 
			R[n]*((double)(m - 1))/((double)(nifc - 1));
		      G[n + nifc - m] = 
			G[n]*((double)(m - 1))/((double)(nifc - 1));
		      B[n + nifc - m] = 
			B[n]*((double)(m - 1))/((double)(nifc - 1));
		    }
		}
	      
	      /* Next material */
	      
	      mat = NextItem(mat);
	    }
	  
	  /*******************************************************************/

	  /***** Adjust colors for track plotting ****************************/

	  /* Check if tracks are recorded for plotting */

	  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
	    {
	      /* Brighten up... */
	      
	      for(n = 4; n < COLOURS; n++)
		{
		  /* Check index */
		  
		  if ((n < 0) || (n > COLOURS - 1))
		    Die(FUNCTION_NAME, "Index error");

		  if ((R[n] = (long)(R[n] + 50)) > 255)
		    R[n] = 255;
		  if ((G[n] = (long)(G[n] + 50)) > 255)
		    G[n] = 255;
		  if ((B[n] = (long)(B[n] + 50)) > 255)
		    B[n] = 255;
		}
	      
	      /* Particle */
	      
	      R[0] = 0;
	      G[0] = 0;
	      B[0] = 0;
	      
	      /* Tail */
      
	      if ((long)RDB[DATA_TRACK_PLOT_ANIM] == YES)
		{
		  for (n = 0; n < TRACK_PLOT_NCOL; n++)
		    {
		      f = ((double)n/((double)TRACK_PLOT_NCOL));
		      
		      /* Hot */
		      
		      if (1 == 2)
			{
			  R[ncol + n] = (long)(2.0*f*255.0);
			  G[ncol + n] = (long)(2.0*f*255.0) - 200;
			  B[ncol + n] = (long)(2.0*f*255.0) - 300;
			}
		      
		      /* Purple */
		      
		      else if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
			{
			  R[ncol + n] = (long)(2.0*f*255.0);
			  G[ncol + n] = (long)(2.0*f*255.0) - 200;
			  B[ncol + n] = (long)(2.0*f*255.0);
			}
		      
		      /* Green */
		      
		      else
			{
			  R[ncol + n] = (long)(2.0*f*255.0) - 200;
			  G[ncol + n] = (long)(2.0*f*255.0);
			  B[ncol + n] = (long)(2.0*f*255.0) - 200;
			}
		      
		      /* Check limits */
		      
		      if (R[ncol + n] > 255)
			R[ncol + n] = 255;
		      else if (R[ncol + n] < 0)
			R[ncol + n] = 0;
		      
		      if (G[ncol + n] > 255)
			G[ncol + n] = 255;
		      else if (G[ncol + n] < 0)
			G[ncol + n] = 0;
		      
		      if (B[ncol + n] > 255)
			B[ncol + n] = 255;
		      else if (B[ncol + n] < 0)
			B[ncol + n] = 0;
		    }
		}
	    }
	  
	  /*******************************************************************/

	  /***** Set colors for importance mesh ******************************/

	  /* Check mode */

	  if ((long)RDB[gpl + GPL_IMP_SCALE] > 0)
	    {
	      /* Override palette */
	      
	      MakePalette(&R[nc0], &G[nc0], &B[nc0], ncol - nc0, 
			  PALETTE_BLUE_RED);
	    }
	  
	  /*******************************************************************/

	  /***** Put colors for source point animation ***********************/

	  /* Check index */

	  if ((ncol < 1) || (n > COLOURS))
	    Die(FUNCTION_NAME, "Index error");
	  
	  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
	    MakePalette(&R[ncol], &G[ncol], &B[ncol], ncol, 
			(long)RDB[DATA_SOURCE_PT_ANIM_PALETTE]);
	  
	  /*******************************************************************/

	  /* Boundary mode (1 = plot cell boundaries, 2 = plot material */
	  /*                boundaries, 3 = plot both) */

	  if ((((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0) || 
	       ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)) &&
	      ((long)RDB[gpl + GPL_IMP_SCALE] == 0))
	    mode = 0;
	  else
	    mode = (long)RDB[gpl + GPL_PLOT_BOUND];

	  /* Print progress */

	  if ((nftot == 1) && ((long)RDB[DATA_PART_PTR_SOURCE] < VALID_PTR))
	    fprintf(out, " %3.0f%% complete\n", 
		    100.0*(count++)/((double)nplot));

	  /* Set file name */
	  
	  if ((long)RDB[DATA_BURN_STEP] > 0)
	    {
	      if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
		{
		  /* Check number of frames */
		  
		  if (nftot < 2)
		    sprintf(fname, "%s_trck%ld_bu%ld.png", 
			    GetText(DATA_PTR_INPUT_FNAME),
			    (long)RDB[gpl + GPL_IDX], 
			    (long)RDB[DATA_BURN_STEP]);
		  else 
		    sprintf(fname, "%s_trck%ld_bu%ld_frame%s.png", 
			    GetText(DATA_PTR_INPUT_FNAME),
			    (long)RDB[gpl + GPL_IDX], (long)RDB[DATA_BURN_STEP],
			    IdxStr(nftot - nf, nftot));
		}
	      else
		sprintf(fname, "%s_geom%ld_bu%ld.png", 
			GetText(DATA_PTR_INPUT_FNAME),
			(long)RDB[gpl + GPL_IDX], (long)RDB[DATA_BURN_STEP]);
	    }
	  else
	    {
	      if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
		{
		  /* Check number of frames */

		  if (nftot < 2)
		    sprintf(fname, "%s_trck%ld.png", 
			    GetText(DATA_PTR_INPUT_FNAME),
			    (long)RDB[gpl + GPL_IDX]);
		  else
		    sprintf(fname, "%s_trck%ld_frame%s.png", 
			    GetText(DATA_PTR_INPUT_FNAME),
			    (long)RDB[gpl + GPL_IDX], 
			    IdxStr(nftot - nf, nftot));
		}
	      else if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
		{
		  sprintf(fname, "%s_spt%ld_frame%s.png", 
			  GetText(DATA_PTR_INPUT_FNAME),
			  (long)RDB[gpl + GPL_IDX], 
			  IdxStr((long)RDB[DATA_CYCLE_IDX],
			(long)RDB[DATA_CRIT_CYCLES] +
				 (long)RDB[DATA_CRIT_SKIP]));
		}
	      else
		sprintf(fname, "%s_geom%ld.png", GetText(DATA_PTR_INPUT_FNAME),
			(long)RDB[gpl + GPL_IDX]);
	    }

	  /* Reset error counter */

	  nerr = 0;

	  /********************************************************************/

	  /***** Set boundaries and image size ********************************/
      
	  /* Get boundaries of geometry */
	  
	  xmin = RDB[DATA_GEOM_MINX];
	  xmax = RDB[DATA_GEOM_MAXX];
	  ymin = RDB[DATA_GEOM_MINY];
	  ymax = RDB[DATA_GEOM_MAXY];
	  zmin = RDB[DATA_GEOM_MINZ];
	  zmax = RDB[DATA_GEOM_MAXZ];
	  
	  /* Check if boundaries are set by user */
	  
	  if (RDB[gpl + GPL_XMIN] != -INFTY)
	    xmin = RDB[gpl + GPL_XMIN];
	  
	  if (RDB[gpl + GPL_XMAX] !=  INFTY)
	    xmax = RDB[gpl + GPL_XMAX];
	  
	  if (RDB[gpl + GPL_YMIN] != -INFTY)
	    ymin = RDB[gpl + GPL_YMIN];
	  
	  if (RDB[gpl + GPL_YMAX] !=  INFTY)
	    ymax = RDB[gpl + GPL_YMAX];
	  
	  if (RDB[gpl + GPL_ZMIN] != -INFTY)
	    zmin = RDB[gpl + GPL_ZMIN];
	  
	  if (RDB[gpl + GPL_ZMAX] !=  INFTY)
	    zmax = RDB[gpl + GPL_ZMAX];
	  
	  /* Check boundaries and swap */
	  
	  if (xmin > xmax)
	    {
	      tmp = xmax;
	      xmax = xmin;
	      xmin = tmp;
	    }
	  
	  if (ymin > ymax)
	    {
	      tmp = ymax;
	      ymax = ymin;
	      ymin = tmp;
	    }
	  
	  if (zmin > zmax)
	    {
	      tmp = zmax;
	      zmax = zmin;
	      zmin = tmp;
	    }
	  
	  /* Move limits away from boundaries to avoid numerical errors in */
	  /* some systems. */

	  xmin = xmin + EXTRAP_L;
	  xmax = xmax - EXTRAP_L;
	  ymin = ymin + EXTRAP_L;
	  ymax = ymax - EXTRAP_L;
	  zmin = zmin + EXTRAP_L;
	  zmax = zmax - EXTRAP_L;
	  
	  /* Set image size */
	  
	  xp = (long)RDB[gpl + GPL_PIX_X];
	  yp = (long)RDB[gpl + GPL_PIX_Y];
	  
	  /* Allocate memory */
	  
	  matrix1 = (long **)Mem(MEM_ALLOC, xp, sizeof(long *));	
	  matrix2 = (long **)Mem(MEM_ALLOC, xp, sizeof(long *));	
	  
	  for(n = 0; n < xp; n++)
	    {
	      matrix1[n] = (long *)Mem(MEM_ALLOC, yp, sizeof(long));
	      matrix2[n] = (long *)Mem(MEM_ALLOC, yp, sizeof(long));
	    }

	  /* Allocate memory for source points */

	  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
	    {
	      spt = (double **)Mem(MEM_ALLOC, xp, sizeof(double *));	
	      
	      for(n = 0; n < xp; n++)
		spt[n] = (double *)Mem(MEM_ALLOC, yp, sizeof(double));
	    }
	  else
	    spt = NULL;

	  /*******************************************************************/

	  /***** Draw materials **********************************************/

	  /* Start parallel timer */
	  
	  StartTimer(TIMER_OMP_PARA);
	  
#ifdef OPEN_MP
#pragma omp parallel private(n, m, x, y, z, u, v, w, cell, mat, ptr, f, T, id)
#endif
	  {
	    /* Get Open MP thread id */

	    id = OMP_THREAD_NUM;
	    
	    /* Avoid compiler warnings */
	    
	    x = 0;
	    y = 0;
	    z = 0;
	    u = 0;
	    v = 0;
	    w = 0;
	    
	    /* Loop over geometry */
	    
#ifdef OPEN_MP
#pragma omp for
#endif
	    for (n = 0; n < xp; n++)
	      {
		for (m = 0; m < yp; m++)
		  {
		    /* Reset pixel value */
		    
		    matrix1[n][m] = 0;
		    matrix2[n][m] = 0;
		    
		    /* Calculate Co-ordinates */
		    
		    switch ((long)RDB[gpl + GPL_GEOM])
		      {		
		      case PLOT_MODE_YZ:
			{
			  /* yz-plot */
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    x = RDB[gpl + GPL_POS];
			  else
			    x = (xmax - xmin)/2.0 + xmin;
			  
			  y = (n/(xp - 1.0))*(ymax - ymin) + ymin;
			  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;
			  
			  u = 0.0;
			  v = 0.0;
			  w = 1.0;
			  
			  break;
			}
		      case PLOT_MODE_XZ:
			{
			  /* xz-plot */
			  
			  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    y = RDB[gpl + GPL_POS];
			  else
			    y = (ymax - ymin)/2.0 + ymin;
			  
			  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;
			  
			  u = 0.0;
			  v = 0.0;
			  w = 1.0;
			  
			  break;
			}
		      case PLOT_MODE_XY:
			{
			  
			  /* xy-plot */
			  
			  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
			  y = (m/(yp - 1.0))*(ymax - ymin) + ymin;
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    z = RDB[gpl + GPL_POS];
			  else
			    z = (zmax - zmin)/2.0 + zmin;
			  
			  u = 0.0;
			  v = 1.0;
			  w = 0.0;
			  
			  break;
			}
		      default:
			Die(FUNCTION_NAME, "Invalid plot mode");
		      }

		    if ((long)RDB[gpl + GPL_IMP_SCALE] > 0)
		      {
			/* Get importance */

			f = WWImportance(x, y, z, -1.0);

			/* Get maximum and maximum */

			fmin = RDB[gpl + GPL_IMP_MIN];
			fmax = RDB[gpl + GPL_IMP_MAX];

			/* Check */

			if (f > 0.0)
			  {
			    /* Interpolate */

			    if ((long)RDB[gpl + GPL_IMP_SCALE] == 1)
			      f = (f - fmin)/(fmax - fmin);
			    else if ((long)RDB[gpl + GPL_IMP_SCALE] == 2)
			      f = (log(f) - log(fmin))/(log(fmax) - log(fmin));
			    else
			      Die(FUNCTION_NAME, "Invalid mode");
			    
			    /* Check */

			    if (f < 0.0)
			      f = 0.0;
			    else if (f > 1.0)
			      f = 1.0;			    

			    /* Put index */

			    matrix1[n][m] = (long)(f*(ncol - nc0 - 1)) + nc0;
			  }
		      }

		    /* Find location */
		    
		    if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
		      BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w, &dummy,
					 id);
		  
		    /* Check return value */
		    
		    if (cell < 0)
		      {
			/* Put error condition */

			nerr = cell;
			
			/* Set colour */
			
			matrix1[n][m] = -cell;
			matrix2[n][m] = -cell;
		      }
		    else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
		      {
			/* Get pointer */
			
			mat = MatPtr(mat, id);
			CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
			
			/* Set colour */
			
			if (matrix1[n][m] == 0)
			  matrix1[n][m] = (long)RDB[mat + MATERIAL_COLOUR_IDX];

			matrix2[n][m] = (long)mat;
			
			/* Reset density and temperature */
			
			f = 1.0;
			T = 0.0;
			
			/* Get point from interface */
			
			IFCPoint(mat, &f, &T, id);

			/* Check for undefined region */
			
			if (f < 0.0)
			  T = -1.0;
			else if ((ptr = (long)RDB[mat + MATERIAL_PTR_IFC]) > 
				 VALID_PTR)
			  {
			    /* Scaling */
			    
			    /* JLe 26.8.2015: tää feilaa kun MATERIAL_ADENS */
			    /* ei enää rajoiteta interfacen maksimiin       */
			    /* version 2.1.24 jälkeen. */

			    /*
			    if (RDB[ptr + IFC_MAX_DENSITY] != 
				RDB[ptr + IFC_MIN_DENSITY])
			      f = (f*RDB[mat + MATERIAL_ADENS] - 
				   RDB[ptr + IFC_MIN_DENSITY])/
				(RDB[ptr + IFC_MAX_DENSITY] - 
				 RDB[ptr + IFC_MIN_DENSITY]);
			    else
			      f = 1.0;
			    */

			    if (RDB[ptr + IFC_MAX_DENSITY] == 
				RDB[ptr + IFC_MIN_DENSITY])
			      f = 1.0;
			  }
			
			/* Check temperature (overrides density) */
			
			if (T > 0.0)
			  {
			    /* Calculate factor */

			    if (RDB[mat + MATERIAL_TMS_TMAX] -
				RDB[mat + MATERIAL_TMS_TMIN] > 0.0)
			      f = 0.7*(T - RDB[mat + MATERIAL_TMS_TMIN])/
				(RDB[mat + MATERIAL_TMS_TMAX] - 
				 RDB[mat + MATERIAL_TMS_TMIN]) + 0.3;
			    else
			      f = 1.0;
			  }		  
			
			/* Check */
			
			if ((f < 0.0) || (f > 1.0))
			  matrix1[n][m] = 4;
			else 
			  matrix1[n][m] = matrix1[n][m] + 
			    (long)((1.0 - f)*(nifc - 1));
		      }
		    
		    /*********************************************************/
		    /***** Plot super-imposed source / detector region *******/
	      
#ifdef mmmmmmmm
		    
		    ptr = (long)RDB[DATA_PTR_SRC0];
		    while (ptr > VALID_PTR)
		      {
			if (InSuperCell((long)RDB[ptr + SRC_PTR_UNIV],
					(long)RDB[ptr + SRC_PTR_CELL], x, y, z)
			    == YES)
			  matrix1[n][m] = 3;
			
			/* Next */
			
			ptr = NextItem(ptr);
		      }
#endif
		    /*********************************************************/
		  }
	      }
	  }      
	  
	  /* Stop parallel timer */
	  
	  StopTimer(TIMER_OMP_PARA);
	  
	  /* Check for geometry errors */
	  
	  if ((nerr < 0) && (ini == YES))
	    {
	      if (nerr == GEOM_ERROR_NO_CELL)
		fprintf(out, "Geometry errors in plot %s (no cell).\n", fname);
	      else if (nerr == GEOM_ERROR_MULTIPLE_CELLS)
		fprintf(out, "Geometry errors in plot %s (overlap).\n", fname);
	      else
		fprintf(out, "Geometry errors in plot %s (unknown).\n", fname);
	    }
	  
	  /********************************************************************/

	  /***** Plot boundaries  *********************************************/

	  /* Start parallel timer */
	  
	  StartTimer(TIMER_OMP_PARA);
	  
#ifdef OPEN_MP
#pragma omp parallel private(n, m, x, y, z, d, cell, mat, pixw, old, ptr, id)
#endif
	  {
	    /* Get Open MP thread id */
	    
	    id = OMP_THREAD_NUM;
	    
	    /* Avoid compiler warnings */
	    
	    pixw = 0;
	    old = -1;
	    
	    /* Sweep 1 */
	    
#ifdef OPEN_MP
#pragma omp for
#endif	
	    
	    for (n = 0; n < xp; n++)
	      {
		/* Reset minimum distance */
		
		d = -1.0;
		
		for (m = 0; m < yp; m++)
		  {
		    /* Calculate Co-ordinates */
		    
		    switch ((long)RDB[gpl + GPL_GEOM])
		      {		
		      case PLOT_MODE_YZ:
			{
			  /* yz-plot */
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    x = RDB[gpl + GPL_POS];
			  else
			    x = (xmax - xmin)/2.0 + xmin;
			  
			  y = (n/(xp - 1.0))*(ymax - ymin) + ymin;
			  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;
			  
			  pixw = (zmax - zmin)/((double)yp);
			  
			  u = 0.0;
			  v = 0.0;
			  w = 1.0;
			  
			  break;
			}
		      case PLOT_MODE_XZ:
			{
			  /* xz-plot */
			  
			  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    y = RDB[gpl + GPL_POS];
			  else
			    y = (ymax - ymin)/2.0 + ymin;
			  
			  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;
			  
			  pixw = (zmax - zmin)/((double)yp);
			  
			  u = 0.0;
			  v = 0.0;
			  w = 1.0;
			  
			  break;
			}
		      case PLOT_MODE_XY:
			{		      
			  /* xy-plot */
			  
			  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
			  y = (m/(yp - 1.0))*(ymax - ymin) + ymin;
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    z = RDB[gpl + GPL_POS];
			  else
			    z = (zmax - zmin)/2.0 + zmin;
			  
			  pixw = (ymax - ymin)/((double)yp);
			  
			  u = 0.0;
			  v = 1.0;
			  w = 0.0;
		      
			  break;
			}
		      default:
			Die(FUNCTION_NAME, "Invalid plot mode");
		      }
		    
		    /* Weight window boundary */

		    if (wwp == YES)
		      if (WWDis(x, y, z, u, v, w) < pixw)
			matrix1[n][m] = 0;

		    /* Check mode */
		    
		    if ((mode == 1) || (mode == 3))
		      {
			/* Find location */
			
			if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
			  BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w, 
					     &dummy, id);

			/* Check pointer */
			
			if (cell > VALID_PTR)
			  {
			    /* Calculate distance to boundary */
			    
			    d = NearestBoundary(id);
			    
			    /* Compare minimum distance to pixel width */
			    
			    if (d < pixw)
			      {
				/* Set colour */
				
				matrix1[n][m] = 0;
			      }
			  }
		      }

		    if ((mode == 2) || (mode == 3))
		      {
			/* Material index */
			
			mat = matrix2[n][m];
			
			/* Compare to previous */
			
			if (old != mat)
			  {
			    /* Put colour */
			    /*
			    if ((old < nc0) || (mat < nc0) || 
				((old > nc0 - 1) && (mat > nc0 - 1) &&
				 (long)((double)(old - nc0)/((double)nifc)) !=
				 (long)((double)(mat - nc0)/((double)nifc))))
			    */
			    if (old != 0)
			      matrix1[n][m] = 0;
			    
			    /* Set previous pointer */
			    
			    old = mat;
			  }
		      }
		    
		    /* Draw border */
		    
		    if ((m == 0) || (m == yp - 1)) 
		      matrix1[n][m] = 0;
		  }
	      }
	    
	    /* Sweep 2 */
	    
#ifdef OPEN_MP	
#pragma omp for
#endif
	    for (m = 0; m < yp; m++)
	      {
		/* Reset minimum distance */
		
		d = -1.0;
		
		for (n = 0; n < xp; n++)
		  {
		    /* Calculate Co-ordinates */
		    
		    switch ((long)RDB[gpl + GPL_GEOM])
		      {		
		      case PLOT_MODE_YZ:
			{
			  /* yz-plot */
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    x = RDB[gpl + GPL_POS];
			  else
			    x = (xmax - xmin)/2.0 + xmin;
			  
			  y = (n/(xp - 1.0))*(ymax - ymin) + ymin;
			  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;
			  
			  pixw = (ymax - ymin)/((double)yp);
			  
			  u = 0.0;
			  v = 1.0;
			  w = 0.0;
			  
			  break;
			}
		      case PLOT_MODE_XZ:
			{
			  /* xz-plot */
			  
			  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    y = RDB[gpl + GPL_POS];
			  else
			    y = (ymax - ymin)/2.0 + ymin;
			  
			  z = ((yp - m - 1)/(yp - 1.0))*(zmax - zmin) + zmin;
			  
			  pixw = (xmax - xmin)/((double)xp);
			  
			  u = 1.0;
			  v = 0.0;
			  w = 0.0;
			  
			  break;
			}
		      case PLOT_MODE_XY:
			{
			  
			  /* xy-plot */
			  
			  x = (n/(xp - 1.0))*(xmax - xmin) + xmin;
			  y = (m/(yp - 1.0))*(ymax - ymin) + ymin;
			  
			  if (RDB[gpl + GPL_POS] > -INFTY)
			    z = RDB[gpl + GPL_POS];
			  else
			    z = (zmax - zmin)/2.0 + zmin;
			  
			  pixw = (xmax - xmin)/((double)xp);
			  
			  u = 1.0;
			  v = 0.0;
			  w = 0.0;
			  
			  break;
			}
		      default:
			Die(FUNCTION_NAME, "Invalid plot mode");
		      }
		    
		    /* Weight window boundary */

		    if (wwp == YES)
		      if (WWDis(x, y, z, u, v, w) < pixw)
			matrix1[n][m] = 0;

		    /* Check mode */
		    
		    if ((mode == 1) || (mode == 3))
		      {		
			/* Find location */
			
			if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
			  BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w, 
					     &dummy, id);
		    
			/* Check pointer */
			
			if (cell > VALID_PTR)
			  {
			    /* Calculate distance to boundary */
			    
			    d = NearestBoundary(id);
			    
			    /* Compare minimum distance to pixel width */
			    
			    if (d < pixw)
			      {
				/* Set colour */
				
				matrix1[n][m] = 0;
			      }
			  }
		      }
		    
		    if ((mode == 2) || (mode == 3))
		      {
			/* Material index */
			
			mat = matrix2[n][m];
			
			/* Compare to previous */
			
			if (old != mat)
			  {
			    /* Avoid double lines */
			    /*
			    if ((old < nc0) || (mat < nc0) ||
				((old > nc0 - 1) && (mat > nc0 - 1) &&
				 (long)((double)(old - nc0)/((double)nifc)) !=
				 (long)((double)(mat - nc0)/((double)nifc))))
			    */
			    if ((old != 0) && (n > 0) && 
				(matrix1[n - 1][m] != 0))
			      matrix1[n][m] = 0;
			    
			    /* Set previous pointer */
			    
			    old = mat;
			  }
		      }
		    
		    /* Draw border */
		    
		    if ((n == 0) || (n == xp - 1)) 
		      matrix1[n][m] = 0;
		  }
	      }
	  }
  
	  /* Stop parallel timer */
	  
	  StopTimer(TIMER_OMP_PARA);
	
	  /********************************************************************/

	  /***** Plot source point distribution *******************************/

	  /* Check pointer to source distribution */

	  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
	    {
	      /* Reset distribution */

	      for (n = 0; n < xp; n++)
		for (m = 0; m < yp; m++)
		  spt[n][m] = 0.0;

	      /* Get pointer */

	      ptr = (long)RDB[DATA_PART_PTR_SOURCE];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Pointer to first after dummy */
	      
	      ptr = NextItem(ptr);
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Loop over source */
	      
	      while (ptr > VALID_PTR)
		{
		  /* Get coordinates */

		  x = RDB[ptr + PARTICLE_X];
		  y = RDB[ptr + PARTICLE_Y];
		  z = RDB[ptr + PARTICLE_Z];

		  /* Check plot type */

		  switch ((long)RDB[gpl + GPL_GEOM])
		    {
		    case PLOT_MODE_YZ:
		      {
			/* yz-plot */
		
			n = xp*(y - ymin)/(ymax - ymin);
			m = yp*(z - zmin)/(zmax - zmin);
	
			break;
		      }
		    case PLOT_MODE_XZ:
			{
			  /* xz-plot */
			  
			  n = xp*(x - xmin)/(xmax - xmin);
			  m = yp*(z - zmin)/(zmax - zmin);

			  break;
			}
		    case PLOT_MODE_XY:
		      {
			/* yz-plot */
			
			n = xp*(x - xmin)/(xmax - xmin);
			m = yp*(y - ymin)/(ymax - ymin);

			break;
		      }
		    default:
		      Die(FUNCTION_NAME, "Invalid plot mode");
		    }
		  
		  /* Put point */

		  if ((n > -1) && (n < xp) && (m > -1) && (m < yp))
		    spt[n][m] = spt[n][m] + 1.0;
		  
		  /* Next */
		  
		  ptr = NextItem(ptr);
		}

	      /* Calculate mean */

	      fmean = 0.0;
	      i = 0;
	      
	      for (n = 0; n < xp; n++)
		for (m = 0; m < yp; m++)
		  if (spt[n][m] > 0.0)
		    {
		      fmean = fmean + spt[n][m];
		      i++;
		    }

	      if (i > 0)
		fmean = fmean/((double)i);

	      /* distribution minimum and maximum */

	      fmin = 0.5*fmean*(1.0 - RDB[DATA_SOURCE_PT_ANIM_F]);
	      fmax = 0.5*fmean*(1.0 + RDB[DATA_SOURCE_PT_ANIM_F]);	

	      /* Put colours */

	      if (fmean > 0)
		for (n = 0; n < xp; n++)
		  for (m = 0; m < yp; m++)
		    if (spt[n][m] > 0.0)
		      {
			/* Calculate factor */
		      
			f = (spt[n][m] - fmin)/(fmax - fmin);

			/* Adjust */

			if (f < 0.0)
			  f = 0.0;
			else if (f > 1.0)
			  f = 1.0;
		      
			/* Put color */
		      
			matrix1[n][m] = (long)(f*((double)ncol - 1.0)) + ncol;
		      }
	    }
	  
	  /*******************************************************************/

	  /***** Draw image **************************************************/

	  /* Check matrix */
	  
	  for (n = 0; n < xp; n++)
	    for (m = 0; m < yp; m++)
	      if ((matrix1[n][m] < 0) || (matrix1[n][m] > COLOURS - 1))
		Die(FUNCTION_NAME, "Invalid colour index %ld (%ld %ld)",
		    matrix1[n][m], n, m);
	  	  
	  /* Create image */
	  
	  im = gdImageCreate(xp, yp);
	  
	  /* Generate palette */
	  
	  for(n = 0; n < COLOURS; n++)
	    palette[n] = gdImageColorAllocate(im, R[n], G[n], B[n]);
	  
	  /* Start parallel timer */
	  
	  StartTimer(TIMER_OMP_PARA);
	  
	  /* Draw image */
	  
#ifdef OPEN_MP      
#pragma omp parallel private(n, m)
#endif
	  {
	    
#ifdef OPEN_MP
#pragma omp for
#endif
	    for (n = 0; n < xp; n++)
	      for (m = 0; m < yp; m++)
		{
		  /* Flip y-axis in xy-plot */
		  
		  if ((long)RDB[gpl +  GPL_GEOM] == PLOT_MODE_XY)
		    gdImageSetPixel(im, n, m, palette[matrix1[n][yp - m - 1]]);
		  else
		    gdImageSetPixel(im, n, m, palette[matrix1[n][m]]);   
		}
	  }
	  
	  /* Stop parallel timer */
	  
	  StopTimer(TIMER_OMP_PARA);

	  /* Plot tracks */

	  if ((long)RDB[DATA_TRACK_PLOTTER_HIS] > 0)
	    {
	      /* Plot tracks */
	      
	      fiss = PlotTracks(gpl, im, ncol, nf, xmin, xmax, ymin, ymax, 
				zmin, zmax, xp, yp);

	      /* Check flag and flash palette */

	      if (fiss == YES)
		FlashPalette(im, R, G, B, 2);
	    }

	  /*******************************************************************/

	  /***** Write in file and cleanup ***********************************/

	  /* Open file for writing */
	  
	  if ((fp = fopen(fname, "w")) == NULL)      
	    Die(FUNCTION_NAME, "Unable to open file for writing"); 
	  
	  /* Write image (png format) */
	  
	  gdImagePng(im, fp);
	  
	  /* Free image */
	  
	  gdImageDestroy(im);
	  
	  /* Free matrix */
	  
	  for(n = 0; n < xp; n++)
	    {
	      Mem(MEM_FREE, matrix1[n]);
	      Mem(MEM_FREE, matrix2[n]);
	    }

	  Mem(MEM_FREE, matrix1);
	  Mem(MEM_FREE, matrix2);

	  /* Free source point distribution */

	  if ((long)RDB[DATA_SOURCE_PT_ANIM] == YES)
	    {
	      for(n = 0; n < xp; n++)
		Mem(MEM_FREE, spt[n]);
	      
	      Mem(MEM_FREE, spt);
	    }

	  /* Close file */
	  
	  fclose(fp);
	  
	  /*******************************************************************/
      
	  /* Next plot */
	  
	  gpl = NextItem(gpl);
	}  
    }

  if (nf > 1)
    fprintf(out, " %3.0f%% complete\n\nOK.\n\n", 100.0);
	      
  /***************************************************************************/

  /* Reset plotter mode */
  
  WDB[DATA_PLOTTER_MODE] = (double)NO;

  /* Done */

  if ((long)RDB[DATA_PART_PTR_SOURCE] < VALID_PTR)
    fprintf(out, " 100%% complete\n\n"); 

  /* Check stop mode */

  if (((long)RDB[DATA_STOP_AFTER_PLOT] == STOP_AFTER_PLOT_GEOM) || 
      (ext == YES))
    exit(-1);

  /****************************************************************************/

#endif
}

/*****************************************************************************/

/***** Flash palette to indicate fission in track plotter mode ***************/

#ifndef NO_GFX_MODE

void FlashPalette(gdImagePtr im, long *r0, long *g0, long *b0, long type)
{
  long n, r, g, b;

  /* Reset palette */

  for(n = 0; n < COLOURS; n++)
    gdImageColorDeallocate(im, n);

  /* Check type */

  if (type == 1)
    {
      /***********************************************************************/

      /***** Invert colors ***************************************************/

      for(n = 0; n < COLOURS; n++)
	{
	  r = 255 - r0[n];
	  g = 255 - g0[n];
	  b = 255 - b0[n];
	  
	  gdImageColorAllocate(im, r, g, b);
	}	
      
      /***********************************************************************/
    }
  else if (type == 2)
    {
      /***********************************************************************/

      /***** Brighten colors *************************************************/

      /* First 5 are unchanged */

      for(n = 0; n < 5; n++)
	gdImageColorAllocate(im, r0[n], g0[n], b0[n]);
	
      /* Adjust remaining */

      for(n = 5; n < COLOURS; n++)
	{
	  if ((r = r0[n] + 10) > 255)
	    r = 255;
	  if ((g = g0[n] + 10) > 255)
	    g = 255;
	  if ((b = b0[n] + 10) > 255)
	    b = 255;
	  
	  gdImageColorAllocate(im, r, g, b);
	}	
      
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid palette flash type");
}

#endif

/*****************************************************************************/
