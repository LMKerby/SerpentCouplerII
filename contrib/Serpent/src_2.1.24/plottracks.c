/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : plottracks.c                                   */
/*                                                                           */
/* Created:       2014/10/08 (JLe)                                           */
/* Last modified: 2015/01/06 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Plots particle tracks                                        */
/*                                                                           */
/* Comments: - Plottaa pelkät neutronit                                      */
/*                                                                           */
/*           - Toi animaatio jättää ilmeisesti aikakatkaistut historiat      */
/*             piirtämättä. Liike myös hidastuu ennen kuin neutroni kuolee   */
/*             pois. Interpolaatiokertoimen asettaminen INFTY:ksi poistaa    */
/*             trackin heti kun se kuolee, mikä näyttää vähän paremmalta.    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PlotTracks:"

/* Local subroutines */

void CoordPt(long, long *, long *, long, long, double, double, double, double, 
	     double, double, double, double, double);

/*****************************************************************************/

#ifndef NO_GFX_MODE

#define MAX_IDX 10000

long PlotTracks(long gpl, gdImagePtr im, long nf, double xmin, double xmax,
		double ymin, double ymax, double zmin, double zmax, long xp,
		long yp)
{
  long nftot, id, evn, n, m, n0, m0, part, type, i, idx, fiss, cyl;
  double t0, t, x0, y0, z0, x1, y1, z1, lmax, d, hl, t1 ,f, r;
  double R[5][MAX_IDX + 1];

  /* Cylindrical symmetry */

  cyl = NO;

  /* Get number of frames */

  if ((nftot = (long)RDB[DATA_TRACK_PLOT_FRAMES]) < 2)
    {
      /***********************************************************************/

      /***** No frames, plot entire track ************************************/

      /* Loop over threads and neutron bank */
  
      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	while ((part = FromTrkBank(PARTICLE_TYPE_NEUTRON, id)) > VALID_PTR)
	  {
	    /* Put particle back to bank */
	
	    ToBank(part, id);

	    /* Reset indexes */

	    n  = -1;
	    m  = -1;

	    /* Loop over events */
	    
	    evn = (long)RDB[part + PARTICLE_PTR_EVENTS];
	    while (evn > VALID_PTR)
	      {
		/* Get type */
		
		type = (long)RDB[evn + EVENT_TYPE];

		/* Skip all but collisions, start, leak, bc and cut-offs */

		if ((type != TRACK_END_COLL) && 
		    (type != TRACK_END_STRT) && 
		    (type != TRACK_END_LEAK) && 
		    (type != TRACK_END_TCUT) && 
		    (type != TRACK_END_ECUT) && 
		    (type != TRACK_END_WCUT) && 
		    (type != TRACK_END_BC))
		  {
		    /* Pointer to next */

		    evn = NextItem(evn);

		    /* Cycle loop */

		    continue;
		  }
		
		/* Get coordinates */
		
		x0 = RDB[evn + EVENT_X];
		y0 = RDB[evn + EVENT_Y];
		z0 = RDB[evn + EVENT_Z];
		
		/* Conversion to cylindrical coordinates */

		if (((long)RDB[gpl + GPL_TYPE] != PLOT_MODE_XY) && 
		    (cyl == YES))
		  {
		    r = sqrt(x0*x0 + y0*y0);
		    x0 = r;
		    y0 = r;
		  }

		/* Get position */
		
		n0 = n;
		m0 = m;

		CoordPt(gpl, &n, &m, xp, yp, x0, y0, z0, xmin, xmax, ymin, 
			ymax, zmin, zmax);

		/* Draw line */

		if (n0 > -1)
		  gdImageLine(im, n0, m0, n, m, 0);	    

		/* Next event */

		evn = NextItem(evn);
	      }
	  }

      /* Flush bank to track plot bank */
  
      FlushBank();

      /* Exit */

      return NO;
      
      /***********************************************************************/
    }
  
  /***************************************************************************/

  /***** Track plot animation ************************************************/

  /* Get number of frames and history lenght */

  nftot = (long)RDB[DATA_TRACK_PLOT_FRAMES];
  lmax = RDB[DATA_TRACK_PLOT_HIS_LENGTH];

  /* Calculate end-of-interval time */

  t = (RDB[DATA_TRACK_PLOT_TMAX] - RDB[DATA_TRACK_PLOT_TMIN])*
    ((double)(nftot - nf - 1)/((double)nftot - 1.0)) + 
    RDB[DATA_TRACK_PLOT_TMIN];

  /* Reset fission flag */

  fiss = NO;
  
  /* Loop over threads and neutron bank */
  
  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    while ((part = FromTrkBank(PARTICLE_TYPE_NEUTRON, id)) > VALID_PTR)
      {
	/* Put particle back to bank */
	
	ToBank(part, id);
	
	/* Reset point index */

	idx = 0;
  
	/* Reset coordinates and time */
	
	x0 = 0.0;
	y0 = 0.0;
	z0 = 0.0;
	t0 = 0.0;
  
	x1 = 0.0;
	y1 = 0.0;
	z1 = 0.0;
	t1 = 0.0;

	/* Loop over events and find last point within interval */
  
	evn = (long)RDB[part + PARTICLE_PTR_EVENTS];
	while (evn > VALID_PTR)
	  {
	    /* Get type */
      
	    type = (long)RDB[evn + EVENT_TYPE];
	    
	    /* Skip all but collisions, start, leak, bc and cut-offs */
	    
	    if ((type != TRACK_END_COLL) && (type != TRACK_END_STRT) && 
		(type != TRACK_END_LEAK) && (type != TRACK_END_TCUT) && 
		(type != TRACK_END_ECUT) && (type != TRACK_END_WCUT) && 
		(type != TRACK_END_BC))
	      {
		/* Pointer to next */
		
		evn = NextItem(evn);
		
		/* Cycle loop */
		
		continue;
	      }
	    
	    /* Update coordinates */

	    x1 = x0;
	    y1 = y0;
	    z1 = z0;
	    
	    x0 = RDB[evn + EVENT_X];
	    y0 = RDB[evn + EVENT_Y];
	    z0 = RDB[evn + EVENT_Z];
	    
	    /* Conversion to cylindrical coordinates */
	    
	    if (((long)RDB[gpl + GPL_TYPE] != PLOT_MODE_XY) && (cyl == YES))
	      {
		r = sqrt(x0*x0 + y0*y0);
		x0 = r;
		y0 = r;
	      }

	    /* Update time */

	    t1 = t0;
	    t0 = RDB[evn + EVENT_T];

	    /* Check interval */

	    if (t0 <= t)
	      break;

	    /* Next event */
	    
	    evn = NextItem(evn);
	  }

	/* Calculate interpolation factor */

	if (t0 < t1)
	  f = (t - t0)/(t1 - t0);
	else
	  f = INFTY;

	/* Interpolate last point */

	x0 = x0 + f*(x1 - x0);
	y0 = y0 + f*(y1 - y0);
	z0 = z0 + f*(z1 - z0);

	/* Store point in array */

	R[0][idx] = x0;
	R[1][idx] = y0;
	R[2][idx] = z0;
	R[3][idx] = t;

	/* Mark first point for extrapolation */

	if (f < 0.0)
	  R[4][idx] = 0.0;
	else
	  R[4][idx] = 1.0;

	/* Update index */

	idx++;

	/* Reset history length */

	hl = 0.0;
	d = -1;

	/* Loop over remaining */

	while (evn > VALID_PTR)
	  {
	    /* Get type */
      
	    type = (long)RDB[evn + EVENT_TYPE];
	    
	    /* Skip all but collisions, start, leak, bc and cut-offs */
	    
	    if ((type != TRACK_END_COLL) && (type != TRACK_END_STRT) && 
		(type != TRACK_END_LEAK) && (type != TRACK_END_TCUT) && 
		(type != TRACK_END_ECUT) && (type != TRACK_END_WCUT) && 
		(type != TRACK_END_BC))
	      {
		/* Pointer to next */
		
		evn = NextItem(evn);
		
		/* Cycle loop */
		
		continue;
	      }
	    
	    /* Update coordinates */

	    x1 = x0;
	    y1 = y0;
	    z1 = z0;
	    
	    x0 = RDB[evn + EVENT_X];
	    y0 = RDB[evn + EVENT_Y];
	    z0 = RDB[evn + EVENT_Z];

	    /* Conversion to cylindrical coordinates */

	    if (((long)RDB[gpl + GPL_TYPE] != PLOT_MODE_XY) && (cyl == YES))
	      {
		r = sqrt(x0*x0 + y0*y0);
		x0 = r;
		y0 = r;
	      }

	    /* Check for coincident point */

	    if ((fabs(x1 - x0) > 1E-2) || (fabs(y1 - y0) > 1E-2) ||
		(fabs(z1 - z0) > 1E-2))
	      {
		/* Store point in array */

		R[0][idx] = x0;
		R[1][idx] = y0;
		R[2][idx] = z0;
		R[3][idx] = t0;
		R[4][idx] = 1.0;

		/* Update index */

		idx++;

		/* Calculate distance */

		d = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) + 
			 (z1 - z0)*(z1 - z0));
		
		/* Update length */
		
		hl = hl + d;
		
		/* Check */
		
		if (hl > lmax)
		  break;
	      }

	    /* Next event */

	    evn = NextItem(evn);
	  }

	/* Interpolate first point */
	
	if ((f = (lmax - hl + d)/d) < 1.0)
	  if (idx > 1)
	    {
	      /* Get previous points */
	      
	      x0 = R[0][idx-1];
	      y0 = R[1][idx-1];
	      z0 = R[2][idx-1];
	      t0 = R[3][idx-1];
	      
	      x1 = R[0][idx-2];
	      y1 = R[1][idx-2];
	      z1 = R[2][idx-2];
	      t1 = R[3][idx-2];
	      
	      /* Calculate new point */
	      
	      R[0][idx - 1] = (x0 - x1)*f + x1;
	      R[1][idx - 1] = (y0 - y1)*f + y1;
	      R[2][idx - 1] = (z0 - z1)*f + z1;
	    }
	
	/* Reset indexes */
	
	n  = -1;
	m  = -1;
	
	/* Loop over points and draw track */
	
	for (i = 0; i < idx; i++)
	  {
	    /* Get coordinates */

	    x0 = R[0][i];
	    y0 = R[1][i];
	    z0 = R[2][i];

	    /* Update indexes */

	    n0 = n;
	    m0 = m;

	    /* Get position in pixels */

	    CoordPt(gpl, &n, &m, xp, yp, x0, y0, z0, xmin, xmax, ymin, 
		    ymax, zmin, zmax);

	    /* Draw line */

	    if ((n0 > -1) && ((long)R[4][i] == 1))
	      gdImageLine(im, n0, m0, n, m, 0);	    
	  }
      }

  /* Flush bank to track plot bank */
  
  FlushBank();
  
  /* Exit subroutine */

  return fiss;

  /***************************************************************************/
}

#endif

/*****************************************************************************/

/***** Coordinates to points *************************************************/

void CoordPt(long gpl, long *n, long *m, long xp, long yp, double x, double y, 
	     double z, double xmin, double xmax, double ymin, double ymax,
	     double zmin, double zmax)
{
  /* Avoid compiler warning */

  *n = -1;
  *m = -1;

  /* Calculate indexes */
	
  switch ((long)RDB[gpl + GPL_TYPE])
    {		
    case PLOT_MODE_YZ:
      {
	/* yz-plot */
	
	*n = (long)((y - ymin)/(ymax - ymin)*(xp - 1.0));
	*m = (long)((z - zmin)/(zmax - zmin)*(yp - 1.0));
	*m = yp - 1 - *m;
	
	break;
      }
    case PLOT_MODE_XZ:
      {
	/* xz-plot */
	
	*n = (long)((x - xmin)/(xmax - xmin)*(xp - 1.0));
	*m = (long)((z - zmin)/(zmax - zmin)*(yp - 1.0));
	*m = yp - 1 - *m;
	
	break;
      }
    case PLOT_MODE_XY:
      {
	/* xy-plot */
	
	*n = (long)((x - xmin)/(xmax - xmin)*(xp - 1.0));
	*m = (long)((y - ymin)/(ymax - ymin)*(yp - 1.0));
	*m = yp - 1 - *m;
	
	break;
      }
    default:
      Die(FUNCTION_NAME, "Invalid plot mode");
    }
}

/*****************************************************************************/
