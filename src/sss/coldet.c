/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : coldet.c                                       */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2015/07/28 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Scores collision flux detectors                              */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ColDet:"

/*****************************************************************************/

void ColDet(long part, long mat0, double flx0, double x0, double y0, 
	    double z0, double u0, double v0, double w0, double E0, double t0, 
	    double wgt0, double g0, long id)
{
  long det0, loc0, idx0, det1, rbin0, ptr, type, mt;
  double f0, val0, u, v, w;


#ifdef OLD_HIST
  
  long hst0, hst, trk, mat1, rbin1, idx1, loc1;
  double flx1, x1, y1, z1, E1, t1, wgt1, g1, val1, f1;

#endif

    /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];
  
  /* Loop over detectors */

  det0 = (long)RDB[DATA_PTR_DET0];
  while (det0 > VALID_PTR)
    {
      /* Compare particle types and check super-imposed */
      
      if ((type != (long)RDB[det0 + DET_PARTICLE]) || 
	  ((long)RDB[det0 + DET_PTR_SBINS] > VALID_PTR))
	{
	  /* Next detector */
      
	  det0 = NextItem(det0);

	  /* Cycle loop */

	  continue;
	}

      /* Get bin index */

      if (1 == 2)
	{
	  /* Tää on collision binnauksen testausta varten (JLe 28.7.2015) */
      
	  if ((idx0 = (long)RDB[part + PARTICLE_COL_IDX]) > 
	      (long)RDB[det0 + DET_N_TOT_BINS] - 1)
	    return;
	  
	  flx0 = 1.0;
	}
      else
	{
	  /* Get bin */
      
	  if ((idx0 = DetBin(det0, mat0, x0, y0, z0, E0, t0, id)) < 0)
	    {
	      /* Next detector */
	      
	      det0 = NextItem(det0);
	      
	      /* Cycle loop */
	      
	      continue;
	    }
	}

      /* Reset response index */

      rbin0 = 0;

      /* Get pointer to response functions */

      loc0 = (long)RDB[det0 + DET_PTR_RBINS];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Loop over responses */
      
      while (loc0 > VALID_PTR)
	{     
	  /* Get mt */

	  mt = (long)RDB[loc0 + DET_RBIN_MT];

	  /* Check special cases */

	  if ((mt == MT_MACRO_RECOILE) || (mt == MT_SOURCE_RATE) ||
	      (mt == MT_PHOTON_PULSE_HEIGHT))
	    {
	      /* Update response index */
		  
	      rbin0++;
	      
	      /* Next response */
	      
	      loc0 = NextItem(loc0);
	      
	      /* Cycle loop */
	      
	      continue;
	    }

	  /* Get response */
	      
	  f0 = DetResponse(det0, loc0, mat0, E0, g0, id);

	  /* Calculate value */
	  
	  if ((val0 = f0*flx0) == 0.0)
	    {
	      /* Update response index */
		  
	      rbin0++;
	      
	      /* Next response */
	      
	      loc0 = NextItem(loc0);
	      
	      /* Cycle loop */
	      
	      continue;
	    }
	  
	  /* Get direction vector */

	  u = RDB[det0 + DET_DIRVEC_U];
	  v = RDB[det0 + DET_DIRVEC_V];
	  w = RDB[det0 + DET_DIRVEC_W];

	  /* Check if non-zero */

	  if ((u != 0.0) || (v != 0.0) || (w != 0.0))
	    {
	      /* Calculate scalar product */

	      val0 = val0*(u*u0 + v*v0 + w*w0);

	      /* Check negative */

	      if (val0 < 0.0)
		val0 = 0.0;
	    }

	  /*******************************************************************/

	  /***** Store result ************************************************/

	  /* Get pointer to statistics */
	      
	  ptr = (long)RDB[det0 + DET_PTR_STAT];
	  CheckPointer(FUNCTION_NAME, "(stat)", DATA_ARRAY, ptr);

	  /* Get pointer to adjoint detector */

	  if ((det1 = (long)RDB[det0 + DET_PTR_ADJOINT]) < VALID_PTR)
	    {
	      /* Check index */

	      CheckValue(FUNCTION_NAME, "idx0", "", idx0, 0, 1000000000);

	      /* Score */

	      AddBuf(val0, wgt0, ptr, id, -1, idx0, rbin0);

	      /* Write to point to source file */
	      
	      WriteSourceFile(det0, x0, y0, z0, u0, v0, w0, E0, wgt0, t0, 
			      flx0, id);
	    }

#ifdef OLD_HIST

	  else
	    {
	      /* Pointer to particle history */
	      
	      hst0 = (long)RDB[part + PARTICLE_PTR_HIST];
		  
	      /* Loop over history */
		  
	      hst = hst0;
	      while ((hst = PrevItem(hst)) != hst0)
		{
		  /* Break at unassigned point */
		  
		  if ((wgt1 = RDB[hst + HIST_WGT]) < 0.0)
		    break;
		      
		  /* Check point type */

		  if (((trk = (long)RDB[hst + HIST_TRK]) != TRACK_END_VIRT) &&
		      (trk != TRACK_END_COLL))
		    continue;
		    
		  /* Get parameters at previous collision point */
		  
		  x1 = RDB[hst + HIST_X];
		  y1 = RDB[hst + HIST_Y];
		  z1 = RDB[hst + HIST_Z];
		  E1 = RDB[hst + HIST_E];
		  t1 = RDB[hst + HIST_T];
		  flx1 = RDB[hst + HIST_FLX];
		  wgt1 = RDB[hst + HIST_WGT];
		  mat1 = (long)RDB[hst + HIST_PTR_MAT];

		  /* Check flux, weight, energy and time */

		  CheckValue(FUNCTION_NAME, "flx1", "", flx1, 0.0, INFTY);
		  CheckValue(FUNCTION_NAME, "wgt1", "", wgt1, ZERO, INFTY);
		  CheckValue(FUNCTION_NAME, "E1", "", E1, ZERO, INFTY);
		  CheckValue(FUNCTION_NAME, "t1", "", t1, ZERO, INFTY);

		  /* Reset response index */

		  rbin1 = 0;

		  /* Get pointer to response functions */

		  loc1 = (long)RDB[det0 + DET_PTR_RBINS];
		  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

		  /* Loop over responses */

		  while (loc1 > VALID_PTR)
		    {
		      /* Get bin index */

		      if ((idx1 = DetBin(det1, mat1, x1, y1, z1, E1, t1, id)) 
			  < 0)
			break;
 
		      /* Get density factor */

		      g1 = DensityFactor(mat1, x1, y1, z1, t1, id);
		      CheckValue(FUNCTION_NAME, "g", "", g1, 0.0, 1.0);

		      /* Get response function */

		      f1 = DetResponse(det1, loc1, mat1, E1, g1, id);

		      /* Calculate value */

		      val1 = f1*flx1;

		      /* Score */

		      AddBuf(val0*val1*wgt1, wgt0, ptr, id, -1, idx0, idx1,
			     rbin0, rbin1);

		      /* Update response index */
	  
		      rbin1++;
		      
		      /* Next response */
		      
		      loc1 = NextItem(loc1);
		    }
		}
	    }

#else

	  else
	    ContribDet(part, det1, ptr, idx0, rbin0, val0, wgt0, -1.0, id);

#endif

	  /*******************************************************************/
	  
	  /* Update response index */
	  
	  rbin0++;

	  /* Next response */

	  loc0 = NextItem(loc0);
	}
      
      /* Next detector */
      
      det0 = NextItem(det0);
    }
}

/*****************************************************************************/
