/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : superdet.c                                     */
/*                                                                           */
/* Created:       2011/09/19 (JLe)                                           */
/* Last modified: 2015/06/07 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Scores super-imposed detectors                               */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SuperDet:"

/*****************************************************************************/

void SuperDet(long part, double x0, double y0, double z0, double u, double v, 
	      double w, double lmax, double E, double t, double wgt, long id)
{
  long det, loc0, idx, surf, type, ptr, np, in0, in1, norm, loc1, rbin;
  double d, l, val, cross, x, y, z, f;

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];
  while (det > VALID_PTR)
    {
      /* Compare particle types and check super-imposed */

      if (((long)RDB[part + PARTICLE_TYPE] != (long)RDB[det + DET_PARTICLE]) ||
	  ((loc0 = (long)RDB[det + DET_PTR_SBINS]) < VALID_PTR))
	{
	  /* Next detector */
      
	  det = NextItem(det);

	  /* Cycle loop */

	  continue;
	}

      /* Loop over surface bins */

      while (loc0 > VALID_PTR)
	{
	  /* Get bin index (loop over response functions becomes later) */
	  /* TODO: tohon ehkä joskus surface bin index */

	  if ((idx = DetBin(det, -1, x0, y0, z0, E, t, 0, id)) < 0)
	    {
	      /* Next bin */

	      loc0 = NextItem(loc0);

	      /* Cycle loop */

	      continue;
	    }

	  /* Check type */

	  if ((long)RDB[loc0 + DET_SBIN_TYPE] == SUPERDET_TYPE_CURRENT)
	    {
	      /***************************************************************/

	      /***** Current detector ****************************************/

	      /* Get normal */

	      norm = (long)RDB[loc0 + DET_SBIN_SURF_NORM];

	      /* Start at previous position */
	      
	      x = x0;
	      y = y0;
	      z = z0;

	      /* Reset distance */

	      l = 0.0;

	      /* Get surface pointer */
	      
	      surf = (long)RDB[loc0 + DET_SBIN_PTR_SURF];
	      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
	      
	      /* Get surface type */
	      
	      type = (long)RDB[surf + SURFACE_TYPE];
	      
	      /* Get number of parameters */
	      
	      np = (long)RDB[surf + SURFACE_N_PARAMS];
	      
	      /* Pointer to parameter list */
	      
	      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	      /* Get initial position */
	      
	      in0 = TestSurface(surf, x, y, z, NO, id);

	      /* Reset result */
      
	      val = 0.0;
      
	      /* Loop over all surfaces in track */
	      
	      while (l < lmax)
		{      
		  /* Get distance */
		  
		  d = SurfaceDistance(surf, &RDB[ptr], type, np, x, y, z, 
				      u, v, w, id);
		  
		  /* Check infinity */
		  
		  if (d == INFTY)
		    break;
		  
		  /* Extrapolate */
		  
		  d = d + EXTRAP_L;
		  
		  /* Update distance */
		  
		  l = l + d;
		  
		  if (l > lmax)
		    break;
		  
		  /* Update coordinates */
		  
		  x = x + d*u;
		  y = y + d*v;
		  z = z + d*w;
		  
		  /* Test position */
		  
		  in1 = TestSurface(surf, x, y, z, NO, id);
		  
		  /* Reset crossing */

		  cross = 0.0;

		  /* Add to number of crossings */
	  
		  if (norm == 0)
		    {
		      /* Net current */
		      
		      if ((in0 == NO) && (in1 == YES))
			cross = 1.0;
		      else if ((in0 == YES) && (in1 == NO))
			cross = -1.0;
		    }
		  else if (norm == 1)
		    {
		      /* Outward current */
		      
		      if ((in0 == YES) && (in1 == NO))
			cross = -1.0;
		    }
		  else if (norm == -1)
		    {
		      /* Inward current */
		      
		      if ((in0 == NO) && (in1 == YES))
			cross = 1.0;
		    }
		  else
		    Die(FUNCTION_NAME, "Invalid normal");
      
		  /* Add to total */

		  val = val + cross;
		  
		  /* Write to point to source file */

		  if (cross != 0.0)
		    WriteSourceFile(det, x, y, z, u, v, w, E, wgt, t, 
				    -1.0, id);
		  
		  /* Put previous position */
		  
		  in0 = in1;
		}

	      /* Check result */

	      if (val != 0.0)
		{
		  /* Get pointer to statistics */
		  
		  ptr = (long)RDB[det + DET_PTR_STAT];
		  CheckPointer(FUNCTION_NAME, "(stat)", DATA_ARRAY, ptr);
		  
		  /* Score */
		  
		  AddBuf1D(val, wgt, ptr, id, idx);
		}
	      
	      /***************************************************************/
	    }
	  else if ((long)RDB[loc0 + DET_SBIN_TYPE] == SUPERDET_TYPE_TLEFLUX)
	    {
	      /***************************************************************/

	      /***** Track-length detector within surface ********************/

	      /* Start at previous position */
	      
	      x = x0;
	      y = y0;
	      z = z0;

	      /* Reset distance */

	      l = 0.0;

	      /* Get surface pointer */
	      
	      surf = (long)RDB[loc0 + DET_SBIN_PTR_SURF];
	      CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
	      
	      /* Get surface type */
	      
	      type = (long)RDB[surf + SURFACE_TYPE];
	      
	      /* Get number of parameters */
	      
	      np = (long)RDB[surf + SURFACE_N_PARAMS];
	      
	      /* Pointer to parameter list */
	      
	      ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	      /* Check initial position */ 

	      in0 = TestSurface(surf, x, y, z, NO, id);

	      /* Reset result */
      
	      val = 0.0;

      	      /* Loop over all surfaces in track */
	      
	      while (l < lmax)
		{	  
		  /* Get distance */
		  
		  d = SurfaceDistance(surf, &RDB[ptr], type, np, x, y, z, 
				      u, v, w, id);
		  
		  /* Check infinity */
		  
		  if (d == INFTY)
		    break;
		  
		  /* When in0 == YES, the neutron is inside the surface */
		  /* and track-length estimator is scored */

		  if (in0 == YES)
		    {
		      /* Check that lmax is not exceeded */ 

		      if( l + d < lmax)
			val = val + d;
		      else
			val = val + (lmax - l);
		    }	  

		  /* Extrapolate */
		  
		  d = d + EXTRAP_L;
		  
		  /* Update distance */
		  
		  l = l + d;
		  		  
		  if (l > lmax)
		    break;
		  
		  /* Update coordinates */
		  
		  x = x + d*u;
		  y = y + d*v;
		  z = z + d*w;
		  
		  /* Test position */
		  
		  in0 = TestSurface(surf, x, y, z, NO, id);
	        }
	    
	      /* Check result */
	      
	      if (val != 0.0)
		{
		  /* Reset response index */

		  rbin = 0;

		  /* Get pointer to response functions */

		  loc1 = (long)RDB[det + DET_PTR_RBINS];
		  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		
		  /* Loop over responses */
      
		  while (loc1 > VALID_PTR)
		    {     
		      /* Get bin index */

		      if ((idx = DetBin(det, -1, x0, y0, z0, E, t, 
					rbin++, id)) < 0)
			break;

		      /* Skip recoil energy and source point types */

		      if (((long)RDB[loc1 + DET_RBIN_MT] != MT_SOURCE_RATE) &&
			  ((long)RDB[loc1 + DET_RBIN_MT] != MT_MACRO_RECOILE))
			{
			  /* Get response function */

			  f = DetResponse(det, loc1, -1, E, 1.0, id);
			  
			  /* Get pointer to statistics */
			  
			  ptr = (long)RDB[det + DET_PTR_STAT];
			  CheckPointer(FUNCTION_NAME, "(stat)", 
				       DATA_ARRAY, ptr);
		  
			  /* Score */
			  
			  AddBuf1D(val*f, wgt, ptr, id, idx);
			}

		      /* Next response */

		      loc1 = NextItem(loc1);		      
		    }
		}
	      
	      /***************************************************************/
	    }
	  else
	    Die(FUNCTION_NAME, "Invalid type");

	  /* Next surface bin */

	  loc0 = NextItem(loc0);
	}

      /* Next detector */

      det = NextItem(det);
    }
}

/*****************************************************************************/
