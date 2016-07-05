/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : samplesrcpoint.c                               */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2016/02/16 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Samples neutron source point                                 */
/*                                                                           */
/* Comments: - Toi cell-haku ei toimi ihan oikein sillä whereami palauttaa   */
/*             alimman cellin. Universumihakua ei oo edes tehty vielä.       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SampleSrcPoint:"
#define MAX_SRC_RESAMPLE 1000000

/*****************************************************************************/

long SampleSrcPoint(long id, long np, long idx)
{
  long n, cell, mat, src, surf, lst, loc0, loc1, rea, ptr, cell0, mat0, ncol;
  long type, stp, part, npmin, npmax, idx0, i;
  unsigned long seed;
  double rnd, mu, x, y, z, u, v, w, E, wgt, xmin, xmax, ymin, ymax, zmin, zmax;
  double t, dummy, wgt0, vol;

  /***************************************************************************/

  /***** Divide source to MPI tasks ******************************************/

#ifdef MPI_MODE1

  /* Check number of tasks */
  
  if (mpitasks > 1)
    {
      /* Calculate number of particles per task */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	npmax = (long)(RDB[DATA_CRIT_POP]/((double)mpitasks));
      else
	npmax = (long)(RDB[DATA_SRC_POP]/((double)mpitasks));

      /* Calculate minimum and maximum particle index for this task */

      npmin = mpiid*npmax;
      
      if (mpiid < mpitasks - 1)
	npmax = (mpiid + 1)*npmax;
      else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	npmax = (long)RDB[DATA_CRIT_POP];
      else
	npmax = (long)RDB[DATA_SRC_POP];

      /* Check with task number */
      
      if ((np < npmin) || (np > npmax - 1))
	return -1;
    }

#endif

  /***************************************************************************/
  
  /***** Source sampling starts here *****************************************/

  /* Reset material pointer */

  mat = -1;

  /* Init random number sequence */

  seed = ReInitRNG(idx);
  SEED[id*RNG_SZ] = seed;

  /* Restore initial seed for debugging */

  SEED0[id*RNG_SZ] = seed;

  /* Store history index for debugging */

  ptr = (long)RDB[DATA_PTR_PRIVA_HIS_IDX];
  PutPrivateData(ptr, idx, id);

  /* Init direction cosines for call to WhereAmI() */

  u = 0.0;
  v = 0.0;
  w = 1.0;

  /* Reset time */

  t = 0.0;

  /* Check if souce definition exist */

  if ((src = (long)RDB[DATA_PTR_SRC0]) < VALID_PTR)
    {
      /***********************************************************************/

      /***** Default uniform fission source **********************************/

      /* Check mode */

      if ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
	Die(FUNCTION_NAME, "No source definition");
      
      /* Loop until neutron is in fissile material */

      for (n = 0; n < MAX_SRC_RESAMPLE; n++)
	{
	  /* Sample coordinates */
		  
	  x = RandF(id)*(RDB[DATA_GEOM_MAXX] - RDB[DATA_GEOM_MINX])
	    + RDB[DATA_GEOM_MINX];
	  y = RandF(id)*(RDB[DATA_GEOM_MAXY] - RDB[DATA_GEOM_MINY])
	    + RDB[DATA_GEOM_MINY];
		  
	  if ((long)RDB[DATA_GEOM_DIM] == 3)
	    z = RandF(id)*(RDB[DATA_GEOM_MAXZ] - RDB[DATA_GEOM_MINZ])
	      + RDB[DATA_GEOM_MINZ];
	  else
	    z = 0.0;

	  /* Find location */
	  
	  if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
	    {
	      /* Apply boundary conditions */

	      BoundaryConditions(&cell, &x, &y, &z, &u, &v, &w, &dummy, id);

	      /* Check cell pointer */

	      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
		
	      /* Check if cell contains fissile material */

	      if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
		if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_FISSILE_MAT) ||
		    ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) ||
		    ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) ||
		    ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN))
		  {
		    /* Sample isotropic direction */
		    
		    IsotropicDirection(&u, &v, &w, id);

		    /* Energy from a maxwellian distribution */

		    E = MaxwellEnergy(1.2895, id);

		    /* Set weight to unity */

		    wgt = 1.0;

		    /* Break loop */

		    break;
		  }
	    }
	}
      
      /* Set type to neutron */

      type = PARTICLE_TYPE_NEUTRON;

      /* Check error */
      
      if (n == MAX_SRC_RESAMPLE)
	Error(0, "Unable to sample fission source - try explicit definition");
      
      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Source from user-defined distribution ***************************/

      /* Sample source definition */

      rnd = RandF(id);
      
      /* loop over sources */

      src = (long)RDB[DATA_PTR_SRC0];
      while (src > VALID_PTR)
	{
	  /* Compare to weight */
	  
	  if ((rnd = rnd - RDB[src + SRC_WGT]) < 0.0)
	    break;
	  
	  /* Next source definition */
	  
	  src = NextItem(src);
	}

      /* Check pointer */

      if (src < VALID_PTR)
	Die(FUNCTION_NAME, "Unable to sample source");
      
      /* Set type */

      type = (long)RDB[src + SRC_TYPE];

      /* Check */

      if ((type != PARTICLE_TYPE_NEUTRON) && (type != PARTICLE_TYPE_GAMMA))
	Die(FUNCTION_NAME, "Invalid particle type");
      
      /***********************************************************************/

      /***** Time ************************************************************/

      /* Check if interval is given */

      if (RDB[src + SRC_TMAX] > RDB[src + SRC_TMIN])
	t = RandF(id)*(RDB[src + SRC_TMAX] - RDB[src + SRC_TMIN]) +
	  RDB[src + SRC_TMIN];

      /***********************************************************************/
      
      /***** Direction *******************************************************/
      
      /* Check if isotropic (this is set to INFTY in processsources.c) */
      
      if (RDB[src + SRC_U0] > 1.0)
	IsotropicDirection(&u, &v, &w, id);
      else
	{
	  /* Set values */
	  
	  u = RDB[src + SRC_U0];
	  v = RDB[src + SRC_V0];
	  w = RDB[src + SRC_W0];
	}
      
      /***********************************************************************/
      
      /***** Energy **********************************************************/

      /* Check source energy */
      
      if ((E = RDB[src + SRC_E]) > ZERO)
	{
	  /* Check particle type */

	  if (type == PARTICLE_TYPE_NEUTRON)
	    {
	      /* Compare to upper boundary */
	      
	      if (E > RDB[DATA_NEUTRON_EMAX])
		Error(src, "Source energy above maximum (%1.2f MeV)", 
		      RDB[DATA_NEUTRON_EMAX]);
	      
	      /* Compare to lower boundary */
	      
	      if (E < RDB[DATA_NEUTRON_EMIN])
		Error(src, "Source energy below minimum (%1.2E MeV)", 
		      RDB[DATA_NEUTRON_EMIN]);
	    }
	  else
	    {
	      /* Compare to upper boundary */
	      
	      if (E > RDB[DATA_PHOTON_EMAX])
		Error(src, "Source energy above maximum (%1.2f MeV)", 
		      RDB[DATA_PHOTON_EMAX]);
	      
	      /* Compare to lower boundary */
	      
	      if (E < RDB[DATA_PHOTON_EMIN])
		Error(src, "Source energy below minimum (%1.2E MeV)", 
		      RDB[DATA_PHOTON_EMIN]);
	    }
	}

      /* Initialize weight */

      wgt0 = 1.0;
            
      /* Check type */
          
      if ((lst = (long)RDB[src + SRC_PTR_EBINS]) > VALID_PTR)
	{
	  /* Bin structure, sample non-zero weight */

	  for (i = 0; i < 10000; i++)
	    {
	      /* Bin Sample bin */
	      
	      n = (long)(RandF(id)*(ListSize(lst) - 1.0)) + 1;

	      /* Get pointer to upper bin */

	      loc1 = ListPtr(lst, n);
	      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	      /* Get weight */

	      if ((wgt0 = RDB[loc1 + SRC_EBIN_WGT]) > 0.0)
		break;
	    }
	  
	  /* Check */
	  
	  if (i == 10000)
	    Die(FUNCTION_NAME, "Unable to sample weight");
	  
	  /* Get bin pointer to lower bin */
	  
	  loc0 = ListPtr(lst, n - 1);
	  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
	  
	  /* Sample energy between boundaries */
	  
	  E = RandF(id)*(RDB[loc1 + SRC_EBIN_EMAX] - RDB[loc0 + SRC_EBIN_EMAX])
	    + RDB[loc0 + SRC_EBIN_EMAX];
	}
      else if ((rea = (long)RDB[src + SRC_PTR_REA]) > VALID_PTR)
	{
	  /* Put incident energy */
	  
	  if ((E = RDB[src + SRC_E]) < ZERO)
	    E = 1.000001*RDB[rea + REACTION_EMIN];
	  
	  /* Init cosine (-1.0, 0.0 and 1.0 case dubious error) */
	  
	  mu = 0.5;
	  
	  /* Sample energy and direction */
	  
	  SampleENDFLaw(rea, -1, E, &E, &mu, id);

	  /* Check if cosine is changed */
	  
	  if (mu != 0.5)
	    {
	      /* Sanity check for mu and direction vectors (for NAN's etc.) */
	      
	      CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
	      CheckValue(FUNCTION_NAME, "u", "", u, -1.01, 1.01);
	      CheckValue(FUNCTION_NAME, "v", "", v, -1.01, 1.01);
	      CheckValue(FUNCTION_NAME, "w", "", w, -1.01, 1.01);

	      /* Rotate */

	      AziRot(mu, &u, &v, &w, id);
	    }
	}
      else if ((E = RDB[src + SRC_E]) < ZERO)
	E = 1.0;
      
      /***********************************************************************/
	  
      /***** Position ********************************************************/

      /* Get boundaries */

      if ((xmin = RDB[src + SRC_XMIN]) == -INFTY)
	xmin = RDB[DATA_GEOM_MINX];
      
      if ((xmax = RDB[src + SRC_XMAX]) == INFTY)
	xmax = RDB[DATA_GEOM_MAXX];
      
      if ((ymin = RDB[src + SRC_YMIN]) == -INFTY)
	ymin = RDB[DATA_GEOM_MINY];
      
      if ((ymax = RDB[src + SRC_YMAX]) == INFTY)
	ymax = RDB[DATA_GEOM_MAXY];
      
      if ((zmin = RDB[src + SRC_ZMIN]) == -INFTY)
	zmin = RDB[DATA_GEOM_MINZ];
      
      if ((zmax = RDB[src + SRC_ZMAX]) == INFTY)
	zmax = RDB[DATA_GEOM_MAXZ];

      /* Calculate sampling volume */

      if ((long)RDB[DATA_GEOM_DIM] == 3)
	vol = (xmax - xmin)*(ymax - ymin)*(zmax - zmin);
      else
	vol = (xmax - xmin)*(ymax - ymin);
      
      /* Pointer to cell and material */
	  
      cell0 = (long)RDB[src + SRC_PTR_CELL];
      mat0 = (long)RDB[src + SRC_PTR_MAT];

      /* Re-sampling loop */
      
      for (n = 0; n < MAX_SRC_RESAMPLE; n++)
	{
	  /* Reset weight */
      
	  wgt = wgt0;

	  /* Add to track counter */

	  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
	  AddPrivateData(ptr, 1.0, id);

	  if ((RDB[src + SRC_X0] > -INFTY) && (RDB[src + SRC_Y0] > -INFTY) &&
	      (RDB[src + SRC_Z0] > -INFTY))
	    {
	      /* Check material and cell pointers */

	      if (mat0 > VALID_PTR)
		Error(src, "Point source not allowed with material source");
	      else if (cell0 > VALID_PTR)
		Error(src, "Point source not allowed with cell source");

	      /* Point source, set coordinates */

	      x = RDB[src + SRC_X0];
	      y = RDB[src + SRC_Y0];
	      z = RDB[src + SRC_Z0];
	    }
	  else if ((surf = (long)RDB[src + SRC_PTR_SURF]) > VALID_PTR)
	    {
	      /* Sample point on surface */
	      
	      SurfaceSrc(src, surf, &x, &y, &z, &u, &v, &w, id);
	    }
	  else if ((long)RDB[src + SRC_READ_FILE_TYPE] == 
		   SRC_FILE_TYPE_FUSION_PLASMA)
	    {
	      /* Fusion plasma source */

	      SamplePlasmaSrc(src, &x, &y, &z, &u, &v, &w, &E, &wgt, &t, id);
	    }
	  else if ((long)RDB[src + SRC_READ_PTR_FILE] > VALID_PTR)
	    {
	      /* Read parameters from file (all variables are overriden) */

	      ReadSourceFile(src, &x, &y, &z, &u, &v, &w, &E, &wgt, &t);
	    }
	  else
	    {
	      /* Sample coordinates */
	      
	      x = RandF(id)*(xmax - xmin) + xmin;
	      y = RandF(id)*(ymax - ymin) + ymin;
	      
	      if ((long)RDB[DATA_GEOM_DIM] == 3)
		z = RandF(id)*(zmax - zmin) + zmin;
	      else
		z = 0.0;
	    }

	  /* Override everything with user-defined subroutine */

	  if ((long)RDB[src + SRC_PTR_USR] > VALID_PTR)
	    UserSrc(src, &x, &y, &z, &u, &v, &w, &E, &wgt, &t, id);
	    
	  /* Check boundaries */

	  if ((x > RDB[src + SRC_XMAX]) || (x < RDB[src + SRC_XMIN]) ||
	      (y > RDB[src + SRC_YMAX]) || (y < RDB[src + SRC_YMIN]) ||
	      (z > RDB[src + SRC_ZMAX]) || (z < RDB[src + SRC_ZMIN]))
	    {
	      /* Not within boundaries */
	      
	      Die(FUNCTION_NAME, "Source point outside boundaries");
	    }

	  /* Source biasing with weight windows (NOTE: ei toimi tässä */
	  /* jos energia mukana */

	  if (WeightWindow(-1, -1, x, y, z, u, v, w, E, &wgt, t, NO, id) ==
	      TRACK_END_WCUT)
	    continue;

	  /* Check time cut-off */

	  if ((t >= RDB[DATA_TIME_CUT_TMIN]) && (t < RDB[DATA_TIME_CUT_TMAX]))
	    if ((cell = WhereAmI(x, y, z, u, v, w, id)) > VALID_PTR)
	      {				  
		/* Check if only surface or point source is used */
		  
		if ((cell0 < VALID_PTR) && (mat0 < VALID_PTR) &&
		    ((long)RDB[DATA_USE_DECAY_SRC] == NO) &&
		    ((long)RDB[cell + CELL_TYPE] != CELL_TYPE_OUTSIDE))
		  break;
		  
		/* Check cell */
		
		if (cell0 > VALID_PTR)
		  {
		    /* Get collision number */
		    
		    ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
		    ncol = (long)GetPrivateData(ptr, id);
		    
		    /* Check collision */
		    
		    if (TestValuePair(cell0 + CELL_COL_COUNT, ncol, id) > 0.0)
		      break;
		  }

		/* Get material pointer */
		
		if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
		  {		 
		    mat = MatPtr(mat, id);
		    CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
		  }
		
		/* Radioactive decay source */
		
		if ((mat > VALID_PTR) && 
		    ((long)RDB[DATA_USE_DECAY_SRC] == YES))
		  {
		    /* Check material and cell pointers */

		    if (mat0 > VALID_PTR)
		      Error(src, 
			    "Decay source not allowed with material source");
		    else if (cell0 > VALID_PTR)
		      Error(src, "Decay source not allowed with cell source");
		    
		    /* Call source routine */
		    
		    if (RadGammaSrc(src, mat, &E, &wgt, vol, id) > VALID_PTR)
		      {
			/* Score source weight */

			ptr = (long)RDB[mat + MATERIAL_SAMPLED_PHOTON_SRC];
			CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
			AddBuf1D(wgt, 1.0, ptr, id, 0);		    

			/* Break loop */

			break;
		      }					      
		  }

		/* Check material */

		if ((mat0 > VALID_PTR) && (mat > VALID_PTR))
		  {
		    /* Check match */

		    if (mat == mat0)
		      break;

		    /* Check match with parent */

		    if ((mat = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > 
			VALID_PTR)
		      if (mat == mat0)
			break;
		  }
	      }
	}

      /* Check failure */
      
      if (n == MAX_SRC_RESAMPLE)
	Error(0, "Source sampling failed because of low efficiency");
      
      /* Check cell pointer */

      if (cell < VALID_PTR)
	Die(FUNCTION_NAME, "Cell pointer is NULL");
      
      /* Get material pointer */
      
      if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
	{		 
	  mat = MatPtr(mat, id);
	  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
	}
      
      /* Score sampling efficiency */

      ptr = (long)RDB[RES_SRC_SAMPLING_EFF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, 1.0, ptr, id, 2 - type);
      AddBuf1D((double)n + 1.0, 1.0, ptr, id, 4 - type);

      /* Score mean weight */

      ptr = (long)RDB[RES_SRC_MEAN_WGT];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(wgt, 1.0, ptr, id, 2 - type);

      /***********************************************************************/ 
    }

  /* Adjust minimum and maximum energy */

  if (type == PARTICLE_TYPE_NEUTRON)
    {
      if (E < 1.0000001*RDB[DATA_NEUTRON_EMIN])
	E = 1.000001*RDB[DATA_NEUTRON_EMIN];
      else if (E > 0.999999*RDB[DATA_NEUTRON_EMAX])
	E = 0.999999*RDB[DATA_NEUTRON_EMAX];
    }
  else
    {
      if (E < 1.0000001*RDB[DATA_PHOTON_EMIN])
	E = 1.000001*RDB[DATA_PHOTON_EMIN];
      else if (E > 0.999999*RDB[DATA_PHOTON_EMAX])
	E = 0.999999*RDB[DATA_PHOTON_EMAX];
    }

  /* Check time (noi karsitaan jo tuolla ylempänä) */

  if (src > VALID_PTR)
    {
      /* Lower limit */

      if (t < RDB[DATA_TIME_CUT_TMIN])
	Error(src, "Source time %1.2E is below time cutoff %1.2E", t, 
	      RDB[DATA_TIME_CUT_TMIN]);

      /* Upper limit */

      if (t >= RDB[DATA_TIME_CUT_TMAX])
	Error(src, "Source time %1.2E is above time cutoff %1.2E", t, 
	      RDB[DATA_TIME_CUT_TMAX]);
    }

  /* Get particle from stack */

  part = FromStack(type, id);
  
  /* Put values */

  WDB[part + PARTICLE_X] = x;
  WDB[part + PARTICLE_Y] = y;
  WDB[part + PARTICLE_Z] = z;

  WDB[part + PARTICLE_U] = u;
  WDB[part + PARTICLE_V] = v;
  WDB[part + PARTICLE_W] = w;

  WDB[part + PARTICLE_E] = E;
  WDB[part + PARTICLE_WGT] = wgt;
  WDB[part + PARTICLE_T0] = t;
  WDB[part + PARTICLE_T] = t;
  WDB[part + PARTICLE_TD] = 0.0;
  WDB[part + PARTICLE_TT] = 0.0;
  WDB[part + PARTICLE_COL_IDX] = 0.0;

  WDB[part + PARTICLE_PTR_MAT] = (double)mat;
  WDB[part + PARTICLE_RNG_IDX] = (double)idx;
  WDB[part + PARTICLE_HISTORY_IDX] = (double)idx;

  /* Get fission matrix index */

  idx0 = FissMtxIndex(mat, id);

  /* Put fission matrix index */
  
  WDB[part + PARTICLE_FMTX_IDX] = (double)idx0;

  /* Score mesh plotter */

  ScoreMesh(part, mat, 0.0, -1.0, x, y, z, E, t, wgt, -1.0, id);

  /* Score source detector */
  
  SrcDet(part, mat, x, y, z, u, v, w, E, t, wgt, id);

  /* Check simulation mode */

  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) ||
      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN) ||
      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN))
    {
      /* Score source rate for fission matrix */
      
      if (idx0 > -1)
	{
	  /* Pointer to matrix */
	  
	  ptr = (long)RDB[DATA_PTR_FMTX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  ptr = (long)RDB[ptr + FMTX_PTR_SRC];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  /* Score total */
	  
	  AddBuf(wgt, 1.0, ptr, id, -1, 0, idx0);
	  
	  /* Score prompt (no delayed in source) */
	  
	  AddBuf(wgt, 1.0, ptr, id, -1, 1, idx0);
	}

      /* Check particle type */

      if (type == PARTICLE_TYPE_GAMMA)
	{
	  /* Score source rate */

	  stp = (long)RDB[RES_TOT_PHOTON_SRCRATE];  
	  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	  AddBuf1D(1.0, wgt, stp, id, 0);
	}
      else
	{
	  /* Score source rate */

	  stp = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
	  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	  AddBuf1D(1.0, wgt, stp, id, 0);
	  
	  /* Score initial source weight */
	  
	  ptr = (long)RDB[RES_INI_SRC_WGT];
	  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
	  AddBuf1D(1.0, wgt, ptr, id, 0);
	  
	  /* Score source rate in fissile and non-fissile materials */
	  
	  if (mat > VALID_PTR)
	    {
	      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
		AddBuf1D(1.0, wgt, stp, id, 1);
	      else
		AddBuf1D(1.0, wgt, stp, id, 2);
	    }
	}
	
      /* Put MPI id */

      WDB[part + PARTICLE_MPI_ID] = (double)mpiid;

      /* Reset generation index */

      WDB[part + PARTICLE_GEN_IDX] = 0.0;

      /* Put particle in que */
      
      ToQue(part, id);
    }
  else if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {
      /* Add to simulated batch size (this is different from DATA_NBATCH in */
      /* MPI mode, and may be different for each parallel task) NOTE: tän   */
      /* voisi tehdä jotenkin fiksumminkin jos jossain vaiheessa katsotaan  */
      /* että miten noi historiat jaetaan eri taskeille. */

#ifdef OPEN_MP
#pragma omp atomic
#endif
      WDB[DATA_SIMUL_BATCH_SIZE]++;

      /* Put particle to bank */
      
      ToBank(part, id);
    }
  else
    Die(FUNCTION_NAME, "Invalid simulation mode");

  /* Return pointer */

  return part;
}

/*****************************************************************************/
