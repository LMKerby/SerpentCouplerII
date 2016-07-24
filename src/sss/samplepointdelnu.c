#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : samplepointdelnu.c                             */
/*                                                                           */
/* Created:       2015/09/29 (VVa)                                           */
/* Last modified: 2015/00/29 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Samples delayed neutrons from precursor concentrations at    */
/*              the beginning of an interval and adds them to que or source  */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SamplePointDelnu:"

/*****************************************************************************/

void SamplePointDelnu(long id, long np, long idx)
{
  long loc0, loc1, part;
  long gbin, seed, npmin, npmax;
  long precptr, reaptr;
  long cell, mat, rea, erg, new, idx0;
  double E, mu, u, v, w, wgt2, rnd;
  double dt, lambda;
  double t0, t1, temit, t, x, y, z;

  /***************************************************************************/

  /* Get pointer to precursor detector */

  loc0 = (long)RDB[DATA_PTR_PREC_DET];

#ifdef MPI_MODE1

  /* Check number of tasks */
  
  if (mpitasks > 1)
    {
      /* Calculate number of particles per task */

      npmax = (long)(RDB[loc0 + PRECDET_N_EMIT]/((double)mpitasks));

      /* Calculate minimum and maximum particle index for this task */

      npmin = mpiid*npmax;
      
      if (mpiid < mpitasks - 1)
	npmax = (mpiid + 1)*npmax;
      else
	npmax = (long)RDB[loc0 + PRECDET_N_EMIT];

      /* Check with task number */
      
      if ((np < npmin) || (np > npmax - 1))
	return;
    }

#endif
  
  /* Get pointer to precursor list */

  precptr = (long)RDB[loc0 + PRECDET_PTR_PREC_ARRAY];
  CheckPointer(FUNCTION_NAME, "(precptr)", DATA_ARRAY, precptr);

  /* Get pointer to reaction list */
 
  reaptr = (long)RDB[loc0 + PRECDET_PTR_REA_ARRAY];
  CheckPointer(FUNCTION_NAME, "(reaptr)", DATA_ARRAY, reaptr);

  /* Get time interval limits */

  t0 = RDB[DATA_TIME_CUT_TMIN];
  t1 = RDB[DATA_TIME_CUT_TMAX];

  /* Calculate time interval length */

  dt = t1 - t0;

  /* Calculate weight of neutrons to emit */

  wgt2 = RDB[loc0 + PRECDET_W_EMIT]/RDB[loc0 + PRECDET_N_EMIT];

  /* Init random number sequence */
  
  seed = ReInitRNG(idx);
  SEED[id*RNG_SZ] = seed;

  /* Sample precursor to emit from */
  /* This is more efficient if these are arranged by emission weight */      
  /* Now they have the same emission weight, so don't have to be sorted */

  rnd = RandF(id);

  /* Get particle to precursor list */

  part = (long)RDB[DATA_PART_PTR_PSOURCE];
      
  /* Get first item after dummy */      

  part = NextItem(part);

  /* Loop over precursors for the sample */

  while (part > VALID_PTR)
    {
      /* Compare to emission fraction */
	  
      if ((rnd = rnd - RDB[part + PARTICLE_U]) < 0.0)
	break;
	  
      /* Next precursor */
	  
      part = NextItem(part);
    }

  /* Check pointer */

  if (part < VALID_PTR)
    Die(FUNCTION_NAME, "Unable to sample precursor");

  /* Precursor was now chosen */
      
  /* Get lambda */

  lambda = RDB[part + PARTICLE_DN_LAMBDA];

  /* Sample emission time */
  /* (time from this moment to emission) */

  temit = -1/lambda*log(1 - (1 - exp(-lambda*dt))*drand48());

  /* Calculate absolute time of emission */

  t = t0 + temit;

  /* Get emission coordinates */

  x = RDB[part + PARTICLE_X];
  y = RDB[part + PARTICLE_Y];
  z = RDB[part + PARTICLE_Z];

  /* Get cell (if geometry boundary is something else than a square */
  /* we can be outside the geometry */

  u = 0.0;
  v = 0.0;
  w = 0.0;

  if ((cell = WhereAmI(x, y, z, u, v, w, id)) < VALID_PTR)
    {
      fprintf(out, " x = %E, y = %E, z = %E\n", x, y, z);

      /* If there is no cell at these coordinates something is wrong */
      /* Maybe the precursor has been read from file and is outside  */
      /* geometry due to limited accuracy when printing coordinates  */
      /* TODO: binary files? */

      Warn(FUNCTION_NAME, "Precursor is not in geometry cell?");
    }
	  
  /* Get material pointer */
	  
  if ((mat = (long)RDB[cell + CELL_PTR_MAT]) < VALID_PTR)
    {
      fprintf(out, " x = %E, y = %E, z = %E\n", x, y, z);
      fprintf(out, " Cell name %s\n", GetText(cell + CELL_PTR_NAME));

      /* If there is no material at these coordinates something is wrong  */
      /* Might be that coordinates written to file fall into the neighbor */
      /* cell due to limited output accuracy */
      /* TODO: binary files? */

      Warn(FUNCTION_NAME, "Precursor is not in material?");
    }

  /* Sample emission energy */
  /* Get precursor group bin */

  gbin = (long)RDB[part + PARTICLE_DN_GROUP];

  /* Get pointer to precursor group*/

  loc1 = (long)RDB[precptr + gbin];
  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

  /* Pointer to reaction */

  rea = (long)RDB[reaptr + gbin];
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to precursor groups energy grid */

  erg = (long)RDB[loc1 + PREC_PTR_ERG];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
	      
  /* Reset mu (scattering cosine) */
		  
  mu = 0.1;

  /* Sample energy */
  /* Using 1e-6 MeV as incoming neutron energy */
  /* Sampling is probably independent from incoming neutron */
  /* energy, but it has to be checked at some point */

  SampleENDFLaw(rea, erg, 1e-6, &E, &mu, id);

  /* Check that direction cosine was not changed */
  /* If it was changed, the delayed neutron emission reaction */
  /* actually has some kind of directional data, which would */
  /* be weird */

  if (mu != 0.1)
    Die(FUNCTION_NAME, "mu was changed by SampleENDFLaw");
		  
  /* Sample emission direction cosines isotropically */
  /* The direction for the delayed neutron should not */
  /* depend on the direction of the neutron that caused */
  /* the fission */

  IsotropicDirection(&u, &v, &w, id);

  /*
    fprintf(output, "Sampled some stuff for delayed neutron\n");
    fprintf(output, "t = %E\n", t);
    fprintf(output, "E = %E\n", E);
    fprintf(output, "(x,y,z) = (%E, %E, %E)\n", x, y, z);
    fprintf(output, "(u,v,w) = (%E %E %E)\n", u, v, w);
  */

  /* Create the neutron to be emitted and set the parameters */
  /* Duplicate incident neutron */
      
  new = FromStack(PARTICLE_TYPE_NEUTRON, id);

  /* Put variables */
      
  WDB[new + PARTICLE_X] = x;
  WDB[new + PARTICLE_Y] = y;
  WDB[new + PARTICLE_Z] = z;
      
  WDB[new + PARTICLE_U] = u;
  WDB[new + PARTICLE_V] = v;
  WDB[new + PARTICLE_W] = w;
      
  WDB[new + PARTICLE_E] = E;
  WDB[new + PARTICLE_WGT] = wgt2;
  WDB[new + PARTICLE_PTR_MAT] = (double)mat;
  WDB[new + PARTICLE_DN_GROUP] = (double)gbin;
  WDB[new + PARTICLE_DN_LAMBDA] = lambda;

  /* Store absolute emission time */
  /* PARTICLE_T0 is the birth time of the neutron*/
  /* PARTICLE_T  is the current time of the neutron */

  WDB[new + PARTICLE_T0] = t;
  WDB[new + PARTICLE_T] = t;

  /* Delayed neutron emission time */
      
  WDB[new + PARTICLE_TD] = temit;

  /* Reset thermalization time */
      
  WDB[new + PARTICLE_TT] = 0.0;      

  /* Get fission matrix index */

  idx0 = FissMtxIndex(mat, id);

  /* Put fission matrix index */
  
  WDB[new + PARTICLE_FMTX_IDX] = (double)idx0;

  /* Put MPI id */

  WDB[new + PARTICLE_MPI_ID] = (double)mpiid;

  /* Generation index from precursor */
  /* Set in precdet.c */

  WDB[new + PARTICLE_GEN_IDX] = RDB[part + PARTICLE_GEN_IDX];

  /* TODO: Set this */

  WDB[part + PARTICLE_RNG_IDX] = (double)idx;
  WDB[part + PARTICLE_HISTORY_IDX] = (double)idx;

#ifdef ScoretheseAfterBuffersAreOpen

  /* Score mesh plotter */

  ScoreMesh(new, mat, 0.0, -1.0, x, y, z, E, t, wgt, 1.0, id);


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

  /* Score initial source weight */
	  
  ptr = (long)RDB[RES_INI_SRC_WGT];
  CheckPointer(FUNCTION_NAME, "(stp)", DATA_ARRAY, stp);
  AddBuf1D(1.0, wgt, ptr, id, 0);

  /* Score source rate */

  ptr = (long)RDB[RES_TOT_NEUTRON_SRCRATE];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, wgt, ptr, id, 0);
	  	  
  /* Score source rate in fissile and non-fissile materials */
	  
  if (mat > VALID_PTR)
    {
      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
	AddBuf1D(1.0, wgt, ptr, id, 1);
      else
	AddBuf1D(1.0, wgt, ptr, id, 2);
    }

#endif	

  /* Score collision detector (for source rates) */
  /*
  ColDet(new, mat, 1.0, 0.0, x, y, z, u, v, w, E, t, wgt2, -2.0, id);
  */
  /* Score mesh plotter (for source plots) */
  /*
  ScoreMesh(new, mat, 0.0, -2.0, x, y, z, E, t, wgt2, 1.0, id);
  */
  /* Put particle in que */

  ToQue(new, id);
      
  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
