/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : precursorpopcontrol.c                          */
/*                                                                           */
/* Created:       2015/09/15 (VVa)                                           */
/* Last modified: 2015/09/15 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Does population control for precursors in PSOURCE, all       */
/*              precursors are normalized for the same emission during       */
/*              upcoming time-interval                                       */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrecursorsPopControl:"

/*****************************************************************************/

void PrecursorPopControl()
{
  long ptr, part, n, m, nsrc, nbatch, mul, N, np, id, loc0;
  double wgt0, wgt2, wgt, P, emit0, emit, emit1, t0, t1, dt, lambda, aveemit;
  

  /***************************************************************************/

  if ((long)RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_NONE)
    return;

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(out, "precursorpopcontrol.c -->\n");  
#endif

  /***** Move banked precursors to source ************************************/

  /* Reset total weight and source size */

  wgt0 = 0.0;
  emit0 = 0.0;
  nsrc = 0;

  /* Get time interval limits */

  t0 = RDB[DATA_TIME_CUT_TMIN];
  t1 = RDB[DATA_TIME_CUT_TMAX];

  dt = t1 - t0;

  /* Loop over source to calculate initial source size and weight */
  /* Pointer to first item after dummy */
      
  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];

  ptr = NextItem(ptr);
      
  /* Loop over source */
      
  while (ptr > VALID_PTR)
    {
	  
      wgt = RDB[ptr + PARTICLE_WGT];
      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];

      /* Add to total weight and source size */
	  
      wgt0 += wgt;
      nsrc += 1;

      /* Add to total activity */

      emit0 += wgt*(1-exp(-lambda*dt));
      
      ptr = NextItem(ptr);
    }

  /* Normalize to an average activity */

  /* Current average activity */

  aveemit = emit0/(double)nsrc;

#ifdef DNPRINT
  fprintf(out, "Average emission was %E neutrons during interval of %E s\n",aveemit, dt);
#endif

  /***************************************************************************/
  /***** Population control **************************************************/

  /* Normalize to the wanted average activity */

  nbatch = (long)RDB[DATA_SRC_POP]*1.5;
  aveemit = emit0/(double)nbatch;

  wgt = wgt0;
  np = nsrc;

  /* Reset id */

  id = 0;

  /* Loop over source */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  ptr = NextItem(ptr);

  /* Loop over population and sample particles for duplication */

  for (n = 0; n < nsrc; n++)
    {
      /* Get decay constant */

      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];

      /* Calculate number of neutrons to emit */

      emit = RDB[ptr + PARTICLE_WGT]*(1-exp(-lambda*dt));

      /* remove initial weight from total weight */

      wgt = wgt - RDB[ptr + PARTICLE_WGT];

      if (emit > aveemit)
	{
	  /* Sample splitting */
	  /* Calculate multiplication */

	  P = emit/aveemit;
	  mul = (long)P;
	  P = P - (double)mul;

	  /* Sample multiplication */

	  if (drand48() < P)
	    N = mul;
	  else
	    N = mul - 1;

	  /* Put new weight of particle */

	  WDB[ptr + PARTICLE_WGT] = RDB[ptr + PARTICLE_WGT]*aveemit/emit;

	  /* Add new weight to total weight */

	  wgt = wgt + RDB[ptr + PARTICLE_WGT];

	  /* Create N particles */

	  for (m = 0; m < N; m++)
	    {

	      /* Duplicate neutron */
	  
	      part = DuplicateParticle(ptr, id++);

	      /* Check id */

	      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
		id = 0;

	      /* Add to weight and population size */

	      wgt = wgt + RDB[ptr + PARTICLE_WGT];
	      np = np + 1;

	      /* Put particle in source */
	      /* Will be sorted later   */
	  
	      AddItem(DATA_PART_PTR_PSOURCE, part);
	    }

	  ptr = NextItem(ptr);

	}
      else if (emit < aveemit)
	{

	  /* Get pointer and next */
	  
	  part = ptr;
	  ptr = NextItem(ptr);

	  /* Sample russian roulette */

	  if (drand48() < emit/aveemit)
	    {
	      /* Passed russian roulette */
	      /* Increase particle weight */

	      WDB[part + PARTICLE_WGT] = RDB[part + PARTICLE_WGT]*aveemit/emit;

	      /* Add new particle weight */

	      wgt = wgt + RDB[part + PARTICLE_WGT];

	    }
	  else
	    {
	      /* Failed russian roulette, remove from source */

	      RemoveItem(part);

	      /* subtract from population size */
	      
	      np = np - 1;

	      /* Put particle back to stack */

	      ToStack(part, id++);

	      /* Check id */

	      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
		id = 0;

	    }

	}
      else
	{
	  /* Put particle weight back to total weight, do nothing to particle */

	  wgt = wgt + RDB[ptr + PARTICLE_WGT];

	  /* Next particle */

	  ptr = NextItem(ptr);
	}
	
    }  

#ifdef DNPRINT
  fprintf(out, "%ld precursors in source (was %ld), %ld wanted, wgt is now %E was %E\n",
	 np, nsrc, nbatch, wgt, wgt0);  
#endif
    
  /* Calculate new average emission during next interval */

  aveemit = aveemit*wgt0/wgt;

  /* Store new average emission during next interval */
  /* Will be used in precdet as a threshold for Russian Roulette */

  WDB[loc0 + PRECDET_AVE_EMIT] = aveemit;

  /***************************************************************************/

  /***** Normalize source ****************************************************/

  /* Pointer to first item after dummy */

  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  ptr = NextItem(ptr);  

  emit1 = 0.0;

  /* Normalize weights */

  wgt2 = 0.0;

  while(ptr > VALID_PTR)
    {

      /* Scale weight */

      WDB[ptr + PARTICLE_WGT] = RDB[ptr + PARTICLE_WGT]*wgt0/wgt;
 
      /* Get decay constant */

      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];

      /* Calculate proportion of all delayed neutrons that will be emitted */
      /* from this precursor during next time interval */

      WDB[ptr + PARTICLE_U] =  RDB[ptr + PARTICLE_WGT]*(1-exp(-lambda*dt))/
	RDB[loc0 + PRECDET_W_EMIT];

      emit1 += RDB[ptr + PARTICLE_WGT]*(1-exp(-lambda*dt));

      wgt2 += RDB[ptr + PARTICLE_WGT];
      
      /* Next */
      
      ptr = NextItem(ptr);
    }

  /* Store total weight to emit */

  WDB[loc0 + PRECDET_W_EMIT] = emit1/RDB[DATA_NORM_COEF_N];

#ifdef DNPRINT
  fprintf(out, "Total weight to emit %E\n", emit1/RDB[DATA_NORM_COEF_N]);
#endif

  /* Check weight */

  if (fabs(wgt2/wgt0 - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Mismatch in weight %E %%", (wgt2/wgt0 - 1.0)*100.0);

  /* Check minimum stack size */

  if (RDB[DATA_OMP_MAX_THREADS]*RDB[DATA_PART_MIN_PSTACK]/
      RDB[DATA_PART_ALLOC_P] < 0.2)
    {
      /* Allow memory allocation */
      
      Mem(MEM_ALLOW);

      /* Allocate 20% more particles */

      AllocParticleStack(PARTICLE_TYPE_PRECURSOR, 
			 (long)(0.2*RDB[DATA_PART_ALLOC_P]));

      /* Disallow memory allocation */

      Mem(MEM_DENY);
    }

  /* Re-distribute stacks */

  /* Stacks will be redistributed at normalizedynsrc.c */
  /*
  ReDistributeStacks();
  */

#ifdef DNPRINT
  fprintf(out, "<-- precursorpopcontrol.c\n\n");  
#endif

  /***************************************************************************/
}

/*****************************************************************************/
