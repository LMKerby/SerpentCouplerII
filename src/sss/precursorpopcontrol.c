#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : precursorpopcontrol.c                          */
/*                                                                           */
/* Created:       2015/09/15 (VVa)                                           */
/* Last modified: 2016/04/04 (VVa)                                           */
/* Version:       2.1.26                                                     */
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
  long ptr, part, n, m, i, nsrc, nbatch, mul, N, np, id, loc0, gbin, ng;
  long gcount[8], *groupnums0, *groupnums1;
  double wgt0, wgt2, wgt1, wgt, P, emit0, emit, emit1, emit2, t0, t1, dt;
  double lambda, aveemit, wcount[8], acount[8], ecount[8];
  double *groupwgts0, *groupwgts1;

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

  /* Get number of delayed neutron groups */

  ng = (long)RDB[loc0 + PRECDET_NG];

  /* Allocate memory for temporary lists */

  groupwgts0  = (double *)Mem(MEM_ALLOC, ng, sizeof(double));
  groupwgts1  = (double *)Mem(MEM_ALLOC, ng, sizeof(double));
  groupnums0  = (long *)Mem(MEM_ALLOC, ng, sizeof(double));
  groupnums1  = (long *)Mem(MEM_ALLOC, ng, sizeof(double));

  /* Loop over source to calculate initial source size and weight */
  /* Pointer to first item after dummy */
      
  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];

  ptr = NextItem(ptr);
      
  /* Loop over source */
      
  while (ptr > VALID_PTR)
    {

      /* Get particle weight */

      wgt = RDB[ptr + PARTICLE_WGT];

      /* Get delayed neutron group bin */

      gbin = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Get decay constant */

      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];

      /* Add to total weight and source size */
	  
      wgt0 += wgt;
      nsrc += 1;

      /* Add to initial group weight */

      groupwgts0[gbin] += wgt;

      /* Add to initial group population */

      groupnums0[gbin] += 1;

      /* Add to total activity */

      emit0 += wgt*(1-exp(-lambda*dt));
      
      /* Next particle */

      ptr = NextItem(ptr);
    }

  /* Normalize to average emission on upcoming time interval */

  /* Current average emission */

  aveemit = emit0/(double)nsrc;

#ifdef DNPRINT
  fprintf(out, "Average emission was %E neutrons during interval of %E s\n",aveemit, dt);
#endif

  /***************************************************************************/
  /***** Population control **************************************************/

  /* Normalize to the wanted average emission */

  /*
  nbatch = (long)RDB[DATA_SRC_POP]*(long)RDB[DATA_PREC_SRC_FACT];
  */

  nbatch = (long)RDB[DATA_SRC_POP]*2;

  /* Average emission for the wanted number of neutrons */

  aveemit = emit0/(double)nbatch;

  /* Store initial emission */

  emit1 = emit0;

  /* Store initial weight */

  wgt1 = wgt0;

  /* Store initial number of precursors */

  np = nsrc;

  /* Store current group weights and populations*/

  for (n = 0; n < ng; n++)
    {
      groupwgts1[n] = groupwgts0[n];
      groupnums1[n] = groupnums0[n];
    }

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

      /* Get DN group */

      gbin = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Calculate number of neutrons to emit */

      emit = RDB[ptr + PARTICLE_WGT]*(1-exp(-lambda*dt));

      /* remove initial weight from total weight */

      wgt1 = wgt1 - RDB[ptr + PARTICLE_WGT];

      /* remove initial weight from group weight */

      groupwgts1[gbin] -= RDB[ptr + PARTICLE_WGT];

      /* Remove from group population */

      groupnums1[gbin] -= 1;

      /* remove initial emission from total emission */

      emit1 = emit1 - emit;

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

	  wgt = aveemit/emit*RDB[ptr + PARTICLE_WGT];

	  WDB[ptr + PARTICLE_WGT] = wgt;

	  /* Add new weight to total weight */

	  wgt1 += wgt;

	  /* Add new weight to group weight */

	  groupwgts1[gbin] += wgt;

	  /* Add to group population */

	  groupnums1[gbin] += 1;

          /* Add new emission to total emission */

          emit1 += wgt*(1-exp(-lambda*dt));

	  /* Create N particles */

	  for (m = 0; m < N; m++)
	    {

	      /* Duplicate neutron */
	  
	      part = DuplicateParticle(ptr, id++);

	      /* Check id */

	      if (id > (long)RDB[DATA_OMP_MAX_THREADS] - 1)
		id = 0;

	      /* Add to weight and population size */

	      wgt1 += wgt;

	      /* Add new weight to group weight */
	      
	      groupwgts1[gbin] += wgt;

	      /* Add to group population */

	      groupnums1[gbin] += 1;

	      /* Add new emission to total emission */

	      emit1 += wgt*(1-exp(-lambda*dt));

	      np = np + 1;

	      /* Put particle in source */
	      /* Will be sorted later   */
	  
	      AddItem(DATA_PART_PTR_PSOURCE, part);
	    }

	  ptr = NextItem(ptr);

	}
      /* Never remove the last precursor from a group, because the weights */
      /* are scaled separately for each group and this would lead to lost  */
      /* weight as the weight in this group could not be conserved */
      else if ((emit < aveemit) && (groupnums1[gbin] != 0))
	{
	  /* Get pointer and next */
	  
	  part = ptr;
	  ptr = NextItem(ptr);

	  /* Sample russian roulette */

	  if (drand48() < emit/aveemit)
	    {
	      /* Passed russian roulette */
	      /* Increase particle weight */
	      wgt = RDB[part + PARTICLE_WGT]*aveemit/emit;

	      /* Store new particle weight */

	      WDB[part + PARTICLE_WGT] = wgt;

	      /* Add to weight and population size */

	      wgt1 += wgt;

	      /* Add new weight to group weight */
	      
	      groupwgts1[gbin] += wgt;

	      /* Add to group population */

	      groupnums1[gbin] += 1;

	      /* Add new emission to total emission */

	      emit1 += wgt*(1-exp(-lambda*dt));

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
	  /* Put particle back as is */

	  wgt = RDB[ptr + PARTICLE_WGT];

	  /* Add to weight and population size */

	  wgt1 += wgt;

	  /* Add new weight to group weight */
	      
	  groupwgts1[gbin] += wgt;

	  /* Add to group population */

	  groupnums1[gbin] += 1;

	  /* Add new emission to total emission */

	  emit1 += wgt*(1-exp(-lambda*dt));

	  /* Next particle */

	  ptr = NextItem(ptr);
	}
	
    }  


#ifdef DNPRINT

  fprintf(out, "%ld precursors in source (was %ld), %ld wanted, emit is now %E was %E\n",
	  np, nsrc, nbatch, emit1, emit0);  
    
  fprintf(out, "Wgt0 %E, wgt1 %E\n", wgt0, wgt1);

  for (n = 0; n < ng; n++)
    {
      if (groupwgts1[n] < 0)
	{
	  groupwgts1[n] = 0;
	  fprintf(out, "group %ld, num0 %ld, num1 %ld, wgt0 %E, wgt1 %E, cannot scale\n", n, groupnums0[n], groupnums1[n], groupwgts0[n], groupwgts1[n]);
	}
      else
	fprintf(out, "group %ld, num0 %ld, num1 %ld, wgt0 %E, wgt1 %E, should scale by %f\n", n, groupnums0[n], groupnums1[n], groupwgts0[n], groupwgts1[n], groupwgts0[n]/groupwgts1[n]);
    }

#endif

  /***************************************************************************/

  /* Reset counters, these are used just to print out some stuff */
  /* Should be removed from the final implementation or maybe    */
  /* turned into some RES-variables or smth. */
      
  for (i = 0 ; i < 8 ; i++)
    {
      gcount[i] = 0;
      wcount[i] = 0;
      ecount[i] = 0;
      acount[i] = 0;
    }

  /***** Normalize source ****************************************************/

  /* Pointer to first item after dummy */
  
  ptr = (long)RDB[DATA_PART_PTR_PSOURCE];
  ptr = NextItem(ptr);  
  
  emit2 = 0.0;
  
  /* Normalize weights */
  
  wgt2 = 0.0;
  
  while(ptr > VALID_PTR)
    {
      /* Get decay constant */
      
      lambda = RDB[ptr + PARTICLE_DN_LAMBDA];
      
      /* Get DN group */

      gbin = (long)RDB[ptr + PARTICLE_DN_GROUP];

      /* Scale weight to conserve weight in each group */

      if (groupwgts1[gbin] <= 0)
	{
	  wgt = 0;
	  Die(FUNCTION_NAME, "Invalid group weight");
	}
      else
	wgt = RDB[ptr + PARTICLE_WGT]*groupwgts0[gbin]/groupwgts1[gbin];
      
      /* Put new weight */

      WDB[ptr + PARTICLE_WGT] = wgt;
      
      /* Calculate proportion of total emission by this precursor */
      /* Total emission should be conserved also */

      WDB[ptr + PARTICLE_U] =  wgt*(1-exp(-lambda*dt))/emit0;

      /* Add to total emission */

      emit2 += wgt*(1-exp(-lambda*dt));
      
      /* Add to total weight */

      wgt2 += wgt;

      /* Score some counters (can be removed at some point) */

      gcount[gbin]++;
      wcount[gbin] += wgt;
      ecount[gbin] += wgt*(1-exp(-lambda*dt));
      acount[gbin] += wgt*lambda;
      
      /* Next precursor */
      
      ptr = NextItem(ptr);
    }

  /* Store total weight to emit */

  WDB[loc0 + PRECDET_W_EMIT] = emit2/RDB[DATA_NORM_COEF_N];

  /* Calculate new average emission during next interval */

  aveemit = emit2/(double)np;

  /* Store new average emission during next interval */
  /* Will be used in precdet as a threshold for Russian Roulette */

  WDB[loc0 + PRECDET_AVE_EMIT] = aveemit;

  /* Free temporary lists */

  Mem(MEM_FREE, groupwgts0);
  Mem(MEM_FREE, groupwgts1);
  Mem(MEM_FREE, groupnums0);
  Mem(MEM_FREE, groupnums1);

#ifdef DNPRINT

  fprintf(out, "Emit0 %E, emit2 %E\n", emit0, emit2);

  fprintf(out, "Precursor distribution in list:\n");

  for (i = 0; i < 8; i++)
    fprintf(out, "%2ld ", gcount[i]);

  fprintf(out, "%% num\n");

  fprintf(out,"Weight distribution in list:\n");

  for (i = 0; i < 8; i++)
    fprintf(out, "%6f ", wcount[i]);

  fprintf(out,"%% wgt\n");

  fprintf(out,"Activity in list:\n");

  for (i = 0; i < 8; i++)
    fprintf(out, "%6f ", acount[i]);

  fprintf(out,"%% act\n");

  fprintf(out,"Emission in list:\n");

  for (i = 0; i < 8; i++)
    fprintf(out, "%6f ", ecount[i]);

  fprintf(out,"%% emit\n");

  fprintf(out, "Neutrons to emit %E\n", emit2);
  fprintf(out, "Total weight to emit %E\n", emit2/RDB[DATA_NORM_COEF_N]);
#endif

  /* Check weight */

  if (fabs(wgt2/wgt0 - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Mismatch in weight %E %%", (wgt2/wgt0 - 1.0)*100.0);

  /* Check emission */

  if (fabs(emit2/emit0 - 1.0) > 1E-6)
    Die(FUNCTION_NAME, "Mismatch in emission %E %%", (emit2/emit0 - 1.0)*100.0);

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
#ifdef __cplusplus 
} 
#endif 
