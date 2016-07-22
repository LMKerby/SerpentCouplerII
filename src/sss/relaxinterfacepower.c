/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : relaxinterfacepower.c                          */
/*                                                                           */
/* Created:       2014/07/07 (VVa)                                           */
/* Last modified: 2016/04/04 (VVa)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Relaxes power solution for multi-physics interface           */
/*                                                                           */
/* Comments: -Add relaxation for interface flux                              */
/*           -Add convergence output for all interface types                 */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RelaxInterfacePower:"

/*****************************************************************************/

void RelaxInterfacePower()
{
  long loc0, loc1, ptr, ptr1, ptr2, n, nz, nr, i, j, k, l, na, nt, type, idx;
  long tme, tb, loc2, nst, uni;
  double alpha, tming, tmaxg, tmin, tmax;
  double sumpow, pow1, pow2, D, mult;
  double momnew, momold, relnew, relold, PL2, rel_PL2;
  double eps, rel_eps, maxeps, rel_maxeps, maxdP, rel_maxdP, dP, rel_dP;
  double L2abs, rel_L2abs, L2rel, rel_L2rel;
  char tmpstr[MAX_STR];
  FILE *fout;

  /* Check that interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  if (WDB[DATA_SOL_REL_ITER] == (double)0)
    {
      sprintf(tmpstr, "%s_Pconv%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	      (long)RDB[DATA_BURN_STEP]);

      fout = fopen(tmpstr, "w");

      fprintf(fout,"\nidx = 1;\n\n");
    }
  else
    {
      sprintf(tmpstr, "%s_Pconv%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	      (long)RDB[DATA_BURN_STEP]);

      fout = fopen(tmpstr, "a");

      fprintf(fout,"\nidx = idx + 1;\n\n");
    }

  fprintf(out, "Relaxing interface powers...\n");
  
  /* Get relaxation alpha from memory */

  alpha = RDB[DATA_SOL_REL_ALPHA];

  /* Reset power counters */

  sumpow = 0.0;
  pow1 = 0.0;
  pow2 = 0.0;

  /* Reset maximum of convergence criterion */

  maxeps = 0.0;
  maxdP = 0.0;
  L2abs = 0.0;
  L2rel = 0.0;
  PL2 = 0.0;

  rel_maxeps = 0.0;
  rel_maxdP = 0.0;
  rel_L2abs = 0.0;
  rel_L2rel = 0.0;
  rel_PL2 = 0.0;

  /* Get relaxation factor (< 1.0 underrelaxes) */

  /* Check if no relaxation is wanted */

  if((D = RDB[DATA_SOL_REL_FACT]) == 0.0)
    D = 1;

  /* Get time bin index */
  if(RDB[DATA_SIMULATION_MODE] != (double)SIMULATION_MODE_CRIT)
    {
      tb = (long)RDB[DATA_DYN_TB];

      /* Get pointer to time bin structure */

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);
      
      /* Get pointer to bins */

      tme = (long)RDB[tme + TME_PTR_BINS];
      CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);    

      /* Get transport time interval */

      tming = RDB[tme + tb];

      tmaxg = RDB[tme + tb + 1];
    }
  else
    {
      tming = -INFTY;
      tmaxg =  INFTY;
    }
  /* Loop over interfaces */

  while(loc0 > VALID_PTR)
    {

      /* Cycle loop if no output is requested */

      if (RDB[loc0 + IFC_CALC_OUTPUT] == (double)NO)
	{
	  /* Next interface */
	  
	  loc0 = NextItem(loc0);

	  /* Cycle loop*/

	  continue;
	}

      /* Get interface type */

      type = (long)RDB[loc0 + IFC_TYPE];

      /* Only handle fuel behavior interfaces atm. */

      if((type == IFC_TYPE_FUEP) || (type == IFC_TYPE_FPIP))
	{

	  /* Get pointer to first rod */

	  loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];

	  /* Loop over pins to relax power */

	  while(loc1 > VALID_PTR)
	    {
	      /* Pointer to universe */

	      uni = (long)RDB[loc1 + IFC_FUEP_PTR_UNI];
	      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

	      /* Pointer to nest */

	      nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
	      CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

	      /* Get pointer to output limits */

	      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
	      CheckPointer(FUNCTION_NAME, "(limptr)", DATA_ARRAY, ptr);	  

	      /* get number of bins*/

	      nz = (long)RDB[ptr + FUEP_NZ];
	  
	      na = (long)RDB[ptr + FUEP_NA];
      
	      nr = (long)RDB[ptr + FUEP_NR];

	      nt = (long)RDB[ptr + FUEP_NT];

	      /* Get number of nests from nest count or FINIX definition */

	      mult = RDB[nst + NEST_COUNT];

	      if ((ptr = (long)RDB[loc1 + IFC_FUEP_PTR_FINIX]) > VALID_PTR)
		if (RDB[ptr + FINIX_N_RODS] > 0.0)
		  mult = RDB[ptr + FINIX_N_RODS];

	      /* Get pointers to power tallies */

	      ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
	      CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

	      ptr1 = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_REL];
	      CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

	      ptr2 = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_PREV];
	      CheckPointer(FUNCTION_NAME, "(Pptr2)", DATA_ARRAY, ptr2);	      

	      for(i=0; i < nz; i++)
		{
		  for(j=0; j < na; j++)
		    {
		      for(k=0; k < nr; k++)
			{
			  for(l=0; l < nt; l++)
			    {


			      /* Find time bin */

			      loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_TB];
			      CheckPointer(FUNCTION_NAME, "(TBptr)", DATA_ARRAY, ptr);

			      tmin = RDB[loc2 + l];

			      tmax = RDB[loc2 + l + 1];

			      if((tmin >= tming-1E-15) && (tmax <= tmaxg+1E-15))
				{
				  /* Calculate index for relaxed tally */

				  idx = i + j*nz + k*nz*na + l*nz*na*nr;

				  /* Add to sum of relaxed power */

				  pow2 = pow2 + RDB[ptr1 + idx]*mult;

				  momnew = Mean(ptr, i, j, k, l);
				  momold = RDB[ptr2 + idx];

				  relold = RDB[ptr1 + idx];
				  relnew = relold - alpha*D*(relold - momnew);

				  /********************/
				  /* Store new values */
				  /********************/

				  /* Store new momentary power */

				  WDB[ptr2 + idx] = momnew;

				  /* Store new relaxed power */

				  /* If this is the first calculation, then use plain value */
				  /* Otherwise relax                                        */

				  if (WDB[DATA_SOL_REL_ITER] == (double)0)
				    WDB[ptr1 + idx] = Mean(ptr,i,j,k,l);				
				  else
				    WDB[ptr1 + idx] = RDB[ptr1 + idx] - 
				      alpha*D*(RDB[ptr1 + idx] - Mean(ptr,i,j,k,l));


				  /**********************************************************/
				  /* Convergence criterions based on momentary distribution */
				  /**********************************************************/

				  /* Calculate pointwise relative convergence criterion */

				  eps = fabs(momnew - momold)/momnew;
				  dP = fabs(momnew - momold);

				  if (eps > maxeps)
				    maxeps = eps;
				  
				  if (dP > maxdP)
				    maxdP = dP;

				  /* This is the 2-norm of the momentary power distribution */

				  PL2 = PL2 + momnew*momnew;

				  /* Add to absolute L2 norm of difference */

				  L2abs = L2abs + (momnew - momold)*(momnew - momold);

				  /* Add to relative L2 norm of difference */

				  L2rel = L2rel + (momnew - momold)*(momnew - momold)/
				    (momnew*momnew);

				  /********************************************************/
				  /* Convergence criterions based on relaxed distribution */
				  /********************************************************/

				  /* Calculate pointwise relative convergence criterion */

				  rel_eps = fabs(relnew - relold)/relnew;
				  rel_dP = fabs(relnew - relold);

				  if (rel_eps > rel_maxeps)
				    rel_maxeps = rel_eps;
				  
				  if (rel_dP > rel_maxdP)
				    rel_maxdP = rel_dP;

				  /* This is the 2-norm of the new power distribution */

				  rel_PL2 = rel_PL2 + relnew*relnew;

				  /* Add to absolute L2 norm of difference */

				  rel_L2abs = rel_L2abs + (relnew - relold)*(relnew - relold);

				  /* Add to relative L2 norm of difference */

				  rel_L2rel = rel_L2rel + (relnew - relold)*(relnew - relold)/
				    (relnew*relnew);

				  /*************************************/
				  /* Add to sums of power or debugging */
				  /*************************************/

				  /* Add to sum of power */

				  sumpow = sumpow + RDB[ptr1 + idx]*mult;

				  /* Add to sum of momentary power */

				  pow1 = pow1 + Mean(ptr,i,j,k,l)*mult;

				}

			    }

			}

		    }

		}

	      /* Next pin */

	      loc1 = NextItem(loc1);
	    }
	}

      else if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH)
	{

	  /* Get pointers to power tallies */

	  ptr = (long)RDB[loc0 + IFC_PTR_STAT];
	  CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

	  ptr1 = (long)RDB[loc0 + IFC_PTR_STAT_REL];
	  CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

	  /* Get number of regions */

	  n = (long)RDB[loc0 + IFC_STAT_NREG];

	  for (i = 0; i < n; i++)
	    {

	      /* Add to sum of relaxed power */

	      pow2 = pow2 + RDB[ptr1 + i];

	      /* If this is the first calculation, then use plain value */
	      /* Otherwise relax                                        */

	      if(WDB[DATA_SOL_REL_ITER] == (double)0)				
		WDB[ptr1 + i] = Mean(ptr, i);
	      else
		WDB[ptr1 + i] = RDB[ptr1 + i] - 
		  alpha*D*(RDB[ptr1 + i] - Mean(ptr, i));

	      /* Add to sum of power */

	      sumpow = sumpow + RDB[ptr1 + i];

	      /* Add to sum of momentary power */

	      pow1 = pow1 + Mean(ptr,i);

	    }

	}
      else
	{

	  /* Get pointers to power tallies */

	  ptr = (long)RDB[loc0 + IFC_PTR_STAT];
	  CheckPointer(FUNCTION_NAME, "(Pptr_other)", DATA_ARRAY, ptr);

	  ptr1 = (long)RDB[loc0 + IFC_PTR_STAT_REL];
	  CheckPointer(FUNCTION_NAME, "(Pptr1_other)", DATA_ARRAY, ptr1);

	  /* Get number of regions */

	  n = (long)RDB[loc0 + IFC_STAT_NREG];

	  /* Get number of regions */

	  nz = (long)RDB[loc0 + IFC_NZ];

	  /* Get number of regions */

	  nr = (long)RDB[loc0 + IFC_NR];

	  for (i = 0; i < n; i++)
	    {

	      for (j = 0; j < nz; j++)
		{

		  for (k = 0; k < nr; k++)
		    {

		      /* Calculate index for relaxed tally */

		      idx = i + j*n + k*n*nz;		      

		      /* If this is the first calculation, then use plain value */
		      /* Otherwise relax                                        */

		      pow2 = pow2 + RDB[ptr1 + idx];

		      if(WDB[DATA_SOL_REL_ITER] == (double)0)				
			WDB[ptr1 + idx] = Mean(ptr, i, j, k);
		      else
			WDB[ptr1 + idx] = RDB[ptr1 + idx] - 
			  alpha*D*(RDB[ptr1 + idx] - Mean(ptr, i, j, k));


		      /* Add to sum of power */

		      sumpow = sumpow + RDB[ptr1 + idx];

		      /* Add to sum of momentary power */

		      pow1 = pow1 + Mean(ptr,i,j,k);
		    }

		}

	    }

	}
	
      /* Next interface */

      loc0 = NextItem(loc0);
    }

  fprintf(fout, "P_eps(idx) = %E;\n", maxeps);
  fprintf(fout, "P_delta(idx) = %E;\n", maxdP);
  fprintf(fout, "P_L2(idx) = %E;\n", sqrt(L2abs));
  fprintf(fout, "P_L2rel(idx) = %E;\n", sqrt(L2abs/PL2));
  fprintf(fout, "P_L2rel2(idx) = %E;\n", sqrt(L2rel));
  
  fprintf(fout, "relaxedP_eps(idx) = %E;\n", rel_maxeps);
  fprintf(fout, "relaxedP_delta(idx) = %E;\n", rel_maxdP);
  fprintf(fout, "relaxedP_L2(idx) = %E;\n", sqrt(rel_L2abs));
  fprintf(fout, "relaxedP_L2rel(idx) = %E;\n", sqrt(rel_L2abs/rel_PL2));
  fprintf(fout, "relaxedP_L2rel2(idx) = %E;\n", sqrt(rel_L2rel));
  
  fclose(fout);

  fprintf(out, "OK.\n\n");
  /*
  fprintf(out, "\nPower = %E (was: %E mom: %E) alpha %E D %E\n", sumpow, pow2, pow1, alpha, D);
  */
  return;

  /***************************************************************************/
}

/*****************************************************************************/
