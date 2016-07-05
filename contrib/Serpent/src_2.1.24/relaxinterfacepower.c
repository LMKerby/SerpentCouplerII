/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : relaxinterfacepower.c                          */
/*                                                                           */
/* Created:       2014/07/07 (VVa)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Relaxes power solution for multi-physics interface           */
/*                                                                           */
/* Comments: -Add relaxation for interface flux                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RelaxInterfacePower:"

/*****************************************************************************/

void RelaxInterfacePower()
{
  long loc0, loc1, ptr, ptr1, n, nz, nr, i, j, k, l, na, nt, type, idx;
  long tme, tb, loc2, nst, uni;
  double alpha, tming, tmaxg, tmin, tmax;
  double sumpow, pow1, pow2, D;

  /* Check that interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  fprintf(out, "Relaxing interface powers...\n");
  
  /* Get relaxation alpha from memory */

  alpha = RDB[DATA_SOL_REL_ALPHA];

  /* Reset power counters */

  sumpow = 0.0;
  pow1 = 0.0;
  pow2 = 0.0;

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

	      /* Get pointers to power tallies */

	      ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
	      CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);

	      ptr1 = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_REL];
	      CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr1);

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

				  pow2 = pow2 + RDB[ptr1 + idx];

				  /* If this is the first calculation, then use plain value */
				  /* Otherwise relax                                        */

				  if(WDB[DATA_SOL_REL_ITER] == (double)0)				
				    WDB[ptr1 + idx] = Mean(ptr,i,j,k,l);				
				  else
				    WDB[ptr1 + idx] = RDB[ptr1 + idx] - 
				      alpha*D*(RDB[ptr1 + idx] - Mean(ptr,i,j,k,l));

				  /* Add to sum of power */

				  sumpow = sumpow + RDB[ptr1 + idx]*RDB[nst + NEST_COUNT];

				  /* Add to sum of momentary power */

				  pow1 = pow1 + Mean(ptr,i,j,k,l)*RDB[nst + NEST_COUNT];

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

  fprintf(out, "OK.\n\n");
  fprintf(out, "\nPower = %E (was: %E mom: %E) alpha %E D %E\n", sumpow, pow2, pow1, alpha, D);

  return;

  /***************************************************************************/
}

/*****************************************************************************/
