/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writefinixifc.c                                */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Writes new values from FINIX to IFC                          */
/*                                                                           */
/* Comments:  - Also calculates convergence criterion                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "UpdateFinixIFC:"

/*****************************************************************************/

void UpdateFinixIFC()
{
  long axi, fib, fpe, ifc, tme, i, j, k;
  long ang, t0ptr, t1ptr, crptr, hrptr, dfptr, tbi;
  long prev, nt, maxt, storeconv;
  double f;
  double rc, rc0, rh, rh0, T, T_prev;
  double TL2;
  Rod *rod;
  Results *results, *results_prev;
  Options *options;
  double eps, maxeps, dT, maxdT, maxT, maxTp, maxr, tming, tmaxg, tmaxnext;
  double L2abs, L2rel;
  char tmpstr[MAX_STR];
  FILE *fout;

  /*printf("Updating FINIX Interface\n");*/

  fib = (long)RDB[DATA_PTR_FIN0];

  /* Get pointer to interface */

  ifc = (long)RDB[fib + FINIX_PTR_IFC];

  /* Pointer to FUEP BLOCK */

  fpe = (long)RDB[ifc + IFC_PTR_FUEP];

  if (WDB[DATA_SOL_REL_ITER] == (double)0)
    {
      sprintf(tmpstr, "%s_Tconv%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	      (long)RDB[DATA_BURN_STEP]);

      fout = fopen(tmpstr, "w");

      fprintf(fout,"\nidx = 1;\n\n");
    }
  else
    {
      sprintf(tmpstr, "%s_Tconv%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	      (long)RDB[DATA_BURN_STEP]);

      fout = fopen(tmpstr, "a");

      fprintf(fout,"\nidx = idx + 1;\n\n");
    }

  /* Reset maximum of convergence criterion */

  maxeps = 0.0;
  maxr  = 0.0;
  maxT  = 0.0;
  maxTp = 0.0;
  maxdT = 0.0;

  L2abs = 0.0;
  L2rel = 0.0;
  TL2 = 0.0;

  /* Loop over pins */
  while (fpe > VALID_PTR)
    {

      /* Get transport time bin limits */

      tmaxg = RDB[DATA_TIME_CUT_TMAX];

      tming = RDB[DATA_TIME_CUT_TMIN];

      /************************************************/
      /* Get end time of next transport time interval */
      /************************************************/

      /* Get current time bin index */

      nt = (long)RDB[DATA_DYN_TB];

      /* Get number of time bins */

      maxt = (long)RDB[DATA_DYN_NB];

      /* Check if this is the last transport time interval */

      if (nt == maxt - 1)
	{
	  /* If this is the last transport time interval, solve only until */
	  /* end of this interval */
      
	  tmaxnext = RDB[DATA_TIME_CUT_TMAX];

	}
      else
	{
	  /* This is not the last transport time interval, solve also for */
	  /* next time interval to have a better initial guess for that   */

	  /* Get pointer to time bin structure */

	  tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];
	  CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);

	  /* Get pointer to bins */

	  tme = (long)RDB[tme + TME_PTR_BINS];
	  CheckPointer(FUNCTION_NAME, "(tme)", DATA_ARRAY, tme);    

	  /* Get time cut-off from next bin  */

	  tmaxnext = RDB[tme + nt + 2];

	}

      /* Reset flag to store convergence parameters */

      storeconv = 1;

      /* Loop over time steps */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

      while (tbi > VALID_PTR)
	{

	  /* If this IFC time bin ends before current transport time interval */
	  /* get next IFC time interval */

	  if (RDB[tbi + IFC_FUEP_T_TMAX] < tming + 1E-15)
	    {
	      /* Next time bin */
	      
	      tbi = NextItem(tbi);

	      /* Cycle loop */

	      continue;
	    }

	  /* If this IFC time bin begins after the end time of current    */
	  /* transport interval, we don't have tallied power for it yet   */
	  /* but we have calculated an initial guess based on the power   */
	  /* we have for current transport interval. Do not take that in  */
	  /* account for the convergence analysis */

	  if (RDB[tbi + IFC_FUEP_T_TMIN] > tmaxg - 1E-15)
	    storeconv = 0;

	  /* If this IFC time bin begins after the end time of next transp. bin */
	  /* We have looped over all relevant ifc bins, break out of loop       */

	  if (RDB[tbi + IFC_FUEP_T_TMIN] > tmaxnext - 1E-15)
	    break;

	  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	    {
	      /* Criticality source simulation get steady state temperature */
	      /* Get pointer to FINIX block based on fpe */

	      fib = FinixPtrFromFpe(fpe);

	      /* Get steady state Finix pointers */

	      rod = (Rod *)((long)RDB[fib + FINIX_PTR_ROD]);
	      options = (Options *)((long)RDB[fib + FINIX_PTR_OPTIONS]);
	      results = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);

	      results_prev = results;

	    }
	  else
	    {
	      /* Get Finix pointers for the time step*/
	      
	      rod = (Rod *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_ROD]);
	      results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	      options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
		
	      /* Get pointer to results at BOI */

	      if ((prev = PrevItem(tbi)) > VALID_PTR)
		{
		  results_prev = 
		    (Results *)((long)RDB[prev + IFC_FUEP_FINIX_PTR_RESULTS]);
		}
	      else
		{
		  /* This is the first step, get initial steady state temperature */
		  /* via steady state calculation*/

		  /* First timestep, get steady state temperature */
		  /* Get pointer to FINIX block based on fpe */

		  fib = FinixPtrFromFpe(fpe);

		  /* Check if found */

		  if (fib < VALID_PTR)
		    Die(FUNCTION_NAME, "Could not find linked Finix block");
	  
		  results_prev = (Results *)((long)RDB[fib + FINIX_PTR_RESULTS]);

		}
	    }
	  /* Get axial bin pointer */

	  axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];
	  CheckPointer(FUNCTION_NAME, "(axi)", DATA_ARRAY, axi);

	  /* Reset axial index */

	  i = 0;
	  
	  /* Loop over axial zones to update temperatures and hot radii */

	  while (axi > VALID_PTR)
	    {


	      if (RDB[axi + IFC_FUEP_AX_N_ANG] > 1)
		Die(FUNCTION_NAME, "Angular zones not supported yet.");

	      ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG];

	      /* Tän seuraavan vois periaatteessa korvata memcpylla jos ei ois */
	      /* noita R2:ia. Vielä simppelimpi systeemi olisi, että IFC ja FIN*/
	      /* pointterit osoittaisi samaan muistialueeseen (yksiköt)*/

	      t0ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP0];
	      t1ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP1];
	      crptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
	      hrptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];
	      dfptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF];

	      k=0;

	      /* Handle central hole */

	      if (rod->pellet_inner_radius != 0.0)
		{
		  /* Get temperature */

		  T = results->temperature[i][0];
		  T_prev = results_prev->temperature[i][0];

		  /* Get radial node positions */
		      
		  rc = results->radial_node_position_cold[i][0];
		  rh = results->radial_node_position[i][0];

		  /* Radius is zero */

		  WDB[crptr + 0 ] = 0.0;
		  WDB[hrptr + 0 ] = 0.0;

		  /* Calculate density factor of central hole       */
		  /* (Should be dependent on the other gas volumes) */

		  f = (rc*rc)/(rh*rh);

		  /* Cutoff for larger than unity density factors */

		  if (f > 1.0)
		    f = 1.0;

		  /* Store density factor */

		  WDB[dfptr + 0 ] =  f;

		  /* Calculate local convergence criterion */

		  if (storeconv)
		    {

		      if(RDB[t0ptr + 0] == 0.0)
			Die(FUNCTION_NAME, "Zero temperature");

		      /* Calculate pointwise relative convergence criterion */

		      eps = fabs((T - RDB[t1ptr + 0])/T);
		      dT = fabs(T - RDB[t1ptr + 0]);

		      /* Add to absolute L2 norm of difference */

		      L2abs = L2abs + (T - RDB[t1ptr + 0])*
			(T - RDB[t1ptr + 0]);

		      /* Add to relative L2 norm of difference */

		      L2rel = L2rel + (T - RDB[t1ptr + 0])*
			(T - RDB[t1ptr + 0])/(T*T);

		      /* Add to absolute L2 norm of temperature distribution */
		      
		      TL2 = TL2 + T*T;

		      /* Compare to maximum and store values */
		      if (eps > maxeps)
			{
			  maxeps = eps;
			  maxr = rc;
			  maxT = T;
			  maxTp = RDB[t1ptr + 0];
			}

		      /* Store new maximum absolute diference */

		      if (dT > maxdT)
			maxdT = dT;

		    }

		  /* Update centerline temperatures */

		  WDB[t0ptr  + 0] = T_prev;
		  WDB[t1ptr  + 0] = T;

		  /* Increment radial node count */

		  k = 1;
		}

	      /* Loop over remaining radial nodes */

	      for (j = 0; j < options->pellet_radial_nodes + 
		     options->clad_radial_nodes ; j++)
		{

		  /* Get temperature */

		  T = results->temperature[i][j];
		  T_prev = results_prev->temperature[i][j];

		  /* Get radial node positions */

		  rc = results->radial_node_position_cold[i][j];
		  rh = results->radial_node_position[i][j];

		  /* Get radial node positions for prev node */

		  if (j > 0)
		    {
		      rc0 = results->radial_node_position_cold[i][j-1];
		      rh0 = results->radial_node_position[i][j-1];
		    }
		  else
		    {
		      rc0 = 0.0;
		      rh0 = 0.0;
		    }

		  /* Write new cold radius */

		  WDB[crptr + k] = (rc*100)
		    *(rc*100);

		  /* Write new hot radius */

		  WDB[hrptr + k] = (rh*100)
		    *(rh*100);

		  /* Calculate density factor */
		  /* Handle pellet centerline densityfactor */

		  if (rh > 0)
		    f = (rc*rc - rc0*rc0)/
		      (rh*rh - rh0*rh0);
		  else
		    f = 1.0;

		  /* Cutoff for larger than unity density factors   */
		  /* TODO: Majorant density can be now increased in */
		  /* mat-card -> Cut-off should only come to play   */
		  /* if majorant density would be exceeded          */

		  if (f > 1.0)
		    f = 1.0;

		  /* Store density factor */

		  WDB[dfptr + k] = f;

		  /* Calculate local convergence criterion */
		  if (storeconv)
		    {
		      eps = fabs((RDB[t1ptr + k]-T)/T);
		      dT =  fabs(RDB[t1ptr + k]-T);

		      /* Add to absolute L2 norm of difference */

		      L2abs = L2abs + (RDB[t1ptr + k]-T)*
			(RDB[t1ptr + k]-T);

		      /* Add to relative L2 norm of difference */

		      L2rel = L2rel + (RDB[t1ptr + k]-T)*(RDB[t1ptr + k]-T)/(T*T);

		      /* Add to L2 norm of distribution */

		      TL2 = TL2 + T*T;

		      /* Compare to maximum and store values*/

		      if(eps > maxeps)
			{
			  maxeps = eps;
			  maxr = rc;
			  maxT = T;
			  maxTp = RDB[t1ptr + k];
			}

		      /* Store new maximum absolute diference */

		      if (dT > maxdT)
			maxdT = dT;
		    }

		  /* Update temperatures */

		  WDB[t0ptr  + k] = T_prev;
		  WDB[t1ptr  + k] = T;
	      
		  /* Increment radial node count */
		  k++;

		}

	      /* Next axial segment */

	      axi = NextItem(axi);
	      i++;

	    }

	  /* Do not handle more time intervals */

	  tbi = NextItem(tbi);
	}

      /* Next rod */

      fpe = NextItem(fpe);
    }

  L2abs = sqrt(L2abs);
  L2rel = sqrt(L2rel);
  /*
  fprintf(out, "Max Eps was %E (%f -> %f) at %E\n", maxeps, maxTp, maxT, maxr);
  fprintf(out, "L2 abs %E rel %E\n", sqrt(L2abs), sqrt(L2rel));
  fprintf(out, "%E %%Epsmax\n", maxeps);
  fprintf(out, "%E %%dTmax\n", maxdT);
  */
  fprintf(fout, "T_eps(idx) = %E;\n", maxeps);
  fprintf(fout, "T_dT(idx) = %E;\n", maxdT);
  fprintf(fout, "T_L2(idx) = %E;\n", sqrt(L2abs));
  fprintf(fout, "T_L2rel(idx) = %E;\n", sqrt(L2abs/TL2));
  fprintf(fout, "T_L2rel2(idx) = %E;\n", sqrt(L2rel));
    
  fclose(fout);

  /*printf("Done...\n");*/

}	  

#endif

/*****************************************************************************/
