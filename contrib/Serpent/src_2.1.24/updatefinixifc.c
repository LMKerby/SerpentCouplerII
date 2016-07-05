/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writefinixifc.c                                */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Writes new values from FINIX to IFC                          */
/*                                                                           */
/* Comments:  - Also calculates convergence criterion                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "UpdateFinixIFC:"

/*****************************************************************************/

void UpdateFinixIFC()
{
  long axi, fib, ptr, fpe, ifc, i, j, k;
  long ang, t0ptr, t1ptr, crptr, hrptr, dfptr, tbi;
  long prev;
  long found, nu, redo;
  double f;
  int *nnodes;
  double **T;
  double **T0;
  double **r_hot;
  double **r_cold;
  double eps, maxeps, dT, maxdT, maxT, maxTp, maxr, tming, tmaxg;
  double L2abs, L2rel;

  /* Get transport time bin limits */

  tmaxg = RDB[DATA_TIME_CUT_TMAX];

  tming = RDB[DATA_TIME_CUT_TMIN];

  /*printf("Updating FINIX Interface\n");*/

  fib = (long)RDB[DATA_PTR_FIN0];

  /* Get pointer to interface */

  ifc = (long)RDB[fib + FINIX_PTR_IFC];

  /* Pointer to FUEP BLOCK */

  fpe = (long)RDB[ifc + IFC_PTR_FUEP];

  /* Reset redo flag */

  redo = 0;

  /* Reset maximum of convergence criterion */

  maxeps = 0.0;
  maxr  = 0.0;
  maxT  = 0.0;
  maxTp = 0.0;
  maxdT = 0.0;

  L2abs = 0.0;

  L2rel = 0.0;

  /* Loop over pins */
  while(fpe > VALID_PTR)
    {

      if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	{
	  /* Steady state simulation */

	  /* Get steady state temperature */

	  fib = (long)RDB[DATA_PTR_FIN0];
	  CheckPointer(FUNCTION_NAME, "(fib)", DATA_ARRAY, fib);

	  /* Find correct pin */

	  while(fib > VALID_PTR)
	    {
	      /* Get Finix pin uni */

	      found=0;

	      /* Pointer to number of rod segments */

	      nu = WDB[fpe + IFC_FUEP_N_UNI];

	      /* Loop over rod segments */

	      ptr = (long)RDB[fpe + IFC_FUEP_PTR_UNI_LIST];

	      for(i=0; i < nu; i++)
		{
		  if(CompareStr(ptr + i, fib + FINIX_PTR_UNI_NAME))
		    found=1;
		}

	      if(found==1)
		break;
	      else
		fib = NextItem(fib);	 

	    }

	  /* Check if found */
	  if(fib < VALID_PTR)
	    Die(FUNCTION_NAME, "Could not find universe");

	  /* Get steady state distributions */

	  nnodes   = (int*)((long)RDB[fib + FINIX_PTR_NNODES]);
	  T        = (double**)((long)RDB[fib + FINIX_PTR_T]);
	  T0       = (double**)((long)RDB[fib + FINIX_PTR_T]);
	  r_hot    = (double**)((long)RDB[fib + FINIX_PTR_R]);
	  r_cold   = (double**)((long)RDB[fib + FINIX_PTR_R_COLD]);

	  /* Get time bin pointer */

	  tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

	  /* Get axial bin pointer */

	  axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];
	  CheckPointer(FUNCTION_NAME, "(axi)", DATA_ARRAY, axi);

	  /* Reset axial index */

	  i = 0;
	  
	  /* Loop over axial zones to update temperatures and hot radii */

	  while(axi > VALID_PTR)
	    {

	      if(RDB[axi + IFC_FUEP_AX_N_ANG] > 1)
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

	      if(r_cold[i][0] != 0.0)
		{
		  /* Radius is zero */

		  WDB[crptr + 0 ] = 0.0;
		  WDB[hrptr + 0 ] = 0.0;

		  /* Calculate density factor of central hole       */
		  /* (Should be dependent on the other gas volumes) */

		  f = (r_cold[i][0]*r_cold[i][0])/(r_hot[i][0]*r_hot[i][0]);

		  /* Cutoff for larger than unity density factors */

		  if (f > 1.0)
		    f = 1.0;

		  /* Store density factor */

		  WDB[dfptr + 0 ] =  f;

		  /* Calculate local convergence criterion */
		  if(RDB[t0ptr + 0] == 0.0)
		    Die(FUNCTION_NAME, "Zero temperature");

		  /* Calculate pointwise relative convergence criterion */

		  eps = fabs((RDB[t1ptr + 0]-T[i][0])/RDB[t1ptr + 0]);
		  dT = fabs(RDB[t1ptr + 0]-T[i][0]);

		  /* Add to absolute L2 norm of difference */

		  L2abs = L2abs + (RDB[t1ptr + 0]-T[i][0])*
		    (RDB[t1ptr + 0]-T[i][0]);

		  /* Add to relative L2 norm of difference */

		  L2rel = L2rel + (RDB[t1ptr + 0]-T[i][0])*
		    (RDB[t1ptr + 0]-T[i][0])/
		    (T[i][0]*T[i][0]);

		  /* Compare to maximum and store values */
		  if(eps > maxeps)
		    {
		      maxeps = eps;
		      maxr = r_cold[i][0];
		      maxT = T[i][0];
		      maxTp = RDB[t1ptr + 0];
		    }

		  /* Store new maximum absolute diference */

		  if (dT > maxdT)
		    maxdT = dT;

		  /* Update centerline temperatures */

		  WDB[t0ptr  + 0 ] = T0[i][0];
		  WDB[t1ptr  + 0 ] = T[i][0];

		  /* Increment radial node count */

		  k = 1;
		}

	      /* Loop over remaining radial nodes */

	      for(j = 0; j < nnodes[1] + nnodes[2]; j++)
		{
		  /* Write new cold radius */

		  WDB[crptr + k ] = (r_cold[i][j]*100)*(r_cold[i][j]*100);

		  /* Write new hot radius */

		  WDB[hrptr + k ] = (r_hot[i][j]*100)*(r_hot[i][j]*100);

		  /* Calculate density factor */

		  if(j == 0)
		    {
		      /* Innermost radial node */

		      /* Calculate density factor */

		      if(r_hot[i][j] > 0)
			f = (r_cold[i][j]*r_cold[i][j])/(r_hot[i][j]*r_hot[i][j]);
		      else
			f = 1.0;
		    }
		  else
		    {
		      /* Not innermost radial node*/

		      
		      /* Calculate density factor */

		      f = (r_cold[i][j]*r_cold[i][j] 
			   - r_cold[i][j-1]*r_cold[i][j-1])/
			(r_hot[i][j]*r_hot[i][j] 
			 - r_hot[i][j-1]*r_hot[i][j-1]);
		    }

		  /* Cutoff for larger than unity density factors */

		  if (f > 1.0)
		    f = 1.0;

		  /* Store density factor */

		  WDB[dfptr + k ] = f;

		  /* Calculate local convergence criterion */

		  eps = fabs((RDB[t1ptr + k]-T[i][j])/RDB[t1ptr + k]);
		  dT =  fabs(RDB[t1ptr + k]-T[i][j]);

		  /* Add to absolute L2 norm of difference */

		  L2abs = L2abs + (RDB[t1ptr + k]-T[i][k])*
		    (RDB[t1ptr + k]-T[i][k]);

		  /* Add to relative L2 norm of difference */

		  L2rel = L2rel + (RDB[t1ptr + k]-T[i][k])*
		    (RDB[t1ptr + k]-T[i][k])/
		    (T[i][k]*T[i][k]);

		  /* Compare to maximum and store values*/

		  if(eps > maxeps)
		    {
		      maxeps = eps;
		      maxr = r_cold[i][j];
		      maxT = T[i][j];
		      maxTp = RDB[t1ptr + k];
		    }

		  /* Store new maximum absolute diference */

		  if (dT > maxdT)
		    maxdT = dT;

		  /* Update temperatures */

		  WDB[t0ptr  + k ] = T0[i][j];
		  WDB[t1ptr  + k ] = T[i][j];
	      
		  /* Increment radial node count */
		  k++;
		}

	      /* Next axial segment */

	      axi = NextItem(axi);
	      i++;
	    }

	}
      else
	{
	  /* Loop over time steps */

	  tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

	  while(tbi > VALID_PTR)
	    {

	      /* Don't update if not in the right time interval */
	      if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
		{

		  if(!((RDB[tbi + IFC_FUEP_T_TMIN] >= tming-1E-15) && (RDB[tbi + IFC_FUEP_T_TMAX] <= tmaxg+1E-15)))
		    {
		      tbi = NextItem(tbi);
		      continue;
		    }

		}

	      if((prev = PrevItem(tbi)) > VALID_PTR)
		{
		  /* Not first timestep, get temperatures from prev step*/

		  T0 = (double**)((long)RDB[prev + IFC_FUEP_FINIX_PTR_T]);

		}
	      else
		{
		  /* First timestep, get steady state temperature */
		  fib = (long)RDB[DATA_PTR_FIN0];
		  CheckPointer(FUNCTION_NAME, "(fib)", DATA_ARRAY, fib);

		  /* Find correct pin */

		  while(fib > VALID_PTR)
		    {
		      /* Get Finix pin uni */

		      found=0;

		      /* Pointer to number of rod segments */

		      nu = WDB[fpe + IFC_FUEP_N_UNI];

		      /* Loop over rod segments */

		      ptr = (long)RDB[fpe + IFC_FUEP_PTR_UNI_LIST];

		      for(i=0; i < nu; i++)
			{
			  if(CompareStr(ptr + i, fib + FINIX_PTR_UNI_NAME))
			    found=1;
			}

		      if(found==1)
			break;
		      else
			fib = NextItem(fib);	 

		    }

		  /* Check if found */
		  if(fib < VALID_PTR)
		    Die(FUNCTION_NAME, "Could not find universe");

		  /* Get steady state temperature distribution */

		  T0 = (double**)((long)RDB[fib + FINIX_PTR_T]);

		}

	      /* Get Finix pointers */
	  
	      nnodes = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	      r_hot = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R]);
	      r_cold = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD]);
	      T = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_T]);

	      /* Get axial bin pointer */
	      axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];
	      CheckPointer(FUNCTION_NAME, "(axi)", DATA_ARRAY, axi);

	      /* Reset axial index */

	      i = 0;
	  
	      /* Loop over axial zones to update temperatures and hot radii */

	      while(axi > VALID_PTR)
		{

		  if(RDB[axi + IFC_FUEP_AX_N_ANG] > 1)
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

		  if(r_cold[i][0] != 0.0)
		    {
		      /* Radius is zero */

		      WDB[crptr + 0 ] = 0.0;
		      WDB[hrptr + 0 ] = 0.0;

		      /* Calculate density factor of central hole       */
		      /* (Should be dependent on the other gas volumes) */

		      f = (r_cold[i][0]*r_cold[i][0])/(r_hot[i][0]*r_hot[i][0]);

		      /* Cutoff for larger than unity density factors */

		      if (f > 1.0)
			f = 1.0;

		      /* Store density factor */

		      WDB[dfptr + 0 ] =  f;

		      /* Update redo-flag */

		      if(fabs(WDB[t1ptr  + 0 ] - T[i][0]) > 3.0)
			redo = 1;

		      /* Calculate local convergence criterion */
		      if(RDB[t0ptr + 0] == 0.0)
			Die(FUNCTION_NAME, "Zero temperature");

		      /* Calculate pointwise relative convergence criterion */

		      eps = fabs((RDB[t1ptr + 0]-T[i][0])/RDB[t1ptr + 0]);
		      dT  = fabs(RDB[t1ptr + 0]-T[i][0]);

		      /* Add to absolute L2 norm of difference */

		      L2abs = L2abs + (RDB[t1ptr + 0]-T[i][0])*
			(RDB[t1ptr + 0]-T[i][0]);

		      /* Add to relative L2 norm of difference */

		      L2rel = L2rel + (RDB[t1ptr + 0]-T[i][0])*
			(RDB[t1ptr + 0]-T[i][0])/
			(T[i][0]*T[i][0]);

		      /* Compare to maximum and store values */
		      if(eps > maxeps)
			{
			  maxeps = eps;
			  maxr = r_cold[i][0];
			  maxT = T[i][0];
			  maxTp = RDB[t1ptr + 0];
			}

		      if (dT > maxdT)
			maxdT = dT;

		      /* Update centerline temperatures */

		      WDB[t0ptr  + 0 ] = T0[i][0];
		      WDB[t1ptr  + 0 ] = T[i][0];

		      /* Increment radial node count */

		      k = 1;
		    }

		  /* Loop over remaining radial nodes */

		  for(j = 0; j < nnodes[1] + nnodes[2]; j++)
		    {
		      /* Write new cold radius */

		      WDB[crptr + k ] = (r_cold[i][j]*100)*(r_cold[i][j]*100);

		      /* Write new hot radius */

		      WDB[hrptr + k ] = (r_hot[i][j]*100)*(r_hot[i][j]*100);

		      /* Calculate density factor */

		      if(j == 0)
			{
			  /* Innermost radial node */

			  /* Calculate density factor */

			  if(r_hot[i][j] > 0)
			    f = (r_cold[i][j]*r_cold[i][j])/(r_hot[i][j]*r_hot[i][j]);
			  else
			    f = 1.0;
			}
		      else
			{
			  /* Not innermost radial node*/

		      
			  /* Calculate density factor */

			  f = (r_cold[i][j]*r_cold[i][j] 
			       - r_cold[i][j-1]*r_cold[i][j-1])/
			    (r_hot[i][j]*r_hot[i][j] 
			     - r_hot[i][j-1]*r_hot[i][j-1]);
			}

		      /* Cutoff for larger than unity density factors */

		      if (f > 1.0)
			f = 1.0;

		      /* Store density factor */

		      WDB[dfptr + k ] = f;

		      /* Calculate local convergence criterion */

		      eps = fabs((RDB[t1ptr + k]-T[i][j])/RDB[t1ptr + k]);

		      dT = fabs(RDB[t1ptr + k]-T[i][j]);

		      /* Add to absolute L2 norm of difference */

		      L2abs = L2abs + (RDB[t1ptr + k]-T[i][k])*
			(RDB[t1ptr + k]-T[i][k]);

		      /* Add to relative L2 norm of difference */

		      L2rel = L2rel + (RDB[t1ptr + k]-T[i][k])*
			(RDB[t1ptr + k]-T[i][k])/
			(T[i][k]*T[i][k]);

		      /* Compare to maximum and store values*/

		      if(eps > maxeps)
			{
			  maxeps = eps;
			  maxr = r_cold[i][j];
			  maxT = T[i][j];
			  maxTp = RDB[t1ptr + k];
			}

		      if (dT > maxdT)
			maxdT = dT;

		      /* Update temperatures */

		      WDB[t0ptr  + k ] = T0[i][j];
		      WDB[t1ptr  + k ] = T[i][j];
	      
		      /* Increment radial node count */
		      k++;
		    }

		  /* Next axial segment */

		  axi = NextItem(axi);
		  i++;
		}

	      /* Also handle upcoming time interval (if it exists) */

	      /* Next time bin */

	      if((tbi = NextItem(tbi)) < VALID_PTR)
		break;

	      T0 = T;

	      r_hot = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R]);
	      r_cold = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD]);
	      T = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_T]);	      

	      /* Get axial bin pointer */
	      axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];
	      CheckPointer(FUNCTION_NAME, "(axi)", DATA_ARRAY, axi);

	      /* Reset axial index */

	      i = 0;
	  
	      /* Loop over axial zones to update temperatures and hot radii */

	      while(axi > VALID_PTR)
		{

		  if(RDB[axi + IFC_FUEP_AX_N_ANG] > 1)
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

		  if(r_cold[i][0] != 0.0)
		    {
		      /* Radius is zero */

		      WDB[crptr + 0 ] = 0.0;
		      WDB[hrptr + 0 ] = 0.0;

		      /* Calculate density factor of central hole       */
		      /* (Should be dependent on the other gas volumes) */

		      f = (r_cold[i][0]*r_cold[i][0])/(r_hot[i][0]*r_hot[i][0]);

		      /* Cutoff for larger than unity density factors */

		      if (f > 1.0)
			f = 1.0;

		      /* Store density factor */

		      WDB[dfptr + 0 ] =  f;

		      /* Update centerline temperatures */
		      /* No time variation in temperature */

		      WDB[t0ptr  + 0 ] = T0[i][0];
		      WDB[t1ptr  + 0 ] = T[i][0];

		      /* Increment radial node count */

		      k = 1;
		    }

		  /* Loop over remaining radial nodes */

		  for(j = 0; j < nnodes[1] + nnodes[2]; j++)
		    {
		      /* Write new cold radius */

		      WDB[crptr + k ] = (r_cold[i][j]*100)*(r_cold[i][j]*100);

		      /* Write new hot radius */

		      WDB[hrptr + k ] = (r_hot[i][j]*100)*(r_hot[i][j]*100);

		      /* Calculate density factor */

		      if(j == 0)
			{
			  /* Innermost radial node */

			  /* Calculate density factor */

			  if(r_hot[i][j] > 0)
			    f = (r_cold[i][j]*r_cold[i][j])/(r_hot[i][j]*r_hot[i][j]);
			  else
			    f = 1.0;
			}
		      else
			{
			  /* Not innermost radial node*/

		      
			  /* Calculate density factor */

			  f = (r_cold[i][j]*r_cold[i][j] 
			       - r_cold[i][j-1]*r_cold[i][j-1])/
			    (r_hot[i][j]*r_hot[i][j] 
			     - r_hot[i][j-1]*r_hot[i][j-1]);
			}

		      /* Cutoff for larger than unity density factors */

		      if (f > 1.0)
			f = 1.0;

		      /* Store density factor */

		      WDB[dfptr + k ] = f;

		      /* Update temperatures */

		      WDB[t0ptr  + k ] = T0[i][j];
		      WDB[t1ptr  + k ] = T[i][j];
	      
		      /* Increment radial node count */
		      k++;
		    }

		  /* Next axial segment */

		  axi = NextItem(axi);
		  i++;
		}	

	      /* Do not handle more time intervals */

	      tbi = -1;
	    }

	}
      /* Next rod */

      fpe = NextItem(fpe);
    }

  L2abs = sqrt(L2abs);
  L2rel = sqrt(L2rel);
  
  fprintf(out, "Max Eps was %E (%f -> %f) at %E\n", maxeps, maxTp, maxT, maxr);
  fprintf(out, "L2 abs %E rel %E\n", L2abs, L2rel);
  fprintf(out, "%E %%Epsmax\n", maxeps);
  fprintf(out, "%E %%dTmax\n", maxdT);
  
  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      if (maxdT > 3.0)
	WDB[DATA_ITERATE] = (double)YES;
      
    }

  /*printf("Done...\n");*/
}	  

#endif

/*****************************************************************************/
