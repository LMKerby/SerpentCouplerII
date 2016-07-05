/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : updatefinixpower.c                             */
/*                                                                           */
/* Created:       2013/03/27 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Updates the power density for Finix nodes                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "UpdateFinixPower:"

/*****************************************************************************/
/* Calculates power density at axial region zi, at radius r*/
double GetPden(double r, long buf, long zi, double h, long nr, double rf0, double rf1, long ti, long ptr)
{
  long ri=0;
  long nz, na, nt, idx;
  double val;

  /* Get correct radial bin index */

  while(r*r >= ri/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0)
    {
      ri++;
    }

  ri--;

  if(ri == nr)
    ri--;

  /* Get number of axial, angular and time bins */

  nz = (long)RDB[ptr + FUEP_NZ];	  
  na = (long)RDB[ptr + FUEP_NA];      
  nt = (long)RDB[ptr + FUEP_NT];

  /* Calculate global index */

  idx = zi + nz*0 + nz*1*ri + nz*1*nr*ti;

  /* Convert power to power density */
  
  val = RDB[buf + idx]/(h*PI*((rf1*rf1-rf0*rf0)/(double)nr));

  /* Return power density */

  return val;
}

void UpdateFinixPower(long fib, long loc0)
{
  long ptr, i,j,k, n, m, limptr;
  long nr, nz, na, nt, pbuf, pbuf2, tb, tbi;
  double r0, r1, sum;
  double rb, rf0, rf1, A, B, C, r_;
  double zmin, zmax, h, tming, tmaxg;
  int *nnodes;
  double **r;
  double **r_cold;
  double **Pden;
  double *Lhr;

  /* Check if found */

  if(loc0 < VALID_PTR)
    Die(FUNCTION_NAME, "Could not find universe");

  /* Pointer to power buffer */

  /* These are private to MPI tasks (haven't been broadcast anywhere) */
  /* Should be identical though */

  /* Created in AllocInterfaceStat() */
  pbuf = (long)RDB[loc0 + IFC_FUEP_PTR_POWER];
  CheckPointer(FUNCTION_NAME, "(ptr buf)", DATA_ARRAY, pbuf);

  /* Created in AllocInterfaceStat() */
  pbuf2 = (long)RDB[loc0 + IFC_FUEP_PTR_POWER_REL];
  CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, pbuf2);

  /* Created in ReadIFCFB */
  ptr = (long)RDB[loc0 + IFC_FUEP_OUT_PTR_LIM];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  limptr = ptr;

  /* Get number of bins in different dimensions */

  nr = (long)RDB[ptr + FUEP_NR];
  nz = (long)RDB[ptr + FUEP_NZ];
  na = (long)RDB[ptr + FUEP_NA];
  nt = (long)RDB[ptr + FUEP_NT];

  /* Reset time bin number */

  tb = 0;

  /* Get transport time bin limits */

  tmaxg = RDB[DATA_TIME_CUT_TMAX];

  tming = RDB[DATA_TIME_CUT_TMIN];

  /* Loop over time bins */

  tbi = (long)RDB[loc0 + IFC_FUEP_PTR_T];

  while(tbi > VALID_PTR)
    {

      /* Don't update if not in the right time interval */
      if ((RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT) && 
	  (!((RDB[tbi + IFC_FUEP_T_TMIN] >= tming-1E-15) && 
	     (RDB[tbi + IFC_FUEP_T_TMAX] <= tmaxg+1E-15))))
	{

	  /* Check that not in infinite time bin (steady state) */

	  if(RDB[tbi + IFC_FUEP_T_TMAX] < INFTY/100)
	    {

	      tb++;
	      tbi = NextItem(tbi);
	      continue;
	    }
	}

      /* Get the Finix pointers for this time bin */

      if (RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
	{
  
	  nnodes = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  Pden = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_PDEN]);
	  Lhr = (double*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_LHR]);
	  r = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R]);
	  r_cold = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD]);

	}
      else
	{

	  nnodes = (int*)((long)RDB[fib + FINIX_PTR_NNODES]);
	  Pden = (double**)((long)RDB[fib + FINIX_PTR_PDEN]);
	  Lhr = (double*)((long)RDB[fib + FINIX_PTR_LHR]);
	  r = (double**)((long)RDB[fib + FINIX_PTR_R]);
	  r_cold = (double**)((long)RDB[fib + FINIX_PTR_R_COLD]);

	}

      /* Loop over axial regions */
      for(i=0; i< nnodes[0]; i++)
	{
	  /* Get rmin and rmax (not squared)*/

	  rf0 = RDB[ptr + FUEP_RMIN];
	  rf1 = RDB[ptr + FUEP_RMAX];

	  /* Check that the r-limits of power estimators match those of */
	  /* Finix nodes */

	  /* Fuel inner radius */

	  if(rf0 - 100*r_cold[i][0] > ZERO)
	    Die(FUNCTION_NAME, 
		"Mismatch in fuel inner radius (%E) (%E) (diff. %E)\n"
		,rf0,100*r_cold[i][0],rf0-100*r_cold[i][0]);

	  /* Fuel outer radius */

	  if(rf1 - 100*r_cold[i][nnodes[1]-1] > ZERO)
	    Die(FUNCTION_NAME, 
		"Mismatch in fuel outer radius (%E) (%E) (diff. %E)\n"
		,rf1,100*r_cold[i][nnodes[1]-1],rf1-100*r_cold[i][nnodes[1]-1]);

	  /* Get number of radial power nodes */

	  nr = (long)RDB[ptr + FUEP_NR];

	  /* This is the radial index of the current finix node */
	  n = 0;

	  /* This is the radial index of the current power tally */
	  j = 0;

	  /* Get the z-limits of the current tally */

	  zmin = RDB[ptr + FUEP_ZMIN];
	  zmax = RDB[ptr + FUEP_ZMAX];

	  /* Check 2D calculation */

	  if ((zmin < INFTY / 100.0) || (zmax > INFTY / 100.0))
	    h = 1.0;
	  else
	    h = (zmax-zmin)/RDB[ptr + FUEP_NZ];

	  /* Calculate total power in this axial zone for LHR */

	  sum=0.0;

	  for(k=0;k < nr; k++)
	    {
	      sum += RDB[pbuf2 + i + 0*nz + k*1*nz + tb*1*nz*nr];
	    } 
	  if(sum==0)
	    sum = 1E-4;

	  /* Calculate and store linear heat rate */

	  /* Relaxation is handled on the power tally in */
	  /* relaxinterfacepower.c */

	  Lhr[i] = sum/h*100;

	  fprintf(out, "LHR = %E\n", sum/h);

	  while(n<nnodes[1]){

	    /* First finix node is a special case */

	    if(n==0){

	      /* Get inner and outer radius */

	      r0=r_cold[i][n]*100;
	      r1=r_cold[i][n+1]*100;

	      /* Get power densities at inner and outer radii */

	      A=GetPden(r0, pbuf2, i, h, nr, rf0, rf1, tb, limptr);
	      B=GetPden(r1, pbuf2, i, h, nr, rf0, rf1, tb, limptr);

	      /* This is the pden in the middle of the nodes (node boundary) */
	      C=GetPden((r0+r1)/2.0, pbuf2, i, h, nr, rf0, rf1, tb, limptr);

	
	      /* If there is a power tally boundary in the first "control volume"*/
	      if(A!=C){
		/* Find the power tally boundary */
		m=0;
		while(r0 >= (rb = sqrt(m/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0))){
		  m++;
		}

		/* If there are two or more boundaries between these nodes things get hard */
		if(C != B)
		  Die(FUNCTION_NAME, "Please increase the number of Finix radial nodes");
	  
		/* Calculate power density by a very special average */
		/* To conserve the power both in axial slice and in this node volume */
		C=(A*(rb*rb-r0*r0) + 
		   B*(-2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)+r1*r1-rb*rb))
		  /(-2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)+r1*r1-r0*r0);  
	      }

	    }else if(n==nnodes[1]-1){

	      /* The last finix node is also a special case */
	      r0=r_cold[i][n]*100;
	      r_=r_cold[i][n-1]*100;

	      /* Get power densities at inner and outer radii */

	      A=GetPden(r_, pbuf2, i, h, nr, rf0, rf1, tb, limptr);
	      B=GetPden(r0, pbuf2, i, h, nr, rf0, rf1, tb, limptr);

	      /* This is the pden in the middle of the nodes */
	      /* +1E-6 because the intervals are (r_,r1]*/
	      C=GetPden((r_+r0)/2.0+1E-6, pbuf2, i, h, nr, rf0, rf1, tb, limptr);

	      /* If there is a power tally boundary in the last "control volume"*/
	      if(B!=C){
		/* Find the power tally boundary */
		m=0;
		while(r_ >=  (rb = sqrt(m/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0))){
		  m++;
		}

		/* If there are two or more boundaries between these nodes things get hard */
		if(C != A)
		  Die(FUNCTION_NAME, "Please increase the number of Finix radial nodes");

		/* Calculate power density by a very special average */
		/* To conserve the power both in axial slice and in this node volume */
		C=(A*(2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_)+rb*rb-r0*r0)+B*(r0*r0-rb*rb))/
		  (2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_));  
	      }

	    }else{

	      /* Get radii of adjacent nodes */

	      r_=r_cold[i][n-1]*100;
	      r0=r_cold[i][n]*100;
	      r1=r_cold[i][n+1]*100;

	      /* Get power densities at inner and outer radii of node volume */
	      /* and at node point */	  
	      A=GetPden((r_+r0)*0.5+1E-6, pbuf2, i, h, nr, rf0, rf1, tb, limptr);
	      B=GetPden((r0+r1)*0.5, pbuf2, i, h, nr, rf0, rf1, tb, limptr);
	      C=GetPden(r0, pbuf2, i, h, nr, rf0, rf1, tb, limptr);
	  
	      /* If there is a power tally boundary in the "control volume" around this node */

	      if(A!=B){
   
		/* Find boundary*/
		m=0;
		while((r_+r0)*0.5+1E-6 >=  (rb = sqrt(m/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0))){
		  m++;
		}
	    
		/* If there are two or more boundaries between these nodes things get hard */
		if((r0+r1)*0.5 > sqrt((m+1)/(double)nr*(rf1*rf1-rf0*rf0)+rf0*rf0))
		  Die(FUNCTION_NAME, "Please increase the number of Finix radial nodes");	      

		/* Get power densities at the adjacent node points */
	       
		A=GetPden(r_, pbuf2, i, h, nr, rf0, rf1, tb, limptr);
		B=GetPden(r1, pbuf2, i, h, nr, rf0, rf1, tb, limptr);	  

		/* Calculate power density by a very special average */
		/* To conserve the power both in axial slice and in this node volume */
	    
		C=(A*(2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_)+rb*rb-r0*r0)+
		   B*(-2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)+r1*r1-rb*rb))/
		  (2.0/3.0*(r0*r0+r0*r_+r_*r_)-r_*(r0+r_)
		   -2.0/3.0*(r1*r1+r1*r0+r0*r0)+r0*(r1+r0)
		   +r1*r1-r0*r0);
	      }

	    }

	    /* Store new power density */

	    /* The thing stored is not precisely power density but some kind of */
	    /* power factor (see FINIX source for more information)             */

	    Pden[i][n]=(C*1E6)/Lhr[i]*PI*r[i][nnodes[1]]*r[i][nnodes[1]];

	    /* Move to the next node */
	    n++;
	  }


	}

      /*printf("\n");*/

      /* Next time bin */

      tb++;
      tbi = NextItem(tbi);
    }


}

#endif

/*****************************************************************************/
