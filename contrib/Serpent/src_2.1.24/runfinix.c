/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : runfinix.c                                     */
/*                                                                           */
/* Created:       2013/03/27 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Runs Finix for steady state or transient                     */
/*                                                                           */
/* Comments:  - Transientti pitää vielä kirjoittaa                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "RunFinix:"

/*****************************************************************************/

void RunFinix(long fib, long fpe)
{
  long tbi, tbi2, tb,i,j, prev;
  double tmin, tmax, tming, tmaxg;
  int *nnodes;
  double **T;
  double **Tprev;
  double **r;
  double **rprev;
  double **r_cold;
  double **Pden;
  double *Lhr;
  double **bu;
  double **params;
  double **bcond;
  double *sresults;
  double **vresults;
  int *options;
  char **err=NULL;

  /* Check steady state */

  if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
    {

      /* Criticality source simulation */

      /* Get Finix pointers for steady state */

      nnodes   = (int*)((long)RDB[fib + FINIX_PTR_NNODES]);
      T        = (double**)((long)RDB[fib + FINIX_PTR_T]);
      r        = (double**)((long)RDB[fib + FINIX_PTR_R]);
      r_cold   = (double**)((long)RDB[fib + FINIX_PTR_R_COLD]);
      Pden     = (double**)((long)RDB[fib + FINIX_PTR_PDEN]);
      Lhr      = (double*)((long)RDB[fib + FINIX_PTR_LHR]);
      bu       = (double**)((long)RDB[fib + FINIX_PTR_BU]);
      params   = (double**)((long)RDB[fib + FINIX_PTR_PARAMS]);
      bcond    = (double**)((long)RDB[fib + FINIX_PTR_BCOND]);
      sresults  = (double*)((long)RDB[fib + FINIX_PTR_SRESULTS]);
      vresults  = (double**)((long)RDB[fib + FINIX_PTR_VRESULTS]);
      options  = (int*)((long)RDB[fib + FINIX_PTR_OPTIONS]);

      /* Solve steady state solution */

      fprintf(out, "Solving FINIX steady state for rod %s\n", 
	      GetText(fib + FINIX_PTR_UNI_NAME));

      err = finix_solve_initial_steady_state(nnodes,T,r,r_cold,Pden,Lhr,bu,params,bcond,sresults,vresults,options);

      /* Print possible errors */

      if(err!=NULL)
	finix_printf_err(err);
	  

      finix_free_err(&err);
      /*
      finix_printf_array(nnodes,T);
	      
      finix_printf_array(nnodes,Pden);
      */
    }
  else
    {
      /* Time dependent simulation */

      /* Get transport time bin limits */

      tmaxg = RDB[DATA_TIME_CUT_TMAX];

      tming = RDB[DATA_TIME_CUT_TMIN];

      /* Loop over time bins */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      tb = 0;
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
 
	  /* Get Finix pointers for the time step*/
	  nnodes   = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  T        = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_T]);
	  r        = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R]);
	  r_cold   = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD]);
	  Pden     = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_PDEN]);
	  Lhr      = (double*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_LHR]);
	  bu       = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BU]);
	  params   = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_PARAMS]);
	  bcond    = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BCOND]);
	  sresults  = (double*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_SRESULTS]);
	  vresults  = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_VRESULTS]);
	  options  = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);

	  /* Get temperature at the beginning of time step if available */

	  if((prev = PrevItem(tbi)) > VALID_PTR)
	    {
	      Tprev = (double**)((long)RDB[prev + IFC_FUEP_FINIX_PTR_T]);
	    }
	  else
	    {
	      /* This is the first step, get initial steady state temperature */
	      /* via steady state calculation*/

	      Tprev = (double**)((long)RDB[fib + FINIX_PTR_T]);

	    }

	  /* Run Finix for the current step */

	  tmin = RDB[tbi + IFC_FUEP_T_TMIN];
	  tmax = RDB[tbi + IFC_FUEP_T_TMAX];

	  fprintf(out,"Solving finix transient from %E to %E (%E s)\n", tmin, tmax, tmax-tmin);

	  /* Solve transient */

	  err = finix_solve_transient(tmax-tmin, nnodes,Tprev,T,r,r_cold,Pden,Lhr,bu,params,bcond,sresults,vresults,options);

	  /* Print errors */

	  if(err!=NULL)
	    finix_printf_err(err);	  

	  finix_free_err(&err);

	  fprintf(out, "Temperatures:\n");
	  finix_printf_array(nnodes,T);	    
	      
	  fprintf(out, "\nPowers:\n");
	  finix_printf_array(nnodes,Pden);
	  fprintf(out, "\n");

	  /* Solve the temperature distribution after next interval also */

	  tbi2 = NextItem(tbi);

	  if (tbi2 < VALID_PTR)
	    {
	      /* Next time bin */

	      tbi = NextItem(tbi);
	      tb++;

	      continue;
	    }

	  Tprev = T;

	  /* Get pointers to temperature and radius on the future */
	  /* interval                                             */
		  
	  T = (double**)((long)RDB[tbi2 + IFC_FUEP_FINIX_PTR_T]);
	  r = (double**)((long)RDB[tbi2 + IFC_FUEP_FINIX_PTR_R]);

	  tmin = RDB[tbi2 + IFC_FUEP_T_TMIN];
	  tmax = RDB[tbi2 + IFC_FUEP_T_TMAX];

	  fprintf(out,"Solving finix transient from %E to %E (%E s)\n", tmin, tmax, tmax-tmin);

	  /* Solve transient */

	  err = finix_solve_transient(tmax-tmin, nnodes,Tprev,T,r,r_cold,Pden,Lhr,bu,params,bcond,sresults,vresults,options);

	  /* Next time bin */

	  tbi = NextItem(tbi);
	  tb++;
	}

    }
}

#endif

/*****************************************************************************/
