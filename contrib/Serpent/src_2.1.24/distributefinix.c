/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : distributefinix.c                              */
/*                                                                           */
/* Created:       2013/10/11 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Sets up FINIX time bin data to interface structures          */
/*              this can not be done in processfinix as the interface        */
/*              structures are not present then                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "DistributeFinix:"

/*****************************************************************************/

void DistributeFinix()
{
  long i, j, fib, ifc;
  long fpe, found, nu, ptr, tbi;
  int *nnodes;
  double **T;
  double **r;
  double **r_cold;
  double **Pden;
  double *Lhr;
  double **bu;
  double **params;
  double **bcond;
  double *sresults;
  double **vresults;
  int *options;
  int *nnodes0;
  double **r0;
  double **r_cold0;
  double **params0;
  double **bcond0;
  int *options0;
  char **err=NULL, **allerr=NULL;

  if ((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  fprintf(out,"Distributing FINIX data to timesteps\n");

  while(fib > VALID_PTR)
    {
      /* Get Finix pointers */
      nnodes0   = (int*)((long)RDB[fib + FINIX_PTR_NNODES]);
      r0        = (double**)((long)RDB[fib + FINIX_PTR_R]);
      r_cold0   = (double**)((long)RDB[fib + FINIX_PTR_R_COLD]);
      params0   = (double**)((long)RDB[fib + FINIX_PTR_PARAMS]);
      bcond0    = (double**)((long)RDB[fib + FINIX_PTR_BCOND]);
      options0  = (int*)((long)RDB[fib + FINIX_PTR_OPTIONS]);

      /* Get pointer to corresponding interface */

      ifc = (long)RDB[fib + FINIX_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

      /* Find correct pin from interface */
      /* Pointer to FUEP BLOCK */
      fpe = (long)RDB[ifc + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(fpe)", DATA_ARRAY, fpe);

      while(fpe > VALID_PTR)
	{
	  /* Get Finix pin uni */

	  found=0;

	  /* Get number of universes inthis pin */

	  nu = WDB[fpe + IFC_FUEP_N_UNI];

	  /* Get pointer to rod segments */

	  ptr = (long)RDB[fpe + IFC_FUEP_PTR_UNI_LIST];

	  /* Compare to the universes of the interface FPE block */

	  for(i=0; i < nu; i++)
	    {
	      if(CompareStr(ptr + i, fib + FINIX_PTR_UNI_NAME))
		found=1;
	    }

	  /* Check if found */

	  if(found==1)
	    break;
	  else
	    fpe = NextItem(fpe);	 

	}

      /* Check if found */
      if(fpe < VALID_PTR)
	Die(FUNCTION_NAME, "Could not find universe");
      
      /* Get time bin pointer */
	  
      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

      while(tbi > VALID_PTR)
	{

	  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	    {
	      /* In criticality source mode, just copy steady state pointers here */
	      
	      WDB[tbi + IFC_FUEP_FINIX_PTR_NNODES] = RDB[fib + FINIX_PTR_NNODES];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_T] = RDB[fib + FINIX_PTR_T] ;
	      WDB[tbi + IFC_FUEP_FINIX_PTR_R] = RDB[fib + FINIX_PTR_R] ;
	      WDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD] = RDB[fib + FINIX_PTR_R_COLD];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_PDEN] = RDB[fib + FINIX_PTR_PDEN] ;
	      WDB[tbi + IFC_FUEP_FINIX_PTR_LHR] = RDB[fib + FINIX_PTR_LHR] ;
	      WDB[tbi + IFC_FUEP_FINIX_PTR_BU] = RDB[fib + FINIX_PTR_BU] ;
	      WDB[tbi + IFC_FUEP_FINIX_PTR_PARAMS] = RDB[fib + FINIX_PTR_PARAMS];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_BCOND] = RDB[fib + FINIX_PTR_BCOND] ;
	      WDB[tbi + IFC_FUEP_FINIX_PTR_SRESULTS] = RDB[fib + FINIX_PTR_SRESULTS];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_VRESULTS] = RDB[fib + FINIX_PTR_VRESULTS];
	      WDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS] = RDB[fib + FINIX_PTR_OPTIONS];

	      /* Next time bin (there should not be more than one) */

	      tbi = NextItem(tbi);

	      /* Cycle loop */

	      continue;
	    }

	  /* Initialize the NNODES array */
	  nnodes = (int*)finix_get_default_nnodes();

	  /* Set number of nodes*/
	  nnodes[0] = (int)nnodes0[0];
	  nnodes[1] = (int)nnodes0[1];
	  nnodes[2] = (int)nnodes0[2];

	  /* Initialize arrays for use */

	  err = finix_initialize_arrays(&nnodes,&T,&r,&r_cold,
					&Pden,&Lhr,&bu,&params,
					&bcond,&sresults, &vresults);
	  finix_append_err(&allerr, &err);

	  /* Get default options */

	  options = finix_get_default_options();
	  finix_append_err(&allerr, &err);

	  /* Get parameteres by rodype or from file */

	  if(RDB[fib + FINIX_RODTYPE] >= 0)
	    {
	      /* Get values from library */

	      err=finix_get_default_params(RDB[fib + FINIX_RODTYPE], params, sresults);
	      finix_append_err(&allerr, &err);

	    }
	  else
	    {
	      /* Get user defined values and replace them with user defined */

	      err = finix_read_params_file(GetText(fib + FINIX_PTR_FNAME), params, sresults);
	      finix_append_err(&allerr, &err);
	    }

	  /* Get default node positions (uniform radial intervals) */

	  err=finix_get_default_positions(&nnodes, r, r_cold, params);
	  finix_append_err(&allerr, &err);

	  /* Get default temperatures */

	  err=finix_get_default_cold_state(&nnodes, T);
	  finix_append_err(&allerr, &err);

	  /* Get default burnups (zero) */
	
	  err=finix_get_default_burnup(&nnodes, bu);
	  finix_append_err(&allerr, &err);

	  /* Get default power density (zero) */

	  err=finix_get_default_power(&nnodes, Pden, Lhr);
	  finix_append_err(&allerr, &err);

	  /* Get default boundary conds (rod outer surface T equal to coolant T */
	  /* that is given in params) */

	  err=finix_get_default_bcond(&nnodes, params, bcond);
	  finix_append_err(&allerr, &err);

	  /*calculate gas amount using cold state values*/

	  err=finix_gap_moles_from_pressure_cold(nnodes, T, r, r_cold, params, sresults, vresults);
	  finix_append_err(&allerr, &err);

	  /* Print collected errors at this point */

	  if(err!=NULL)
	    finix_printf_err(err);
    
	  finix_free_err(&err);

	  /* Set the node positions */

	  for(i=0;i<nnodes[0];i++){
	    for(j=0;j<nnodes[1]+nnodes[2];j++){
	      r[i][j]=r0[i][j];
	      r_cold[i][j]=r_cold0[i][j];
	    }
	  }

	  /* Set axial length */
	  params[0][4] = params0[0][4];
	  params[0][5] = params0[0][5];
	  params[0][6] = params0[0][6];

	  /* Set boundary conditions */
	  options[0] = options0[0];

	  for(i=0;i<4; i++)
	    for(j=0; j<nnodes[0];j++)
	      bcond[i][j] = bcond0[i][j];

	  /* Solve initial steady state (HZP) */
 
	  err = finix_solve_initial_steady_state(nnodes,T,r,r_cold,Pden,Lhr,bu,
						 params,bcond,sresults,vresults,
						 options);

	  /* Print out any errors */

	  if(err!=NULL)
	    finix_printf_err(err);

	  finix_free_err(&err);

	  /* Store pointers */
	  WDB[tbi + IFC_FUEP_FINIX_PTR_NNODES] = (long)nnodes;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_T] = (long)T;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_R] = (long)r;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD] = (long)r_cold;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_PDEN] = (long)Pden;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_LHR] = (long)Lhr;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_BU] = (long)bu;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_PARAMS] = (long)params;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_BCOND] = (long)bcond;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_SRESULTS] = (long)sresults;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_VRESULTS] = (long)vresults;
	  WDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS] = (long)options;
	 
	  /* Next time bin */
	  
	  tbi = NextItem(tbi);

	}

      /* Next fuel rod */

      fib = NextItem(fib);
    }

  fprintf(out, "OK.\n\n");

}

#endif

/*****************************************************************************/
