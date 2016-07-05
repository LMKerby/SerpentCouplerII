/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : freefinix.c                                    */
/*                                                                           */
/* Created:       2014/11/05 (VVa)                                           */
/* Last modified: 2014/11/05 (VVa)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Frees all finix arrays                                       */
/*                                                                           */
/* Comments:  -                                                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "FreeFinix:"

/*****************************************************************************/

void FreeFinix()
{
  long fib, fpe, ifc, tbi;
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

  /* Get pointer to first finix pin */

  fib = (long)RDB[DATA_PTR_FIN0];

  /* Get pointer to interface */

  ifc = (long)RDB[fib + FINIX_PTR_IFC];

  /* Pointer to FUEP BLOCK */

  fpe = (long)RDB[ifc + IFC_PTR_FUEP];

  /* Loop over pins */
  while(fpe > VALID_PTR)
    {

      /* Loop over time steps */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];

      while(tbi > VALID_PTR)
	{

	  /* Get pointers to arrays */

	  nnodes   = (int*)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  T        = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_T]);
	  r        = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_R]);
	  r_cold   = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD]);
	  Pden     = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_PDEN]);
	  Lhr      = (double*)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_LHR]);
	  bu       = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_BU]);
	  params   = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_PARAMS]);
	  bcond    = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_BCOND]);
	  sresults = (double*)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_SRESULTS]);
	  vresults = (double**)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_VRESULTS]);
	  options  = (int*)((long)WDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);

	  /* Free arrays */

	  finix_free_arrays(&nnodes, &T, &r, &r_cold, &Pden, 
			    &Lhr, &bu, &params, &bcond, 
			    &sresults, &vresults,&options);

	  /* Next time bin */

	  tbi = NextItem(tbi);
	}

      /* Next rod */

      fpe = NextItem(fpe);
    }

  /* Loop over finix pins */

  while(fib > VALID_PTR)
    {
      /* Get Finix pointers */
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

      /* Free arrays */

      finix_free_arrays(&nnodes, &T, &r, &r_cold, &Pden, 
			&Lhr, &bu, &params, &bcond, 
			&sresults, &vresults,&options);

      /* Next finix pin*/

      fib = NextItem(fib);
    }

}	  

#endif

/*****************************************************************************/
