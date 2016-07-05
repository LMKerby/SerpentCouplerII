/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processfinix.c                                 */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Sets up FINIX data for temperature feedback                  */
/*              Runs all initialization FINIX routines and                   */
/*              prepares all arrays                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "ProcessFinix:"

/*****************************************************************************/

void ProcessFinix()
{
  long mat,  axi, loc0, i,j;
  long fib, nf, nc, nz, prev, npins;
  double rf0, rf1,  rc0, rc1;
  double TFBU[101];
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
  char **err=NULL, **allerr=NULL;

TFBU[1] = 9.4817;
TFBU[2] = 9.4817;
TFBU[3] = 9.4817;
TFBU[4] = 9.4817;
TFBU[5] = 9.4817;
TFBU[6] = 9.4817;
TFBU[7] = 9.4817;
TFBU[8] = 9.4817;
TFBU[9] = 9.4817;
TFBU[10] = 9.4817;
TFBU[11] = 9.4817;
TFBU[12] = 9.4817;
TFBU[13] = 9.4817;
TFBU[14] = 9.4817;
TFBU[15] = 9.4817;
TFBU[16] = 9.4817;
TFBU[17] = 9.4817;
TFBU[18] = 9.4817;
TFBU[19] = 9.4817;
TFBU[20] = 9.4817;
TFBU[21] = 9.4817;
TFBU[22] = 9.4817;
TFBU[23] = 9.5077;
TFBU[24] = 9.5077;
TFBU[25] = 9.5077;
TFBU[26] = 9.5077;
TFBU[27] = 9.5077;
TFBU[28] = 9.5077;
TFBU[29] = 9.5077;
TFBU[30] = 9.5077;
TFBU[31] = 9.5077;
TFBU[32] = 9.5077;
TFBU[33] = 9.5077;
TFBU[34] = 9.5077;
TFBU[35] = 9.5077;
TFBU[36] = 9.5077;
TFBU[37] = 9.5672;
TFBU[38] = 9.5672;
TFBU[39] = 9.5672;
TFBU[40] = 9.5672;
TFBU[41] = 9.5672;
TFBU[42] = 9.5672;
TFBU[43] = 9.5672;
TFBU[44] = 9.5672;
TFBU[45] = 9.5672;
TFBU[46] = 9.5672;
TFBU[47] = 9.5672;
TFBU[48] = 9.5672;
TFBU[49] = 9.6417;
TFBU[50] = 9.6417;
TFBU[51] = 9.6417;
TFBU[52] = 9.6417;
TFBU[53] = 9.6417;
TFBU[54] = 9.6417;
TFBU[55] = 9.6417;
TFBU[56] = 9.6417;
TFBU[57] = 9.6417;
TFBU[58] = 9.7252;
TFBU[59] = 9.7252;
TFBU[60] = 9.7252;
TFBU[61] = 9.7252;
TFBU[62] = 9.7252;
TFBU[63] = 9.7252;
TFBU[64] = 9.7252;
TFBU[65] = 9.7252;
TFBU[66] = 9.7252;
TFBU[67] = 9.8264;
TFBU[68] = 9.8264;
TFBU[69] = 9.8264;
TFBU[70] = 9.8264;
TFBU[71] = 9.8264;
TFBU[72] = 9.8264;
TFBU[73] = 9.8264;
TFBU[74] = 9.8264;
TFBU[75] = 9.9529;
TFBU[76] = 9.9529;
TFBU[77] = 9.9529;
TFBU[78] = 9.9529;
TFBU[79] = 9.9529;
TFBU[80] = 9.9529;
TFBU[81] = 9.9529;
TFBU[82] = 10.1278;
TFBU[83] = 10.1278;
TFBU[84] = 10.1278;
TFBU[85] = 10.1278;
TFBU[86] = 10.1278;
TFBU[87] = 10.1278;
TFBU[88] = 10.1278;
TFBU[89] = 10.4449;
TFBU[90] = 10.4449;
TFBU[91] = 10.4449;
TFBU[92] = 10.4449;
TFBU[93] = 10.4449;
TFBU[94] = 10.4449;
TFBU[95] = 11.7249;
TFBU[96] = 11.7249;
TFBU[97] = 11.7249;
TFBU[98] = 11.7249;
TFBU[99] = 11.7249;
TFBU[100] = 11.7249;
TFBU[101] = 11.7249;

  /* Check that some finix pins are defined */

  if ((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  /* Close list to allow omp-loops over it later */

  CloseList(fib);
  
  fprintf(out, "Processing Finix pins\n");

  /* Avoid compiler warning*/

  prev = 0;
  loc0 = 0;
  mat = 0;

  /* reset number of pins */

  npins = 0;

  /* Loop over feedbacks */
  while (fib > VALID_PTR)
    {

    /* Get fuel in and out radii */

    rf0 = RDB[fib + FINIX_RMIN];

    rf1 = RDB[fib + FINIX_RMAX];

    /* Get clad in and out radii */

    rc0 = RDB[fib + FINIX_CMIN];

    rc1 = RDB[fib + FINIX_CMAX];

    /* Get number of fuel segments */

    nf = (long)RDB[fib + FINIX_NR_FUEL];

    /* Get number of clad segments*/

    nc = (long)RDB[fib + FINIX_NR_CLAD];

    /* Get number of axial segments */

    nz = (long)RDB[fib + FINIX_NZ];

    /* Get pointer to axial segment list */

    axi = (long)RDB[fib + FINIX_PTR_AX];

    /* Initialize the NNODES array */
    nnodes = (int*)finix_get_default_nnodes();

    /* Set number of nodes*/
    nnodes[0] = (int)nz;
    nnodes[1] = (int)nf+1;
    nnodes[2] = (int)nc+1;

    /* Initialize arrays for use */

    err = finix_initialize_arrays(&nnodes,&T,&r,&r_cold,&Pden,&Lhr,&bu,&params,&bcond,&sresults, &vresults);
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
	/* Get user defined values from file */

	fprintf(out, "Reading parameters for FINIX rod from \"%s\"\n",GetText(fib + FINIX_PTR_FNAME));

	err = finix_read_params_file(GetText(fib + FINIX_PTR_FNAME), params, sresults);
	finix_append_err(&allerr, &err);
      }

    /* Check that the radial dimensions are consistent */

    /* Check pellet inner radius */
    if(fabs(params[0][0] - rf0*0.01) > 1e-5)
      Error(0, "Inconsistent pellet inradius for FINIX pin %s (%E vs %E)\n", 
	    GetText(fib + FINIX_PTR_UNI_NAME), params[0][0], rf0*0.01);

    /* Check pellet outer radius */
    if(fabs(params[0][1] - rf1*0.01) > 1e-5)
      Error(0, "Inconsistent pellet outradius for FINIX pin %s (%E vs %E) \n", 
	    GetText(fib + FINIX_PTR_UNI_NAME), params[0][1], rf1*0.01);

    /* Check cladding inner radius */
    if(fabs(params[0][2] - rc0*0.01) > 1e-5)
      Error(0, "Inconsistent cladding inradius for FINIX pin %s (%E vs %E)\n", 
	    GetText(fib + FINIX_PTR_UNI_NAME), params[0][2], rc0*0.01);

    /* Check cladding outer radius */
    if(fabs(params[0][3] - rc1*0.01) > 1e-5)
      Error(0, "Inconsistent cladding outradius for FINIX pin %s (%E vs %E)\n", 
	    GetText(fib + FINIX_PTR_UNI_NAME), params[0][3], rc1*0.01);

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

    /*finix_printf_params(params);*/

    /* Set the node positions */
    /* Equally spaced nodes for finix */

    for(i=0; i< nnodes[0]; i++)
      {

	/* Set pellet node positions */
	for(j=0; j< nnodes[1] ; j++)
	  {
	    r[i][j]=(rf0 + (rf1-rf0)/(double)(nnodes[1]-1)*j)*0.01;
	    r_cold[i][j]=(rf0 + (rf1-rf0)/(double)(nnodes[1]-1)*j)*0.01;
	  }
	
	/* Set cladding node positions */
	for(j=nnodes[1]; j < nnodes[1] + nnodes[2]; j++)
	  {
	    r[i][j]=(rc0 + (rc1-rc0)/(double)(nnodes[2]-1)*
		     (double)(j-nnodes[1]))*0.01;
	    r_cold[i][j]=(rc0 + (rc1-rc0)/(double)(nnodes[2]-1)*
			  (double)(j-nnodes[1]))*0.01;
	  }
      }

    /* Set burnups */

    for(i=0; i< nnodes[0]; i++)
      {

	/* Set pellet node burnups */
	for(j=0; j< nnodes[1] ; j++)
	  {
	    bu[i][j] = TFBU[j];
	  }
	
      }


    /* Set axial length of whole rod */
    if(nz>1){
      params[0][4] = 0.01*(RDB[fib + FINIX_ZMAX] - RDB[fib + FINIX_ZMIN]);
      params[0][5] = 0.01*(RDB[fib + FINIX_ZMAX] - RDB[fib + FINIX_ZMIN]);
    }

    /**************************************************************************/

    /* Set boundary conditions  */

    /* Get BC type */

    options[0] = (long)RDB[fib + FINIX_BCTYPE];

    /* Get pointer to axial zones */

    axi=(long)RDB[fib + FINIX_PTR_AX];

    /* Sort axial zones (probably already sorted) */

    SortList(axi,FINIX_AX_ZMIN,SORT_MODE_ASCEND);

    /* Loop over the zones to set local BC */

    axi=FirstItem(axi);

    /* Cladding surface temperature given */

    if((long)RDB[fib + FINIX_BCTYPE] == FINIX_BCTYPE_CLADT)
      {
	for(i=0;i<nnodes[0];i++)
	  {
	    bcond[0][i] = RDB[axi + FINIX_AX_CLADT];
	    axi=NextItem(axi);
	  }
      }
    /* Heat flux from cladding to coolant given */
    else if((long)RDB[fib + FINIX_BCTYPE] == FINIX_BCTYPE_HFLUX)
      {
	for(i=0;i<nnodes[0];i++)
	  {
	    bcond[3][i] = RDB[axi + FINIX_AX_HFLUX];
	    axi=NextItem(axi);
	  }
    }
    /* Heat transfer coefficient between cladding and coolant given */
    else if((long)RDB[fib + FINIX_BCTYPE] == FINIX_BCTYPE_HTCOE)
      {
	for(i=0;i<nnodes[0];i++)
	  {
	    bcond[2][i] = RDB[axi + FINIX_AX_HTCOE];
	    bcond[1][i] = RDB[axi + FINIX_AX_COOLT];
	    axi=NextItem(axi);
	  }
      }
    /* Coolant temperature given (FINIX calculates heat transfer) */
    else if((long)RDB[fib + FINIX_BCTYPE] == FINIX_BCTYPE_COOLT)
      {
	for(i=0;i<nnodes[0];i++)
	  {
	    bcond[1][i] = RDB[axi + FINIX_AX_COOLT];
	    axi=NextItem(axi);
	  }      
      }
    else
      {
	Die(FUNCTION_NAME,"Unknown boundary condition for temperature feedback\n");
      }

    /*
    finix_printf_params(params);
    
    printf("\n");
    finix_printf_array(nnodes, r);
    finix_printf_array(nnodes, r_cold);    
    printf("\n");
    */

    /* Solve initial steady state (HZP) */
 
    err = finix_solve_initial_steady_state(nnodes,T,r,r_cold,Pden,Lhr,bu,
					   params,bcond,sresults,vresults,
					   options);

    /* Print out any errors */

    if(err!=NULL)
      finix_printf_err(err);

    finix_free_err(&err);

    /*
    finix_printf_array(nnodes,bu);
    finix_printf_array(nnodes,Pden);
    finix_printf_array(nnodes,T);
    */
    /*
    printf("\n");
    */

    /* Store pointers */
    WDB[fib + FINIX_PTR_NNODES] = (long)nnodes;
    WDB[fib + FINIX_PTR_T] = (long)T;
    WDB[fib + FINIX_PTR_R] = (long)r;
    WDB[fib + FINIX_PTR_R_COLD] = (long)r_cold;
    WDB[fib + FINIX_PTR_PDEN] = (long)Pden;
    WDB[fib + FINIX_PTR_LHR] = (long)Lhr;
    WDB[fib + FINIX_PTR_BU] = (long)bu;
    WDB[fib + FINIX_PTR_PARAMS] = (long)params;
    WDB[fib + FINIX_PTR_BCOND] = (long)bcond;
    WDB[fib + FINIX_PTR_SRESULTS] = (long)sresults;
    WDB[fib + FINIX_PTR_VRESULTS] = (long)vresults;
    WDB[fib + FINIX_PTR_OPTIONS] = (long)options;

    /* Increment the number of pins */

    npins++;

    /* Process next pin */

    fib = NextItem(fib);
  }

  /* Store total number of pins (needed for interface file) */

  fib = RDB[DATA_PTR_FIN0];

  while(fib > VALID_PTR)
    {
      WDB[fib + FINIX_N_PIN] = (double)npins;
      fib = NextItem(fib);
    }

  /* Set Density Factor usage on */

  WDB[DATA_USE_DENSITY_FACTOR] = (double)YES;

  fprintf(out, "OK.\n\n");

  /* Read initial conditions */

  ReadFinixIFC();

  /* Write pins to an interface file */

  WriteFinixIFC();

  CreateFinixIFC();

}

#endif

/*****************************************************************************/
