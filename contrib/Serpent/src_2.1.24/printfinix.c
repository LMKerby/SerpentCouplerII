/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printfinix.c                                   */
/*                                                                           */
/* Created:       2013/10/12 (VVa)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Prints some interesting stuff from transient Finix           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "PrintFinix:"

/*****************************************************************************/

void PrintFinix()
{
  long fib, ifc, fpe, tbi, ptr;
  long tb,i,j;
  int *nnodes;
  double **T;
  double **r;
  double **r_cold;
  double **Pden;
  double *Lhr;
  double **bcond;
  char *name, outfile[MAX_STR];
  FILE* fp;

  /* Return if no finix rods are defined */

  if((fib = (long)RDB[DATA_PTR_FIN0]) < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Get pointer to interface block */
  
  ifc = (long)RDB[fib + FINIX_PTR_IFC];
  CheckPointer(FUNCTION_NAME, "(ifc)", DATA_ARRAY, ifc);

  /* Print output file name */

  sprintf(outfile, "%s_fin.m", GetText(DATA_PTR_INPUT_FNAME));

  /* Open output file */

  fp = fopen(outfile,"w");

  /* Loop over pins */

  fpe = (long)RDB[ifc + IFC_PTR_FUEP];

  while(fpe > VALID_PTR)
    {

      /* Get pointer to time bins */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      /* Get pointer to pin names */

      ptr = (long)RDB[fpe + IFC_FUEP_PTR_UNI_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get first uni name and use it for this pin */
      
      name = GetText(ptr);
      
      /* Print timestep end values */

      fprintf(fp,"t=[");

      while(tbi > VALID_PTR)
	{
	  /* Print value */

	  fprintf(fp,"%E ",RDB[tbi + IFC_FUEP_T_TMAX]);

	  /* Next time bin */
	  tbi = NextItem(tbi);
	}

      fprintf(fp,"];\n");

      /* Reset time bin index  */

      tb = 0;

      /* Loop over time bins to print temperatures */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      while(tbi > VALID_PTR)
	{
	  /* Print array starting part */

	  fprintf(fp,"T(%ld,%s,:,:)=[",tb+1,name);

	  /* Get Finix pointers for the time step*/
	  nnodes   = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  T        = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_T]);

	  /* Loop over axial segments */

	  for(i = 0; i < nnodes[0]; i++)
	    {
	      /* Loop over radial nodes to print temperature */

	      for(j=0; j < nnodes[1] + nnodes[2]; j++)
		fprintf(fp,"%E ",T[i][j]);

	      fprintf(fp,"\n");
	    }
	  fprintf(fp,"];\n");

	  /* Increment time bin counter */

	  tb++;

	  /* Next time bin */
	  
	  tbi = NextItem(tbi);
	}

      fprintf(fp,"\n");

      /* Reset time bin index  */

      tb = 0;

      /* Loop over time bins to print powers */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      while (tbi > VALID_PTR)
	{
	  /* Print array starting part */

	  fprintf(fp,"P(%ld,%s,:,:)=[",tb+1, name);

	  /* Get Finix pointers for the time step*/
	  nnodes   = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  r        = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R]);
	  Pden     = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_PDEN]);
	  Lhr      = (double*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_LHR]);

	  /* Loop over axial segments */

	  for (i = 0; i < nnodes[0]; i++)
	    {
	      /* Loop over radial nodes */

	      for (j=0; j < nnodes[1] + nnodes[2]; j++)
		fprintf(fp,"%E ",Pden[i][j]*Lhr[i]/PI/(r[i][nnodes[1]]*r[i][nnodes[1]])/1E6);
	      fprintf(fp,"\n");
	    }

	  fprintf(fp,"];\n");

	  /* Increment time bin counter */

	  tb++;

	  /* Next time bin */
	  
	  tbi = NextItem(tbi);
	}
      fprintf(fp,"\n");

      /* Reset time bin index  */

      tb = 0;

      /* Loop over time bins to print hot radii */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      while (tbi > VALID_PTR)
	{

	  /* Open array */

	  fprintf(fp,"RH(%ld,%s,:,:)=[",tb+1,name);

	  /* Get Finix pointers for the time step*/
	  nnodes   = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  r        = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R]);

	  /* Loop over axial segments */

	  for(i = 0; i < nnodes[0]; i++)
	    {
	      
	      /* Loop over radial nodes */

	      for(j=0; j < nnodes[1] + nnodes[2]; j++)
		fprintf(fp,"%E ",r[i][j]);

	      fprintf(fp,"\n");
	    }

	  fprintf(fp,"];\n");

	  /* Increment time bin counter */

	  tb++;

	  /* Next time bin */
	  
	  tbi = NextItem(tbi);

	}

      fprintf(fp,"\n");

      /* Reset time bin index  */

      tb = 0;

      /* Loop over time bins to print cold radii */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      while (tbi > VALID_PTR)
	{

	  /* Open array */

	  fprintf(fp,"RC(%ld,%s,:,:)=[",tb+1,name);

	  /* Get Finix pointers for the time step*/
	  nnodes   = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  r_cold   = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_R_COLD]);

	  /* Loop over axial segments */

	  for(i = 0; i < nnodes[0]; i++)
	    {
	      
	      /* Loop over radial nodes */

	      for(j=0; j < nnodes[1] + nnodes[2]; j++)
		fprintf(fp,"%E ",r_cold[i][j]);

	      fprintf(fp,"\n");
	    }

	  fprintf(fp,"];\n");

	  /* Increment time bin counter */

	  tb++;

	  /* Next time bin */
	  
	  tbi = NextItem(tbi);

	}

      fprintf(fp,"\n");

      /* Reset time bin index  */

      tb = 0;

      /* Loop over time bins to print heat flux to coolant */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      /* Open array */

      fprintf(fp,"HFlux(%s,:,:)=[",name);

      while(tbi > VALID_PTR)
	{
	  /* Get Finix pointers for the time step*/
	  nnodes   = (int*)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_NNODES]);
	  bcond    = (double**)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BCOND]);

	  /* Loop over axial segments */

	  for(i = 0; i < nnodes[0]; i++)
	    {
	      fprintf(fp,"%E ",bcond[3][i]);
	      fprintf(fp,"\n");
	    }
	  
	  /* Increment time bin counter */

	  tb++;

	  /* Next time bin */

	  tbi = NextItem(tbi);

	}
      fprintf(fp,"];\n\n");

      /* Next rod */

      fpe = NextItem(fpe);
    }

  fclose(fp);
}

#endif

/*****************************************************************************/
