#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printfinix.c                                   */
/*                                                                           */
/* Created:       2013/10/12 (VVa)                                           */
/* Last modified: 2016/02/17 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Prints some interesting stuff from transient Finix           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/finixdata.h"

#define FUNCTION_NAME "PrintFinix:"

/*****************************************************************************/

void PrintFinix()
{
  long fib, ifc, fpe, tbi, ptr;
  long tb,i,j;
  Boundary_conditions *bc;
  Results *results;
  Options *options;
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

  while (fpe > VALID_PTR)
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

      fprintf(fp,"t = [");

      while(tbi > VALID_PTR)
	{
	  /* Print value */

	  fprintf(fp,"%E ",RDB[tbi + IFC_FUEP_T_TMAX]);

	  /* Next time bin */
	  tbi = NextItem(tbi);
	}

      fprintf(fp,"];\n\n");

      /* Reset time bin index  */

      tb = 0;

      /* Loop over time bins to print temperatures */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      while (tbi > VALID_PTR)
	{
	  /* Print array starting part */

	  fprintf(fp,"ROD_%s_T(%ld,:,:) = [",name,tb+1);

	  /* Get Finix pointers for the time step*/
	  
	  results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	  options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
	  bc = (Boundary_conditions *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BC]);

	  /* Loop over axial segments */

	  for (i = 0; i < options->axial_nodes; i++)
	    {
	      /* Loop over radial nodes to print temperature */

	      for ( j = 0 ; j < options->pellet_radial_nodes + 
		      options->clad_radial_nodes ; j++)
		fprintf(fp,"%E ",results->temperature[i][j]);

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

      /* Loop over time bins to print power density */

      tbi = (long)RDB[fpe + IFC_FUEP_PTR_T];
      CheckPointer(FUNCTION_NAME, "(tbi)", DATA_ARRAY, tbi);

      while (tbi > VALID_PTR)
	{
	  /* Print array starting part */

	  fprintf(fp,"ROD_%s_P(%ld,:,:) = [",name,tb+1);

	  /* Get Finix pointers for the time step*/
	  
	  results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	  options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
	  bc = (Boundary_conditions *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BC]);

	  /* Loop over axial segments */

	  for (i = 0; i < options->axial_nodes; i++)
	    {
	      /* Loop over radial nodes  */

	      for ( j = 0 ; j < options->pellet_radial_nodes + 
		      options->clad_radial_nodes ; j++)
		fprintf(fp,"%E ",bc->linear_power[i]/
			(PI*results->radial_node_position[i][options->pellet_radial_nodes]*
			 results->radial_node_position[i][options->pellet_radial_nodes])*
			bc->power_distr[i][j]/1E6);

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

	  fprintf(fp,"ROD_%s_RH(%ld,:,:) = [",name,tb+1);

	  /* Get Finix pointers for the time step*/
	  
	  results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	  options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
	  bc = (Boundary_conditions *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BC]);

	  /* Loop over axial segments */

	  for (i = 0; i < options->axial_nodes; i++)
	    {
	      /* Loop over radial nodes  */

	      for ( j = 0 ; j < options->pellet_radial_nodes + 
		      options->clad_radial_nodes ; j++)
		fprintf(fp,"%E ",results->radial_node_position[i][j]);

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

	  fprintf(fp,"ROD_%s_RC(%ld,:,:) = [",name,tb+1);

	  /* Get Finix pointers for the time step*/
	  
	  results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	  options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
	  bc = (Boundary_conditions *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BC]);
	  
	  /* Loop over axial segments */

	  for (i = 0; i < options->axial_nodes; i++)
	    {
	      /* Loop over radial nodes  */

	      for ( j = 0 ; j < options->pellet_radial_nodes + 
		      options->clad_radial_nodes ; j++)
		fprintf(fp,"%E ",results->radial_node_position_cold[i][j]);

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
      
      fprintf(fp,"ROD_%s_HFlux(:,:) = [",name);

      while(tbi > VALID_PTR)
	{

	  /* Get Finix pointers for the time step*/
	  
	  results = (Results *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_RESULTS]);
	  options = (Options *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_OPTIONS]);
	  bc = (Boundary_conditions *)((long)RDB[tbi + IFC_FUEP_FINIX_PTR_BC]);
	  
	  /* Loop over axial segments */

	  for (i = 0; i < options->axial_nodes; i++)
	    {
	      fprintf(fp,"%E ",bc->heat_flux[i]);
	    }
	  fprintf(fp,"\n");
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
#ifdef __cplusplus 
} 
#endif 
