/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printtfboutput.c                               */
/*                                                                           */
/* Created:       2012/01/24 (JLe)                                           */
/* Last modified: 2013/02/05 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Prints output from temperature feedback calculation          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintTFBOutput:"

/*****************************************************************************/

void PrintTFBOutput()
{
  long tfb, nst, n, m;
  char tmpstr[MAX_STR];
  FILE *fp;
  
  /* Check if feedback is in use */

  if ((long)RDB[DATA_USE_TFB] == NO)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* File name */

  sprintf(tmpstr, "%s_tfb%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	  (long)RDB[DATA_BURN_STEP]);

  /* Open file for writing (append files if burnup calculation) */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) ||
      (((long)RDB[DATA_BURN_STEP] == 0) &&
       ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)))
    fp = fopen(tmpstr, "w");
  else if ((long)RDB[DATA_SIMULATION_COMPLETED] == YES)
    fp = fopen(tmpstr, "a");
  else
    return;

  /* Check pointer */

  if (fp == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Loop over time intervals */

  for (m = 0; m < (long)RDB[DATA_DYN_NB]; m++)
    {
      /* increase counter */
      
      fprintf(fp, "\n%% Increase counter:\n\n");
      
      fprintf(fp, "if (exist(\'idx\', \'var\'));\n");
      fprintf(fp, "  idx = idx + 1;\n");
      fprintf(fp, "else;\n");
      fprintf(fp, "  idx = 1;\n"); 
      fprintf(fp, "end;\n");
      
      /* Print total power and power density */
      
      if (m == 0)
	{
	  fprintf(fp, "\n");
	  
	  fprintf(fp, "%% Total power and power density:\n\n");
	  
	  PrintValues(fp, "TOT_POWER", RES_TOT_NEUTRON_POWER, 1, -1, -1, 0, 0);
	  PrintValues(fp, "TOT_POWDENS", RES_TOT_POWDENS, 1, -1, -1, 0, 0);
	}
      
      /* Loop over temperature feedbacks */

      tfb = (long)RDB[DATA_PTR_TFB0];
      while (tfb > VALID_PTR)
	{
	  /* Pointer to nest */
	  
	  nst = (long)RDB[tfb + TFB_PTR_NST];
	  CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);
	  
	  fprintf(fp, "\n");
	  
	  fprintf(fp, "%% Nest %s:\n\n", GetText(nst + NEST_PTR_NAME));
	  
	  /* Get number of regions */
	  
	  n = (long)RDB[tfb + TFB_N_REG];
	  
	  /* Print power, temperatures, densities and radii */
	  
	  sprintf(tmpstr, "TFB_POWER_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MEAN_POW, n, 1, -1, 0, m);
	  
	  sprintf(tmpstr, "TFB_VTEMP_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MEAN_VTEMP, n, 1, -1, 0, m);
	  
	  sprintf(tmpstr, "TFB_FTEMP_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MEAN_FTEMP, n, 1, -1, 0, m);
	  
	  sprintf(tmpstr, "TFB_TMAX_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MAX_TEMP, n, 1, -1, 0, m);
	  
	  sprintf(tmpstr, "TFB_TMIN_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MIN_TEMP, n, 1, -1, 0, m);
	  
	  sprintf(tmpstr, "TFB_MDENS_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MEAN_MDENS, n, 1, -1, 0, m);
	  
	  sprintf(tmpstr, "TFB_HC_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MEAN_HC, n, 1, -1, 0, m);
	  
	  sprintf(tmpstr, "TFB_RAD_%s", GetText(nst + NEST_PTR_NAME));
	  PrintValues(fp, tmpstr, tfb + TFB_PTR_MEAN_RAD, n + 1, 1, -1, 0, m);
	  
	  /* Next feedback */
	  
	  tfb = NextItem(tfb);
	}
    }

  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
