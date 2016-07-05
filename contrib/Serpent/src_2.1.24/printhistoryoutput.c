/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printhistoryoutput.c                           */
/*                                                                           */
/* Created:       2011/05/18 (JLe)                                           */
/* Last modified: 2012/10/25 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Prints history of statistical variables                      */
/*                                                                           */
/* Comments: - Tää toimii nyt vaan yksiarvoisilla muuttujilla                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintHistoryOutput:"

/*****************************************************************************/

void PrintHistoryOutput()
{
  long loc0, ptr, i, n0, n1, n2, n3, n4, dim, bins[5];
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Test print option */

  if ((long)RDB[DATA_OPTI_PRINT_HIS] == NO)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* Set file name */

  sprintf(tmpstr, "%s_his%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	  (long)RDB[DATA_BURN_STEP]);
	  
  /* Open file */
  
  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Loop over statistics */
  
  loc0 = (long)RDB[DATA_PTR_SCORE0];
  while(loc0 > VALID_PTR)
    {
      /* Check history pointer */

      if ((long)RDB[loc0 + SCORE_PTR_HIS] > 0)
	{
	  /* Get dimension */
	  
	  if ((dim = (long)RDB[loc0 + SCORE_DIM]) > 5)
	    Die(FUNCTION_NAME, "Need more dimensions");

	  /* Pointer to list of maximum values */
	  
	  ptr = (long)RDB[loc0 + SCORE_PTR_NMAX];
	  
	  /* Read bin sizes */
	  
	  for (i = 0; i < dim; i++)
	    {
	      /* Get number of bins */
	      
	      bins[i] = (long)RDB[ptr++];
	    }

	  /* Print variable name */

	  sprintf(tmpstr, "HIS_%s", GetText(loc0 + SCORE_PTR_NAME));
	  fprintf(fp, "%s = [\n", tmpstr);

	  /* Loop over cycles */

	  for (i = 0; i < (long)RDB[DATA_CYCLE_IDX]; i++)
	    {
	      /* Check skip cycles */

	      if (i == (long)RDB[DATA_CRIT_SKIP])
		fprintf(fp, "\n%% ----- Begin active cycles -----\n\n");

	      if (i < (long)RDB[DATA_CRIT_SKIP])
		fprintf(fp, "%4ld  ", i + 1);
	      else
		fprintf(fp, "%4ld  ", i - (long)RDB[DATA_CRIT_SKIP] + 1);

	      /* Check number of dimensions and print data */

	      if (dim == 1)
		{
		  for (n0 = 0; n0 < bins[0]; n0++)
		    fprintf(fp, "%1.5E %1.5E %1.5f  ", 
			    HisVal(loc0, i, n0),
			    HisMean(loc0, i, n0), 
			    HisRelErr(loc0, i, n0));
		}
	      else if (dim == 2)
		{
		  for (n0 = 0; n0 < bins[0]; n0++)
		    for (n1 = 0; n1 < bins[1]; n1++)
		      fprintf(fp, "%1.5E %1.5E %1.5f  ", 
			      HisVal(loc0, i, n0, n1),
			      HisMean(loc0, i, n0, n1), 
			      HisRelErr(loc0, i, n0, n1));
		}
	      else if (dim == 3)
		{
		  for (n0 = 0; n0 < bins[0]; n0++)
		    for (n1 = 0; n1 < bins[1]; n1++)
		      for (n2 = 0; n2 < bins[2]; n2++)
			fprintf(fp, "%1.5E %1.5E %1.5f  ", 
				HisVal(loc0, i, n0, n1, n2),
				HisMean(loc0, i, n0, n1, n2), 
				HisRelErr(loc0, i, n0, n1, n2));
		}
	      else if (dim == 4)
		{
		  for (n0 = 0; n0 < bins[0]; n0++)
		    for (n1 = 0; n1 < bins[1]; n1++)
		      for (n2 = 0; n2 < bins[2]; n2++)
			for (n3 = 0; n3 < bins[3]; n3++)
			  fprintf(fp, "%1.5E %1.5E %1.5f  ", 
				  HisVal(loc0, i, n0, n1, n2, n3),
				  HisMean(loc0, i, n0, n1, n2, n3), 
				  HisRelErr(loc0, i, n0, n1, n2, n3));
		}
	      else if (dim == 5)
		{
		  for (n0 = 0; n0 < bins[0]; n0++)
		    for (n1 = 0; n1 < bins[1]; n1++)
		      for (n2 = 0; n2 < bins[2]; n2++)
			for (n3 = 0; n3 < bins[3]; n3++)
			  for (n4 = 0; n4 < bins[4]; n4++)
			    fprintf(fp, "%1.5E %1.5E %1.5f  ", 
				    HisVal(loc0, i, n0, n1, n2, n3, n4),
				    HisMean(loc0, i, n0, n1, n2, n3, n4), 
				    HisRelErr(loc0, i, n0, n1, n2, n3, n4));
		}
	      else
		Die(FUNCTION_NAME, "wtf?");

	      fprintf(fp, "\n");
	    }

	  fprintf(fp, "];\n\n");
	}
      
      /* Next */
      
      loc0 = NextItem(loc0);
    }
  
  /* Close file */

  fclose(fp);

  /***************************************************************************/
}

/*****************************************************************************/
