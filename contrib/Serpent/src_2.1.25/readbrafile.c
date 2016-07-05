/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readbrafile.c                                  */
/*                                                                           */
/* Created:       2011/01/28 (JLe)                                           */
/* Last modified: 2014/02/21 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Reads energy-dependent isomeric branching ratios             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadBRAFile:"

/*****************************************************************************/

void ReadBRAFile()
{
  long pta, ZA, ZAI, LIS, NS, LFS[20], NR, NP, INTT[20], mt, n, m, i, NR0, NP0;
  long loc0, loc1, ptr;
  double *E[20], *f[20], sum;
  char line[MAX_STR], *eof;
  FILE *fp;

  /* Check pointer */
  
  if ((pta = (long)RDB[DATA_PTR_BRADATA_FNAME_LIST]) < 1)
    return;
 
  fprintf(out, "Reading isomeric branching ratios...\n\n");

  /* Loop over files */

  while ((long)RDB[pta] > 0)
    {      
      /* Test format */

      WDB[DATA_DUMMY] = RDB[pta];
      TestDOSFile(GetText(DATA_DUMMY));

      /* Open file for reading */
  
      fp = OpenDataFile(pta, "isomeric branching ratio data");
        
      /* Read data */

      while (1 != 2)
	{
	  /*******************************************************************/

	  /***** Read data from file *****************************************/

	  /* Loop until line 1 of file 9 */
      
	  do 
	    eof = fgets(line, 82, fp);
	  while (((line[71] != '9') || (line[78] != ' ') || (line[79] != '1')) 
		 && (eof != NULL));

	  /* Check EOF */
	  
	  if (eof == NULL)
	    break;

	  /* Read ZA, LIS and NS */
	  
	  ZA = (long)rint(ENDFColF(1, line));
	  LIS = ENDFColI(3, line);
	  NS = ENDFColI(5, line);
	  
	  /* Get mt */
	  
	  mt = atoi(&line[72]);
	  
	  /* Calculate ZAI */

	  ZAI = 10*ZA + LIS;

	  /* Avoid compiler warning */

	  NP = -1;

	  /* Reset previous */

	  NR0 = -1;
	  NP0 = -1;

	  /* Loop over states */
	  
	  for (n = 0; n < NS; n++)
	    {
	      /* Next line */
	      
	      ENDFNewLine(FUNCTION_NAME, line, fp);
	      
	      /* Read LFS, NR and NP */
	      
	      LFS[n] = ENDFColI(4, line);
	      NR = ENDFColI(5, line);
	      NP = ENDFColI(6, line);
	      
	      /* Allocate memory for data */

	      if (NS > 20)
		Die(FUNCTION_NAME, "NS > 20");
	      else if (NP0 < 0)
		for (m = 0; m < NS; m++)
		  {
		    E[m] = Mem(MEM_ALLOC, NP, sizeof(double));
		    f[m] = Mem(MEM_ALLOC, NP, sizeof(double));
		  }

	      /* Compare to previous */

	      if ((NR0 > -1) && (NR0 != NR))
		Die(FUNCTION_NAME, "Mismatch in NR");
	      else
		NR0 = NR;
	    
	      if ((NP0 > -1) && (NP0 != NP))
		{
		  Warn(FUNCTION_NAME, "Multiple interpolation regions");
		  break;
		}
	      else
		NP0 = NP;

	      /* Check NR */
	      
	      if (NR > 3)
		Die(FUNCTION_NAME, "NR = %ld (sotkee ton pisteiden luvun)", 
		    NR);

	      /* Next line */
	      
	      ENDFNewLine(FUNCTION_NAME, line, fp);

	      /* Check number of points */

	      if (NP != ENDFColI(1,line))
		Die(FUNCTION_NAME, "Interpolation regions");

	      /* Read interpolation scheme */

	      INTT[n] = ENDFColI(2, line);

	      /* Next line */

	      ENDFNewLine(FUNCTION_NAME, line, fp);

	      /* Loop over data */
	      
	      i = 1;
	      for (m = 0; m < NP; m++)
		{
		  /* Read energy and fraction */

		  E[n][m] = ENDFColF(i++, line)/1E+6;
		  f[n][m] = ENDFColF(i++, line);
		  
		  /* Check values */

		  CheckValue(FUNCTION_NAME, "E", "", E[n][m], 1E-12, 1000.0);
		  CheckValue(FUNCTION_NAME, "f", "", f[n][m], 0.0, 1.0);

		  if ((i > 6) && (m < NP - 1))
		    {
		      ENDFNewLine(FUNCTION_NAME, line, fp);
		      i = 1;
		    }
		}
	    }

	  /* Check if all regions were processed */

	  if ((n != NS) || (NS < 2))
	    {
	      /* Free memory */

	      for (n = 0; n < NS; n++)
		{
		  Mem(MEM_FREE, E[n]);
		  Mem(MEM_FREE, f[n]);
		}

	      /* Cycle loop */

	      continue;
	    }
	
	  /*******************************************************************/

	  /***** Check and normalize *****************************************/

	  /* Check interpolation schemes */

	  for (n = 1; n < NS; n++)
	    if (INTT[n] != INTT[n - 1])
	      Die(FUNCTION_NAME, "Mismatch in interpolation schemes");

	  /* Check energy arrays */
	      
	  for (n = 1; n < NS; n++)
	    for (m = 0; m < NP; m++)
	      if (fabs(E[n][m]/E[n - 1][m] - 1.0) > 0.2)
		Warn(FUNCTION_NAME, "Mismatch in energies: %E %E",
		     E[n][m], E[n - 1][m]);
	  
	  /* Check probabilities */
	  
	  for (m = 0; m < NP; m++)
	    {
	      /* Calculate sum */
	      
	      sum = 0.0;
	      for (n = 0; n < NS; n++)
		sum = sum + f[n][m];
	      
	      /* Check */
	      
	      if (sum == 0.0)
		Die(FUNCTION_NAME, "Sum of fractions is zero");
	      else if (fabs(sum - 1.0) > 1.1E-3)
		Warn(FUNCTION_NAME, "Sum of fractions: %1.3f", sum);
	      
	      /* Normalize */
	      
	      for (n = 0; n < NS; n++)
		f[n][m] = f[n][m]/sum;
	    }
      
	  /*******************************************************************/

	  /***** Store data **************************************************/

	  /* Handle only (n,gamma) reactions to single excited state for now */

	  if ((mt == 102) && (NS == 2))
	    {
	      /* Print */
	      
	      fprintf(out, "Nuclide %7s, reaction %s\n", 
		      ZAItoIso(ZAI,1), ReactionMT(mt));

	      /* Allocate memory for structure */
	      
	      loc0 = NewItem(DATA_PTR_BRA_LIST, BRA_LIST_BLOCK_SIZE);
	      
	      /* Put data */
	      
	      WDB[loc0 + BRA_LIST_ZAI] = (double)ZAI;
	      WDB[loc0 + BRA_LIST_MT] = (double)mt;
	      WDB[loc0 + BRA_LIST_NS] = (double)NS;
	      
	      /* Loop over states */
	      
	      for (n = 0; n < NS; n++)
		{
		  /* Allocate memory for structure */
		  
		  loc1 = NewItem(loc0 + BRA_LIST_PTR_STATES, 
				 BRA_STATE_BLOCK_SIZE);
	      
		  /* Put data */
		  
		  WDB[loc1 + BRA_STATE_NP] = (double)NP;
		  WDB[loc1 + BRA_STATE_INTT] = (double)INTT[n];
		  WDB[loc1 + BRA_STATE_LFS] = (double)LFS[n];
		  
		  /* Store energy grid */
		  
		  ptr = ReallocMem(DATA_ARRAY, NP);
		  WDB[loc1 + BRA_STATE_PTR_ERG] = (double)ptr;
		  
		  for (m = 0; m < NP; m++)
		    WDB[ptr++] = E[n][m];
		  
		  /* Store fractions */
		  
		  ptr = ReallocMem(DATA_ARRAY, NP);
		  WDB[loc1 + BRA_STATE_PTR_FRAC] = (double)ptr;
		  
		  for (m = 0; m < NP; m++)
		    WDB[ptr++] = f[n][m];
		}
	    }
	  
	  /* Free memory */
	  
	  for (n = 0; n < NS; n++)
	    {
	      Mem(MEM_FREE, E[n]);
	      Mem(MEM_FREE, f[n]);
	    }

	  /*******************************************************************/
	}      

      /* Close file */
      
      fclose(fp);

      /* Next file */

      pta++;
    }

  fprintf(out, "\n");
}

/*****************************************************************************/
