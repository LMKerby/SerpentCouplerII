#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readdecayfile.c                                */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2014/03/07 (JLe)                                           */
/* Version:       2.1.19                                                     */
/*                                                                           */
/* Description: Reads ENDF format decay data into ACE block                  */
/*                                                                           */
/* Comments: - Muuta ton ACE-blokin nimi kun se ei sovi                      */
/*           - Converted from Serpent 1.1.14 readdirectoryfile.c             */
/*           - ENDF decay modes are converted into MT = 10000 + RTYP         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadDecayFile:"

/*****************************************************************************/

void ReadDecayFile()
{
  long n, m, i, ptr, NC, NDK, RFS, RTYP1, Z, A, ZA, LIS, I, NSP, STYP, LCON;
  long NER, NT, LCOV, NR, NP, loc0, loc1, loc2, loc3, ace, pta, count, RTYP2;
  long RTYP3, RTYP4, RTYP5, rt;
  char tmpstr[MAX_STR], line[MAX_STR], *eof;
  double AW, AWR, Q, Qtot, BR, FD, ER, FC, lambda, E, RI;
  FILE *fp;

  /* Check pointer */
  
  if ((pta = (long)RDB[DATA_PTR_DECDATA_FNAME_LIST]) < 1)
    return;
 
  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
    fprintf(out, "Reading decay data files...\n");
  
  /* Reset data pointer */

  ace = -1;

  /* Loop over files */

  while ((long)RDB[pta] > 0)
    {      
      /* Reset counter */

      count = 0;

      /* Test format */

      WDB[DATA_DUMMY] = RDB[pta];
      TestDOSFile(GetText(DATA_DUMMY));

      /* Open file for reading */
  
      fp = OpenDataFile(pta, "decay data");

      /* Read data */
  
      do 
	{
	  /*******************************************************************/

	  /***** Loop over endf file *****************************************/

	  /* Loop until in comment section */
      
	  do 
	    eof = fgets(line, 82, fp);
	  while ((atoi(&line[71]) != 1451) && (eof != NULL));
	  
	  /* Check EOF */
	  
	  if (eof == NULL)
	    break;
	  
	  /* Loop until in decay data */
	  
	  do 
	    eof = fgets(line, 82, fp);
	  while ((atoi(&line[71]) != 8457) && (eof != NULL));
	  
	  /* Check EOF */
	  
	  if (eof == NULL)
	    break;
	  
	  /* Allocate memory for data */
	  
	  ace = ReallocMem(ACE_ARRAY, ACE_BLOCK_SIZE);
	  
	  /* Set null pointer to next */
	  
	  ACE[ace + ACE_PTR_NEXT] = NULLPTR;
	  
	  /* Check if previous exists (ei voi k�ytt�� VALID_PTR) */
	  
	  if ((ptr = (long)RDB[DATA_PTR_ACE0]) < 1)
	    {
	      /* First definition, set pointer */
	      
	      WDB[DATA_PTR_ACE0] = (double)ace;
	    }
	  else
	    {
	      /* Find last block  (tohon ei VALID_PTR) */
	      
	      while ((long)ACE[ptr + ACE_PTR_NEXT] > 0)
		ptr = (long)ACE[ptr + ACE_PTR_NEXT];
	      
	      /* Set pointer to new */
	      
	      ACE[ptr + ACE_PTR_NEXT] = (double)ace;
	    }
	  
	  /* Set pointer to decay ACE data */
	  
	  if ((long)RDB[DATA_PTR_DECAY_ACE0] < 1)
	    WDB[DATA_PTR_DECAY_ACE0] = (double)ace;
	  
	  /* Set type */
	  
	  ACE[ace + ACE_TYPE] = (double)NUCLIDE_TYPE_DECAY;

	  /* Add count */

	  count++;
	  
	  /* Read ZA */
	  
	  ZA = (long)rint(ENDFColF(1, line));
	  
	  /* Check (some libraries have neutron, ZA = 1) */
	  
	  CheckValue(FUNCTION_NAME, "ZA", "", ZA, 1, 120000);
	  
	  /* Put value */
	  
	  ACE[ace + ACE_ZA] = (double)ZA;
	  
	  /* Separate Z and A (not used but needed for bug check) */
	  
	  Z = (long)((double)ZA/1000.0);
	  A = ZA - 1000*Z;
	  
	  /* Check values */
	  
	  CheckValue(FUNCTION_NAME, "Z", "", Z, 0, 120);
	  CheckValue(FUNCTION_NAME, "A", "", A, 0, 290);
	  
	  /* Read atomic weight ratio. */
	  
	  AWR = ENDFColF(2, line);
	  CheckValue(FUNCTION_NAME, "AWR", "", AWR, 0.9, 290.0);
	  
	  /* Put value */
	  
	  ACE[ace + ACE_AWR] = AWR;
	  
	  /* Atomic weight */
	  
	  AW = AWR*M_NEUTRON;
	  
	  /* Check for some known library bugs */
	  
	  if ((Z == 1) && ((AW/A - 1.0) < 1E-3))
	    {
	      Warn(FUNCTION_NAME, "Nuclide %ld mass should not be %1.5f!!!\n\nNOTE: Some versions of publicly available decay libraries (such as\n      ENDF/B-VI and JEFF-3.1) have erroneous AWR data. Check your\n      data source.", ZA, AW);
	      
	      /* Set warning flag */
	      
	      WDB[DATA_WARN_ERROR_AWR] = 1.0;
	    }
	  
	  /* Read excited state number (not used) */
	  
	  LIS = ENDFColI(3, line);
	  
	  /* Read isomeric state number */
	  
	  I = ENDFColI(4, line);
	  ACE[ace + ACE_I] = (double)I;
	  
	  /* Set ZAI */
	  
	  ACE[ace + ACE_ZAI] = (double)(10*ZA + I);

	  /* Check that LIS == 0 if I == 0 */
	  
	  if ((I == 0) && (LIS != 0))
	    Die(FUNCTION_NAME, "I = %ld, LIS = %ld", I, LIS);
	  
	  /* Set name and alias */
	  
	  sprintf(tmpstr, "%ld", 10*ZA + I);
	  ACE[ace + ACE_PTR_NAME] = (double)PutText(tmpstr);
	  ACE[ace + ACE_PTR_ALIAS] = ACE[ace + ACE_PTR_NAME];
	  
	  /* Temperature */
	  
	  ACE[ace + ACE_TEMP] = -1.0;
	  
	  /* Library id */
	  
	  ACE[ace + ACE_PTR_LIB_ID] = PutText("dec");
	  
	  /* Read number of radiation types */
	  
	  NSP = ENDFColI(6, line);
	  
	  /* Read half-life */
	  
	  ENDFNewLine(FUNCTION_NAME, line, fp);
	  lambda = ENDFColF(1, line);
	  
	  /* Convert to decay constant */
	  
	  if (lambda > 0.0)
	    lambda = LOG2/lambda;
	  else
	    lambda = 0.0;
	  
	  /* Check value */
	  
	  CheckValue(FUNCTION_NAME, "lambda", "", lambda, 0.0, 1E+30);
	  
	  /* Set value */

	  ACE[ace + ACE_LAMBDA] = lambda;
	  
	  /* Set Xe-135, I-135 and Pm-149 decay constants */
	  
	  if ((long)ACE[ace + ACE_ZAI] == 541350)
	    WDB[DATA_XE135_DC] = lambda;
	  else if ((long)ACE[ace + ACE_ZAI] == 531350)
	    WDB[DATA_I135_DC] = lambda;
	  else if ((long)ACE[ace + ACE_ZAI] == 611490)
	    WDB[DATA_PM149_DC] = lambda;
	  
	  /* Read number of decay energies */
	  
	  NC = ENDFColI(5, line)/2;
	  
	  if ((NC == 3) || (NC == 17))
	    {
	      /* Read total decay heat */
	      
	      ENDFNewLine(FUNCTION_NAME, line, fp);
	      
	      E = ENDFColF(1, line);
	      E = E + ENDFColF(3, line);
	      E = E + ENDFColF(5, line);
	      
	      /* Convert to MeV */
	      
	      E = E*1E-6;
	      CheckValue(FUNCTION_NAME, "decay energy", "", E, 0.0, 200.0);
	      
	      /* Set value */
	      
	      ACE[ace + ACE_DECAY_E] = E;
	      
	      /* Skip lines */
	      
	      if (NC == 17)
		{
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		}
	    }
	  else
	    Die(FUNCTION_NAME, "Error in decay energies: NC = %ld (ZAI = %ld)", 
		NC, (long)ACE[ace + ACE_ZAI]);
	  
	  /* Read number of decay modes */
	  
	  ENDFNewLine(FUNCTION_NAME, line, fp);
	  NDK = ENDFColI(6, line);
	  CheckValue(FUNCTION_NAME, "NDK", "", NDK, 0, 5);
	  
	  if (NDK > 0)
	    {
	      /* Allocate memory for pointer array */
	      
	      loc0 = ReallocMem(ACE_ARRAY, NDK + 1);
	      
	      /* Set pointer */
	      
	      ACE[ace + ACE_PTR_DECAY_LIST] = (double)loc0;
	      
	      /* Reset total Q-value */
	      
	      Qtot = 0.0;
	      
	      /* Loop over decay modes */
	      
	      m = 0;
	      
	      for (n = 0; n < NDK; n++)
		{
		  /* Read line */
		  
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		  
		  /* Read decay mode */

		  rt = (long)(100000.0*ENDFColF(1, line));

		  /* Get first type */

		  RTYP1 = (long)((double)rt/100000.0);
		  CheckValue(FUNCTION_NAME, "RTYP1", "", RTYP1, 1, 7);

		  /* Parse second */

		  rt = rt - 100000*RTYP1;

		  /* Get second type */

		  RTYP2 = (long)((double)rt/10000.0);
		  CheckValue(FUNCTION_NAME, "RTYP2", "", RTYP2, 0, 7);

		  /* Parse third */

		  rt = rt - 10000*RTYP2;

		  /* Get third type */

		  RTYP3 = (long)((double)rt/1000.0);
		  CheckValue(FUNCTION_NAME, "RTYP3", "", RTYP3, 0, 7);

		  /* Parse fourth */

		  rt = rt - 1000*RTYP3;

		  /* Get fourth type */

		  RTYP4 = (long)((double)rt/100.0);
		  CheckValue(FUNCTION_NAME, "RTYP4", "", RTYP4, 0, 7);

		  /* Parse fifth */

		  rt = rt - 100*RTYP4;

		  /* Get fifth type */

		  RTYP5 = (long)((double)rt/10.0);
		  CheckValue(FUNCTION_NAME, "RTYP5", "", RTYP5, 0, 7);
		  
		  /* Check remaining */

		  if (rt - 10*RTYP5 > 0)
		    Note(0, "Nuclide %s has more than 5 successive decays %ld",
			 ZAItoIso(10*ZA + I, 1), rt);
		    
		  /* Set delayed neutron precursor flag */

		  if (RTYP2 == 5)
		    ACE[ace + ACE_DELNU_PREC] = (double)YES;

		  /* Read isomeric state flag */

		  RFS = ENDFColI(2, line);
		  CheckValue(FUNCTION_NAME, "RFS", "", RFS, 0, 3);
		  
		  /* Read Q-value */
		  
		  Q = 1E-6*ENDFColF(3, line);
		  CheckValue(FUNCTION_NAME, "Q", "", Q, 0.0, 200.0);
		  
		  /* Add to total */
		  
		  Qtot = Qtot + Q;
		  
		  /* Read branching ratio */
		  
		  BR = ENDFColF(5, line);
		  CheckValue(FUNCTION_NAME, "BR", "", BR, 0.0, 1.0);
		  
		  /* Store spontaneous fission branching ratio */
		  
		  if (RTYP1 == 6)
		    ACE[ace + ACE_SFBR] = BR;
		

		  if (1 == 2)
		    {
		  
		  if (RTYP2 > 0)
		    Warn(FUNCTION_NAME, "Tassa tehdaan vakivaltaa summaamalla moodit"); 

		  /* Loop over existing reaction modes */
		  
		  for (i = 0; i < m; i++)
		    {
		      /* Get pointer */
		      
		      loc1 = (long)ACE[loc0 + i];
		      
		      /* Compare RTYP and RFS */
		      
		      if ((10000 + RTYP1 == (long)ACE[loc1 + DECAY_RTYP1]) &&
			  (RFS == (long)ACE[loc1 + DECAY_RFS]))
			break;
		    }
		    }
		  else
		    i = m;

		  /* Check if found */
		  
		  if (i < m)
		    {
		      /* Found existing. Add branching ratio */
		      
		      ACE[loc1 + DECAY_BR] = ACE[loc1 + DECAY_BR] + BR;
		    }
		  else
		    {
		      /* Allocate memory for data */
		      
		      loc1 = ReallocMem(ACE_ARRAY, DECAY_BLOCK_SIZE);
		      
		      /* Set pointer */
		      
		      ACE[loc0 + m++] = (double)loc1;
		      
		      /* Check isomeric transition */
		      
		      /* Set values */
		      
		      ACE[loc1 + DECAY_RTYP1] = (double)(10000 + RTYP1);

		      ACE[loc1 + DECAY_RTYP2] = (double)(10000 + RTYP2);
		      ACE[loc1 + DECAY_RTYP3] = (double)(10000 + RTYP3);
		      ACE[loc1 + DECAY_RTYP4] = (double)(10000 + RTYP4);
		      ACE[loc1 + DECAY_RTYP5] = (double)(10000 + RTYP5);
		      ACE[loc1 + DECAY_RFS] = (double)RFS;
		      ACE[loc1 + DECAY_Q] = Q;
		      ACE[loc1 + DECAY_BR] = BR;
		    }
		}
	      
	      /* Set null pointer */
	      
	      ACE[loc0 + m] = NULLPTR;
	      
	      /* Allocate memory for radiation spectra list */
	      
	      loc0 = ReallocMem(ACE_ARRAY, NSP + 1);
	      
	      /* Put pointer */
	      
	      ACE[ace + ACE_PTR_RAD_SPEC] = (double)loc0;
	      
	      /* Loop over radiation types */
	      
	      for (n = 0; n < NSP; n++)
		{
		  /* Allocate memory for radiation type */
		  
		  loc1 = ReallocMem(ACE_ARRAY, RAD_SPEC_BLOCK_SIZE);
		  
		  /* Put pointer */
		  
		  ACE[loc0++] = (double)loc1;
		  
		  /* Next line */
		  
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		  
		  /* Read radiation type */
		  
		  STYP = ENDFColI(2, line);
		  CheckValue(FUNCTION_NAME, "STYP", "", STYP, 0, 9);
		  
		  /* Put type */
		  
		  ACE[loc1 + RAD_SPEC_TYPE] = (double)STYP;
		  
		  /* Read continuum spectrum flag */
		  
		  LCON = ENDFColI(3, line);
		  
		  /* Read number of discrete energies */
		  
		  NER = ENDFColI(6, line);
		  
		  /* Put number of energies */
		  
		  ACE[loc1 + RAD_SPEC_DISC_NE] = (double)NER;
		  
		  /* Next line */
		  
		  ENDFNewLine(FUNCTION_NAME, line, fp);
		  
		  /* Read discrete spectrum normalization factor */
		  
		  FD = ENDFColF(1, line);
		  ACE[loc1 + RAD_SPEC_DISC_NORM] = (double)FD;
		  
		  /* Read average radiation energy */
		  
		  ER = ENDFColI(3, line);
		  ACE[loc1 + RAD_SPEC_AVG_E] = (double)ER;
		  
		  /* Read continuum spectrum normalization factor */
		  
		  FC = ENDFColI(5, line);
		  ACE[loc1 + RAD_SPEC_CONT_NORM] = (double)FC;
		  
		  /* Check number of discrete lines */
		  
		  if (NER > 0)
		    {
		      /* Allocate memory */
		      
		      loc2 = ReallocMem(ACE_ARRAY, NER);
		      loc3 = ReallocMem(ACE_ARRAY, NER);
		      
		      /* Put pointers */
		      
		      ACE[loc1 + RAD_SPEC_PTR_DISC_E] = (double)loc2;
		      ACE[loc1 + RAD_SPEC_PTR_DISC_RI] = (double)loc3;
		    }
		  
		  /* Loop over discrete energies */
		  
		  for (m = 0; m < NER; m++)
		    {
		      /* Next line */
		      
		      ENDFNewLine(FUNCTION_NAME, line, fp);
		      
		      /* Read energy */
		      
		      E = ENDFColF(1, line);
		      
		      /* Check energy */
		      
		      if ((E < 1.0) || (E > 25E+6))
			Warn(FUNCTION_NAME, 
			     "Suspicious line spectrum energy %E (ZAI = %ld)",
			     E, 10*ZA + I);

		      /* Convert to MeV */
		      
		      ACE[loc2++] = E*1E-6;
		      
		      /* Number of values */
		      
		      NT = ENDFColI(5, line);
		      
		      /* Next line */
		      
		      ENDFNewLine(FUNCTION_NAME, line, fp);
		      
		      /* Read relative intensity */
		      
		      RI = ENDFColF(3, line);
		      
		      /* Put value */
		      
		      ACE[loc3++] = RI;
		      
		      /* Check for extra lines */
		      
		      if ((NT == 12) || (NT == 10) || (NT == 8))
			ENDFNewLine(FUNCTION_NAME, line, fp);
		      else if ((NT != 6) && (NT != 4))
			Die(FUNCTION_NAME, "Invalid number of entries");
		    }
		  
		  /* Check for continuous spectrum */
		  
		  if (LCON != 0)
		    {
		      /* Next line */
		      
		      ENDFNewLine(FUNCTION_NAME, line, fp);
		      
		      /* Read covariance data flag */
		      
		      LCOV = ENDFColI(4, line);
		      
		      /* Check covariance parameter */
		      
		      if (LCOV != 0)
			Die(FUNCTION_NAME, "Nuclide %ld has covariance data",
			    (long)ACE[ace + ACE_ZAI]);
		      
		      /* Read number of interpolation ranges */
		      
		      NR = ENDFColI(5, line);
		      
		      /* Check */
		      
		      if (NR != 1)
			Die(FUNCTION_NAME, 
			    "Nuclide %ld has interpolation ranges",
			    (long)ACE[ace + ACE_ZAI]);
		      
		      /* Read number of points */
		      
		      NP = ENDFColI(6, line);
		      CheckValue(FUNCTION_NAME, "NP", "", NP, 2, INFTY);
		      
		      /* Allocate memory */
		      
		      loc2 = ReallocMem(ACE_ARRAY, NP);
		      loc3 = ReallocMem(ACE_ARRAY, NP);
		      
		      /* Put pointers */
		      
		      ACE[loc1 + RAD_SPEC_PTR_CONT_E] = (double)loc2;
		      ACE[loc1 + RAD_SPEC_PTR_CONT_RI] = (double)loc3;
		      
		      /* Next line */
		      
		      ENDFNewLine(FUNCTION_NAME, line, fp);
		      ENDFNewLine(FUNCTION_NAME, line, fp);
		      
		      /* Loop over values */
		      
		      m = 0;
		      i = 1;
		      
		      while (1 != 2)
			{
			  /* Read values */
			  
			  E = ENDFColF(i++, line);
			  RI = ENDFColF(i++, line);
			  
			  /* Check energy */
			  
			  CheckValue(FUNCTION_NAME, "E", "", E, 0.0, 35E+6);
			  
			  /* Put values */
			  
			  ACE[loc2++] = E*1E-6;
			  ACE[loc3++] = RI;
			  
			  /* Increase counter */
			  
			  m++;
			  
			  /* Check number of values */
			  
			  if (m == NP)
			    break;
			  
			  /* Check column index */
			  
			  if (i > 6)
			    {
			      /* Read new line */
			      
			      ENDFNewLine(FUNCTION_NAME, line, fp);
			      
			      /* Reset counter */
			      
			      i = 1;
			    }
			}
		    }
		}
	      
	      /* Put null pointer */
	      
	      ACE[loc0] = NULLPTR;
	      
	      /* Set number of decay modes */
	      
	      ACE[ace + ACE_DECAY_NDK] = (double)m;
	      
	      /* Check total branching ratio */
	      
	      BR = 0.0;
	      
	      /* Loop over reactions */
	      
	      loc0 = (long)ACE[ace + ACE_PTR_DECAY_LIST];
	      
	      while ((loc1 = (long)ACE[loc0++]) > 0)
		BR = BR + ACE[loc1 + DECAY_BR];
	      
	      if (fabs(BR - 1.0) > 5.0E-2)
		{	      
		  Warn(FUNCTION_NAME, 
		       "Sum of decay branching ratios is %f (ZAI = %ld)",
		       BR, 10*ZA + I);
		  
		  /* Add warning counter */
		  
		  WDB[DATA_WARN_ERROR_BRANCH] =
		    RDB[DATA_WARN_ERROR_BRANCH] + 1.0;
		}
	    }
	  
	  /*******************************************************************/
	}
      while (1 != 2);
      
      /* Close file */
      
      fclose(fp);

      /* Check count */

      if (count == 0)
	Error(0, "No decay data available in file \"%s\"", GetText(pta));

      /* Next file */

      pta++;
    }

  if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
    fprintf(out, "OK.\n\n");
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
