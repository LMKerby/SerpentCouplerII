/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : aresoutput.c                                   */
/*                                                                           */
/* Created:       2013/03/01 (JLe)                                           */
/* Last modified: 2013/03/04 (JLe)                                           */
/* Version:       2.1.13                                                     */
/*                                                                           */
/* Description: Writes output for the ARES reactor simulator code            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ARESOutput:"

/*****************************************************************************/

void ARESOutput()
{
#ifdef mmmmmmmmmmmmmmmmmmmmm
  long loc0, ptr, n, m, rep;
  double val;
  FILE *fp;
  char outfile[MAX_STR];

  return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* Open file for writing */

  sprintf(outfile, "%s.cs", GetText(DATA_PTR_INPUT_FNAME));
  if ((fp = fopen(outfile, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Print nominal values */

  fprintf(fp, "         DEN          TFU          TMO\n");
  fprintf(fp, "%12.5E %12.5E %12.5E\n\n", 0.0, 0.0, 0.0);

  /* Tän tarvii ehkä kikkailuun */

  rep = 1;

  /***************************************************************************/

  /****** Absorption cross sections ******************************************/

  /* Fast */

  fprintf(fp, "   SIGA1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal */

  fprintf(fp, "   SIGA2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Fission neutron production cross sections **************************/

  /* Fast */

  fprintf(fp, "   NUSF1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_NSF];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_NSF];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal */

  fprintf(fp, "   NUSF2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_NSF];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_NSF];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Diffusion coefficients *********************************************/

  /* Fast */

  fprintf(fp, "   DIFF1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_DIFF];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_DIFF];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal */

  fprintf(fp, "   DIFF2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_DIFF];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_DIFF];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Down-scattering cross section **************************************/

  /* Single value */

  fprintf(fp, "   SIG12\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_SCATTXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_SCATTXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 2]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Fission energy production cross section ****************************/

  /* Fast */

  fprintf(fp, "   EPSF1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Get fission cross section */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];
	  
	  val = RDB[ptr];

	  /* Multiply by fission energy */
	  
	  ptr = (long)RDB[loc0 + SIM_INF_FISSE];
	  val = val*RDB[ptr]*MEV;

	  /* Print base value */

	  fprintf(fp, "%12.5E ", val);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal */

  fprintf(fp, "   EPSF2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);

	  /* Get fission cross section */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];
	  
	  val = RDB[ptr + 1];

	  /* Multiply by fission energy */
	  
	  ptr = (long)RDB[loc0 + SIM_INF_FISSE];
	  val = val*RDB[ptr + 1]*MEV;

	  /* Print base value */

	  fprintf(fp, "%12.5E ", val);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Fission cross sections *********************************************/

  /* Fast */

  fprintf(fp, "   SIGF1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal */

  fprintf(fp, "   SIGF2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Surface discontinuity factors **************************************/

  /* West, fast */

  fprintf(fp, "   SDF1W\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* West, thermal */

  fprintf(fp, "   SDF2W\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* South, fast */

  fprintf(fp, "   SDF1S\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr + 2]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* South, thermal */

  fprintf(fp, "   SDF2S\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr + 3]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* East, fast */

  fprintf(fp, "   SDF1E\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr + 4]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* East, thermal */

  fprintf(fp, "   SDF2E\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr + 5]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* North, fast */

  fprintf(fp, "   SDF1N\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr + 6]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* North, thermal */

  fprintf(fp, "   SDF2N\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFS];
	  fprintf(fp, "%12.5E ", RDB[ptr + 7]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Corner discontinuity factors ***************************************/

  /* South-West, fast */

  fprintf(fp, "   CD1SW\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* South-West, thermal */

  fprintf(fp, "   CD2SW\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Southe-East, fast */

  fprintf(fp, "   CD1SE\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr + 2]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Southe-East, thermal */

  fprintf(fp, "   CD2SE\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr + 3]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* North-East, fast */

  fprintf(fp, "   CD1NE\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr + 4]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* North-East, thermal */

  fprintf(fp, "   CD2NE\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr + 5]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* North-West, fast */

  fprintf(fp, "   CD1NW\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr + 6]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* North-West, thermal */

  fprintf(fp, "   CD2NW\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_ADFC];
	  fprintf(fp, "%12.5E ", RDB[ptr + 7]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** I-135 production and absorption ************************************/

  /* Fast absorption */

  fprintf(fp, "   SIGI1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_I135ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_I135ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal absorption */

  fprintf(fp, "   SIGI2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_I135ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_I135ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Fission yield */

  fprintf(fp, "   SIGI3\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);

	  /* Get fission cross section (currently yield has no energy */
	  /* dependence, so use thermal value) */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];

	  val = RDB[ptr];
	  
	  /* Print base value */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_I135PRODXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_I135PRODXS];

	  if (val > 0.0)
	    fprintf(fp, "%12.5E ", RDB[ptr]/val);
	  else
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Xe-135 production and absorption ***********************************/

  /* Fast absorption */

  fprintf(fp, "   SIGX1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_XE135ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_XE135ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal absorption */

  fprintf(fp, "   SIGX2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_XE135ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_XE135ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Fission yield */

  fprintf(fp, "   SIGX3\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);

	  /* Get fission cross section (currently yield has no energy */
	  /* dependence, so use thermal value) */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];

	  val = RDB[ptr];
	  
	  /* Print base value */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_XE135PRODXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_XE135PRODXS];

	  if (val > 0.0)
	    fprintf(fp, "%12.5E ", RDB[ptr]/val);
	  else
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Pm-149 production and absorption ***********************************/

  /* Fast absorption */

  fprintf(fp, "   SIGP1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_PM149ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_PM149ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal absorption */

  fprintf(fp, "   SIGP2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_PM149ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_PM149ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Fission yield */

  fprintf(fp, "   SIGP3\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);

	  /* Get fission cross section (currently yield has no energy */
	  /* dependence, so use thermal value) */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];

	  val = RDB[ptr];
	  
	  /* Print base value */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_PM149PRODXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_PM149PRODXS];

	  if (val > 0.0)
	    fprintf(fp, "%12.5E ", RDB[ptr]/val);
	  else
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Sm-149 production and absorption ***********************************/

  /* Fast absorption */

  fprintf(fp, "   SIGS1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_SM149ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_SM149ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal absorption */

  fprintf(fp, "   SIGS2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_SM149ABSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_SM149ABSXS];
	  
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Fission yield */

  fprintf(fp, "   SIGS3\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);

	  /* Get fission cross section (currently yield has no energy */
	  /* dependence, so use thermal value) */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_FISSXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_FISSXS];

	  val = RDB[ptr];
	  
	  /* Print base value */

	  if ((long)RDB[DATA_B1_CALC] == NO)
	    ptr = (long)RDB[loc0 + SIM_INF_SM149PRODXS];
	  else
	    ptr = (long)RDB[loc0 + SIM_B1_SM149PRODXS];

	  if (val > 0.0)
	    fprintf(fp, "%12.5E ", RDB[ptr]/val);
	  else
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Beta-eff ***********************************************************/
  
  /* Total */

  fprintf(fp, "   BETA\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_BETA_EFF];
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 1. precursor group  */

  fprintf(fp, "   BETA1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_BETA_EFF];
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 2. precursor group  */

  fprintf(fp, "   BETA2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_BETA_EFF];
	  fprintf(fp, "%12.5E ", RDB[ptr + 2]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 3. precursor group  */

  fprintf(fp, "   BETA3\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_BETA_EFF];
	  fprintf(fp, "%12.5E ", RDB[ptr + 3]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 4. precursor group  */

  fprintf(fp, "   BETA4\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_BETA_EFF];
	  fprintf(fp, "%12.5E ", RDB[ptr + 4]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 5. precursor group  */

  fprintf(fp, "   BETA5\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_BETA_EFF];
	  fprintf(fp, "%12.5E ", RDB[ptr + 5]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 6. precursor group  */

  fprintf(fp, "   BETA6\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_BETA_EFF];
	  fprintf(fp, "%12.5E ", RDB[ptr + 6]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /****** Inverse neutron velocity *******************************************/

  /* Fast */

  fprintf(fp, "   INVV1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_RECIPVEL];
	  fprintf(fp, "%12.5E ", RDB[ptr]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* Thermal */

  fprintf(fp, "   INVV2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */

	  ptr = (long)RDB[loc0 + SIM_RECIPVEL];
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */

	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/






  
  /****** Delayed neutron precursor decay constants **************************/
  
  fprintf(fp, "   LBDA1\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_LAMBDA];
	  fprintf(fp, "%12.5E ", RDB[ptr + 1]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 2. precursor group  */

  fprintf(fp, "   LBDA2\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_LAMBDA];
	  fprintf(fp, "%12.5E ", RDB[ptr + 2]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 3. precursor group  */

  fprintf(fp, "   LBDA3\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_LAMBDA];
	  fprintf(fp, "%12.5E ", RDB[ptr + 3]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 4. precursor group  */

  fprintf(fp, "   LBDA4\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_LAMBDA];
	  fprintf(fp, "%12.5E ", RDB[ptr + 4]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 5. precursor group  */

  fprintf(fp, "   LBDA5\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_LAMBDA];
	  fprintf(fp, "%12.5E ", RDB[ptr + 5]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /* 6. precursor group  */

  fprintf(fp, "   LBDA6\n");

  /* Loop over data */

  loc0 = (long)RDB[DATA_PTR_SIM0];
  while (loc0 > VALID_PTR)
    {
      /* Tän tarvii ehkä kikkailuun */

      for (m = 0; m < rep; m++)
	{
	  /* Print parameters */

	  fprintf(fp, "%2ld %5.1f %5.1f ", (long)RDB[loc0 + SIM_CR], 
		  RDB[loc0 + SIM_VHI], RDB[loc0 + SIM_BURNUP]);
	  
	  /* Print base value */
	  
	  ptr = (long)RDB[loc0 + SIM_LAMBDA];
	  fprintf(fp, "%12.5E ", RDB[ptr + 6]);
	  
	  /* Print interpolation coefficient */
	  
	  for (n = 0; n < 23; n++)
	    fprintf(fp, "%12.5E ", 0.0);

	  /* Print newline */
	  
	  fprintf(fp, "\n");
	}

      /* Next */

      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /* Close file */

  fclose(fp);

#endif
}

/*****************************************************************************/

