/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdephis.c                                */
/*                                                                           */
/* Created:       2011/06/14 (JLe)                                           */
/* Last modified: 2014/01/23 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Links normalization, sets  pcc modes, etc. for depletion     */
/*              histories                                                    */
/*                                                                           */
/* Comments: TODO: jos normalisaatio on asetettu niin, että neutronivuo on   */
/*                 nolla, askel pitäisi tulkita hajoamisaskeleeksi           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessDepHis:"

/*****************************************************************************/

void ProcessDepHis()
{
  long dep, norm, n, ptr, ntot;

  /* Check burnup mode */

 if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
   return;

 /* Check pointer to depletion histories */

 if ((long)RDB[DATA_BURN_PTR_DEP] < VALID_PTR)
   return;

  /* Reset counter */

  n = 0;

  /* Loop over depletion intervals */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Add to number of burnup intervals */
      
      if (((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_DEC_TOT) &&
	  ((long)RDB[dep + DEP_HIS_STEP_TYPE] != DEP_STEP_DEC_STEP))
	n++;
      
      /* Next */
      
      dep = NextItem(dep);
    }

  /* Avoid compiler warning */

  norm = -1;

  /* Set decay only mode and check normalization (toi restart file */
  /* testataan että saataisiin tehtyä restartit hajoamislaskuina)  */

  if ((n == 0) && ((long)RDB[DATA_READ_RESTART_FILE] == NO))
    WDB[DATA_BURN_DECAY_CALC] = (double)YES;
  else if ((norm = (long)RDB[DATA_PTR_NORM]) < VALID_PTR)
    Error(0, "Normalization must be defined in burnup mode");

  /* Tää on viritys jolla saadaan matlaboutput.c tulostamaan jos */
  /* restartataan decay-laskuun */

  if ((n == 0) && ((long)RDB[DATA_READ_RESTART_FILE] == YES))
    WDB[DATA_CYCLE_IDX] = RDB[DATA_CRIT_SKIP] + 1.0;

  /* Set pointers to normalization */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Set pointer to normalization */

      WDB[dep + DEP_HIS_PTR_NORM] = (double)norm;

      /* Next normalization */

      if (norm > VALID_PTR)
	if (NextItem(norm) > VALID_PTR)
	  norm = NextItem(norm);
      
      /* Next */

      dep = NextItem(dep);
    }

  /* Check steps */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Get pointer to steps */

      ptr = (long)RDB[dep + DEP_HIS_PTR_STEPS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get total number of steps */
      
      ntot = (long)RDB[dep + DEP_HIS_N_STEPS];

      /* Check step type */

      if (((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_BU_TOT) ||
	  ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DAY_TOT) ||
	  ((long)RDB[dep + DEP_HIS_STEP_TYPE] == DEP_STEP_DEC_TOT))
	{
	  /* Loop over steps and check cumulative */
	  
	  for (n = 1; n < ntot; n++)
	    if (RDB[ptr + n] <= RDB[ptr + n - 1])
	      Error(dep, "Steps not in ascending order after %1.5E", 
		    RDB[ptr + n - 1]);
	}

      /* Check values */
      
      for (n = 0; n < ntot; n++)
	if (RDB[ptr + n] <= 0.0)
	  Error(dep, "Invalid step size %E", RDB[ptr + n]);

      /* Next */

      dep = NextItem(dep);
    }
}

/*****************************************************************************/
