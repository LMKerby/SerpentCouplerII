#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : coefcycle.c                                    */
/*                                                                           */
/* Created:       2014/04/15 (JLe)                                           */
/* Last modified: 2014/11/15 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Cycle over burnup points in coefficient calculations         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CoefCycle:"

/*****************************************************************************/

void CoefCycle()
{
  long loc0, loc1, n;
  char tmpstr[MAX_STR];

  /* Get pointer to coefficient matrix */

  loc0 = (long)RDB[DATA_PTR_COEF0];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Set number of burnup points */
  
  WDB[DATA_TOT_COEF_BU] = RDB[loc0 + COEF_N_BU];

  /* Pointer to burnup points */

  loc1 = (long)RDB[loc0 + COEF_PTR_BU_PTS];
  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

  /* Remove existing file */

  if ((long)RDB[DATA_COEF_CALC_IDX] == 1)
    {
      sprintf(tmpstr, "%s.coe", GetText(DATA_PTR_INPUT_FNAME));
      remove(tmpstr);
    }

  /* Loop over points */

  for (n = 0; n < (long)RDB[loc0 + COEF_N_BU]; n++)
    {
      /* Reset some timers */

      ResetTimer(TIMER_BURNUP);
      ResetTimer(TIMER_BATEMAN);

      /* Put index */

      WDB[DATA_COEF_CALC_BU_IDX] = (double)(n + 1);

      /* Add to current run index */

      WDB[DATA_COEF_CALC_RUN_IDX] = RDB[DATA_COEF_CALC_RUN_IDX] + 1;
	
      /* Put restart point */

      WDB[DATA_PTR_COEF_BU_PT] = RDB[loc1];
      if ((WDB[DATA_RESTART_READ_POINT] = atof(GetText(loc1++))) != 0.0)
	{
	  /* Set flag */
	  
	  WDB[DATA_READ_RESTART_FILE] = YES;
	  
	  /* Read restart file */
	  
	  ReadRestartFile(NO);
	}

      /* Prepare transport cycle */

      PrepareTransportCycle();

      /* Run transportcycle(s) */

      do
	{

	  /* Prepare coupled calculation iteration if needed */

	  PrepareCCIter();

	  /* Run transport cycle */

	  TransportCycle();

	  /* Iterate coupled calculation routines */

	  IterateCC();

	  /* Repeat if needed */
	}
      while(RDB[DATA_ITERATE] == (double)YES);

      /* Signal externally coupled program to end calculation */

      SignalExternal(SIGUSR2);

      /* Print coefficient output */

      CoefOutput();
    }
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
