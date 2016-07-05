/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scorecapture.c                                 */
/*                                                                           */
/* Created:       2011/07/30 (JLe)                                           */
/* Last modified: 2014/11/04 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Scores capture parameters                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreCapture:"

/*****************************************************************************/

void ScoreCapture(long mat, long rea, double wgt, long id)
{
  long ptr, nuc, ZAI;

  /* Check reaction pointer (called without pointer to handle implicit */
  /* capture without implicit reaction rates) */

  if (rea > VALID_PTR)
    {
      /* Get nuclide pointer */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
      
      /* Get nuclide ZAI */
      
      ZAI = (long)RDB[nuc + NUCLIDE_ZAI];
      
      /* Check ZAI and MT */
      
      if (((ZAI == 902320) || (ZAI == 922380)) && 
	  ((long)RDB[rea + REACTION_MT] == 102))
	{
	  /* Score analog fertile capture rate */
	  
	  ptr = (long)RDB[RES_ANA_CONV_RATIO];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(1.0, wgt, ptr, id, 1);
	}
      
      /* Check neutron number and MT */
      
      if ((((long)RDB[nuc + NUCLIDE_A] - (long)RDB[nuc + NUCLIDE_Z]) % 2) &&
	  ((long)RDB[rea + REACTION_MT] == 102))
	{
	  /* Score analog fissile fission rate */
	  
	  ptr = (long)RDB[RES_ANA_CONV_RATIO];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(1.0, wgt, ptr, id, 2);
	}

      /* Total and nuclide-wise rates */

      ptr = (long)RDB[RES_ANA_CAPT_FRAC];
      CheckPointer(FUNCTION_NAME, "(ptr 4)", DATA_ARRAY, ptr);
      AddBuf(1.0, wgt, ptr, id, -1, 0, 0);
      
      if ((long)RDB[nuc + NUCLIDE_ZAI] == 902320)
	AddBuf(1.0, wgt, ptr, id, -1, 1, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 922330)
	AddBuf(1.0, wgt, ptr, id, -1, 2, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 922350)
	AddBuf(1.0, wgt, ptr, id, -1, 3, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 922380)
	AddBuf(1.0, wgt, ptr, id, -1, 4, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 942390)
	AddBuf(1.0, wgt, ptr, id, -1, 5, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 942400)
	AddBuf(1.0, wgt, ptr, id, -1, 6, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 942410)
	AddBuf(1.0, wgt, ptr, id, -1, 7, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 541350)
	AddBuf(1.0, wgt, ptr, id, -1, 8, 0);
      else if ((long)RDB[nuc + NUCLIDE_ZAI] == 621490)
	AddBuf(1.0, wgt, ptr, id, -1, 9, 0);
    }

  /* Score analog reaction rate if not implicit mode */

  if ((long)RDB[DATA_OPTI_IMPLICIT_RR] == NO)
    {
      /* Total reaction rate */

      ptr = (long)RDB[RES_TOT_NEUTRON_RR];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, wgt, ptr, id, 0);

      /* Total capture rate */

      ptr = (long)RDB[RES_TOT_CAPTRATE];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(1.0, wgt, ptr, id, 0);
      
      /* Check material fissile flag */
  
      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
	AddBuf1D(1.0, wgt, ptr, id, 1);
      else
	AddBuf1D(1.0, wgt, ptr, id, 2);
    }
}

/*****************************************************************************/
