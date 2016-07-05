/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reinittfb.c                                    */
/*                                                                           */
/* Created:       2012/11/04 (JLe)                                           */
/* Last modified: 2012/11/04 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Prepares TFB data structures for the next first time step    */
/*              in a dynamic calculation.                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReinitTFB:"

/*****************************************************************************/

void ReinitTFB()
{
  long tfb, reg;

  /* Loop over temperature feedbacks */

  tfb = (long)RDB[DATA_PTR_TFB0];
  while(tfb > VALID_PTR)
    {
      /* Loop over regions */

      reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
      while (reg > VALID_PTR)
	{
	  /* Set inital guesses */
	  
	  WDB[reg + TFB_REG_ITER_R0] = RDB[reg + TFB_REG_INI_R0];
	  WDB[reg + TFB_REG_ITER_R1] = RDB[reg + TFB_REG_INI_R1];
	  WDB[reg + TFB_REG_ITER_C0] = RDB[reg + TFB_REG_INI_C0];
	  WDB[reg + TFB_REG_ITER_C1] = RDB[reg + TFB_REG_INI_C1];
	  WDB[reg + TFB_REG_ITER_C2] = RDB[reg + TFB_REG_INI_C2];

	  /* Next region */

	  reg = NextItem(reg);
	}

      /* Next feedback */
      
      tfb = NextItem(tfb);
    }
}

/*****************************************************************************/
