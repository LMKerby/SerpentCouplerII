/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : collecttfb.c                                   */
/*                                                                           */
/* Created:       2012/10/29 (JLe)                                           */
/* Last modified: 2012/11/03 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Calculates linear powers for TFB                             */
/*                                                                           */
/* Comments: - Nää vois ehkä kuitenkin tehdä jossain muualla?                */
/*           - Ton tehon tallentaminen AddStat():lla on vähän huono tapa     */
/*             hoitaa homma, mutta tässä yksinkertaisin kun se pitää jakaa   */
/*             aikavälin pituudella (joka saattaa muuttua)                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CollectTfb:"

/*****************************************************************************/

void CollectTfb()
{
  long tfb, nst, reg, ptr, n, i, tb;
  double norm, val, dt;

  /* Reduce scoring buffer */

  ReduceBuffer();

  /* Get normalization coefficient */

  norm = NormCoef(PARTICLE_TYPE_NEUTRON);

  /* Get time bin index */

  tb = (long)RDB[DATA_DYN_TB];

  /* Get time interval */

  if (RDB[DATA_TIME_CUT_TMAX] < INFTY)
    dt = RDB[DATA_TIME_CUT_TMAX] - RDB[DATA_TIME_CUT_TMIN];	  
  else
    dt = 1.0;

  /* Loop over feedbacks */
      
  tfb = (long)RDB[DATA_PTR_TFB0];
  while (tfb > VALID_PTR)
    {
      /* Pointer to nest */

      nst = (long)RDB[tfb + TFB_PTR_NST];
      CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

      /* Get the number of nests */

      n = (long)RDB[nst + NEST_COUNT];

      /* Reset total power*/

      WDB[tfb + TFB_ITER_POW] = 0.0;

      /* Reset index */
      
      i = 0;
      
      /* Loop over regions */
      
      reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
      while (reg > VALID_PTR)
	{
	  /* Power */
	  
	  ptr = (long)RDB[tfb + TFB_PTR_MEAN_POW];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, i, tb)/dt;
	  AddStat(val*norm, ptr, i, tb);
	  
	  /* Set power for next iteration */
	  
	  WDB[reg + TFB_REG_ITER_POW] = val*norm/((double)n);

	  /* The power flowing to this region from inner regions */

	  WDB[reg + TFB_REG_ITER_POWIN] = RDB[tfb + TFB_ITER_POW];
	  
	  /* Calculate contribution to total power */
	  
          WDB[tfb + TFB_ITER_POW] = RDB[tfb + TFB_ITER_POW] 
	    + val*norm/((double)n);

	  /* Update index */
	  
	  i++;
	  
	  /* Next */
	  
	  reg = NextItem(reg);
	}
      
      /* Next */
      
      tfb = NextItem(tfb);
    }    
}

/*****************************************************************************/
