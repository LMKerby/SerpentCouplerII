/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processbradata.c                               */
/*                                                                           */
/* Created:       2013/09/03 (JLe)                                           */
/* Last modified: 2014/05/03 (JLe)                                           */
/* Version:       2.1.21                                                     */
/*                                                                           */
/* Description: Reconstructs energy-dependent isomeric branching ratios      */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessBraData:"

/*****************************************************************************/

void ProcessBraData(long nuc)
{
  long rea, erg, ne, np, loc0, loc1, ptr, n, INTT, i0;
  long  dum1, dum2;
  double *f;

  /* Check pointer to data */

  if ((long)RDB[nuc + NUCLIDE_BRA_TYPE] != BRA_TYPE_ENE)
    return;

#ifdef DEBUG

  Warn(FUNCTION_NAME, "Lineaariset interpoloinnit kaikissa");

#endif

  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while(rea > VALID_PTR)
    {
      /* Check pointer to state */

      if ((loc1 = (long)RDB[rea + REACTION_PTR_BRA_STATE]) < VALID_PTR)
	{
	  /* Next reaction */

	  rea = NextItem(rea);

	  /* Cycle loop */

	  continue;
	}

#ifdef DEBUG

      printf("%s %ld %ld\n", GetText(nuc + NUCLIDE_PTR_NAME),
	     (long)RDB[rea + REACTION_MT],
	     (long)RDB[rea + REACTION_RFS]);

#endif
      
      /* Get pointer to energy grid */

      erg = (long)RDB[rea + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get number of points */

      ne = (long)RDB[erg + ENERGY_GRID_NE];

      /* Get pointer to data */

      erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Allocate memory for branching ratios */

      f = (double *)Mem(MEM_ALLOC, ne, sizeof(double));
      
      /* Pointer to branching ratio energy grid */

      loc0 = (long)RDB[loc1 + BRA_STATE_PTR_ERG];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Get number of points */

      np = (long)RDB[loc1 + BRA_STATE_NP];      

      /* Get interpolation law */

      INTT = (long)RDB[loc1 + BRA_STATE_INTT];

      /* Pointer to branching ratios */
      
      loc1 = (long)RDB[loc1 + BRA_STATE_PTR_FRAC];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Check data */

      for (n = 0; n < np; n++)
	if ((RDB[loc1 + n] < 0.0) || (RDB[loc1 + n] > 1.0))
	  Die(FUNCTION_NAME, "Error in data (1)");

      /* Interpolate (NOTE: lin-lin) */

      n = InterpolateData(&RDB[erg], f, ne, &RDB[loc0],
			  &RDB[loc1], np, 0, &dum1, &dum2); 

      /* Check data */

      for (n = 0; n < ne; n++)
	if ((f[n] < 0.0) || (f[n] > 1.0))
	  Die(FUNCTION_NAME, "Error in data (2)");

      /* Get pointer to cross section data */

      ptr = (long)RDB[rea + REACTION_PTR_XS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get number of points and index to first point */

      ne = (long)RDB[rea + REACTION_XS_NE];
      i0 = (long)RDB[rea + REACTION_XS_I0];

      /* Multiply */

      for (n = 0; n < ne; n++)
	WDB[ptr + n] = WDB[ptr + n]*f[n + i0];

      /* Free allocated memory */

      Mem(MEM_FREE, f);

      /* Next reaction */

      rea = NextItem(rea);
    }
}

/*****************************************************************************/
