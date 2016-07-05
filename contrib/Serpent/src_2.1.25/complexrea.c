/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : complexrea.c                                   */
/*                                                                           */
/* Created:       2014/10/16 (JLe)                                           */
/* Last modified: 2015/07/05 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Handles reaction MT 5, which includes several channels and   */
/*              neutron multiplication                                       */
/*                                                                           */
/* Comments: - Weight multiplication handled implicitly                      */
/*           - Contribution to energy deposition may be incorrect            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ComplexRea:"

/*****************************************************************************/

void ComplexRea(long rea, double *E, double *u, double *v, double *w, 
		double wgt1, double *wgt2, double *dE, long id)
{
  long ptr, ne, erg, i;
  double Emin, Emax, f, r;
      
  /* Get pointer to data */

  ptr = (long)RDB[rea + REACTION_PTR_MULT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Get number of energies */

  ne = (long)RDB[ptr++];
  CheckValue(FUNCTION_NAME, "(ne)", "", ne, 2, 10000);

  /* Get pointer to energy grid */
  
  erg = (long)RDB[ptr++];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Check number of energy points */

  if (ne != (long)RDB[erg + ENERGY_GRID_NE])
    Die(FUNCTION_NAME, "Mismatch in number of energy groups");

  /* Get minimum and maximum energy */

  Emin = RDB[erg + ENERGY_GRID_EMIN];
  Emax = RDB[erg + ENERGY_GRID_EMAX];

  if ((r = GridFactor(erg, *E, id)) < 0.0)
   {
      /* Avoid compiler warning */

      i = 0;
      r = 0.0;

      /* Check if energy is above or below limits */

      if (*E > Emax)
	{
	  i = ne - 2;
	  r = 1.0;
	}
      else if (*E < Emin)
	{
	  i = 0;
	  r = 0.0;
	}
      else
	Die(FUNCTION_NAME, "wtf?!?");
    }
  else
    {
      /* Get bin index and adjust factor */

      i = (long)r;
      r = r - (double)i;
    }

  /* Check values */

  CheckValue(FUNCTION_NAME, "i", "", i, 0, ne - 2);
  CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);  

  /* Get multiplication */

  f = RDB[ptr + i]*(1.0 - r) + RDB[ptr + i + 1]*r;
  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 25.0);

  /* Adjust weight */

  if ((*wgt2 = f*wgt1) > 0.0)
    {
      /* Perform inelastic scattering on incident neutron */

      InelasticScattering(rea, E, u, v, w, id);

      /* Calculate change in energy for deposition */

      *dE = *dE - *E*f;
    }
}

/*****************************************************************************/
