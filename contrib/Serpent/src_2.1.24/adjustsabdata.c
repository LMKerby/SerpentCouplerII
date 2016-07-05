/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : adjustsabdata.c                                */
/*                                                                           */
/* Created:       2011/01/23 (JLe)                                           */
/* Last modified: 2015/05/20 (TVi)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: - Adjusts cross sections of nuclides with S(a,b) data        */
/*                                                                           */
/* Comments: - Ton nuklidin uus nimi pitää miettiä (muuta c s:ksi?)          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AdjustSabData:"

/*****************************************************************************/

void AdjustSabData(long nuc)
{
  long rea, pte, ptr, ne, n, i0;
  double Emax, f;

  /* Check S(a,b) flag */

  if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA))
    return;  

  /* Skip for SAB nuclides (frac set to zero -> problems)      */
  /* Nailla nuklideilla ei ole myoskaan elastista sirontaa,    */
  /* vaan ainoastaan S(a,b) -reaktiot joten mitaan ei tarvitse */ 
  /* tehda. (nämä siis .xxt -nuklideita) */

  if((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_SAB)
    return;
  
  /* Reset maximum energy */

  Emax = 0.0;

  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check mt */

      if (((long)RDB[rea + REACTION_MT] == 1002) ||
	  ((long)RDB[rea + REACTION_MT] == 1004))
	{
	  /* Compare maximum energy */

	  if (RDB[rea + REACTION_EMAX] > Emax)
	    Emax = RDB[rea + REACTION_EMAX];
	  
	  /* Number of points */
	  
	  ne = (long)RDB[rea + REACTION_XS_NE];

	  /* Pointer to data */

	  ptr = (long)RDB[rea + REACTION_PTR_XS];

	  /* Get factor */

	  f = RDB[rea + REACTION_SAB_FRAC];

	  /* Check */

	  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
	  
	  /* Loop over points and adjust */

	  for (n = 0; n < ne; n++)
	    WDB[ptr + n] = RDB[ptr + n]*f;
	}

      /* Next reaction */

      rea = NextItem(rea);
    }

  /* In case of OTF S(a,b) interpolation, do not set elastic xs */
  /* to zero. --> skip following */
  
  if( (long)RDB[nuc + NUCLIDE_PTR_SAB] > VALID_PTR)
    return;
    
  /* Pointer to elastic */  
  
  rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
  CheckPointer(FUNCTION_NAME, "(ela)", DATA_ARRAY, rea);
  
  /* Get first energy point and number of points */
  
  i0 = (long)RDB[rea + REACTION_XS_I0];
  ne = (long)RDB[rea + REACTION_XS_NE];
  
  /* Pointer to energy data */
  
  pte = (long)RDB[rea + REACTION_PTR_EGRID];
  pte = (long)RDB[pte + ENERGY_GRID_PTR_DATA];
  
  /* Pointer to data */
    
  ptr = (long)RDB[rea + REACTION_PTR_XS];  
  
  /* Reset cross section */
  
  for (n = 0; n < ne; n++)
    if (RDB[pte + i0 + n] < Emax)
      WDB[ptr + n] = 0.0; 
}

/*****************************************************************************/
