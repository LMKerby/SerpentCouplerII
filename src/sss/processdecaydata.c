#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdecaydata.c                             */
/*                                                                           */
/* Created:       2010/09/11 (JLe)                                           */
/* Last modified: 2015/05/21 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes data read from ENDF format decay library and       */
/*              and adds reactions in NUCLIDE block                          */
/*                                                                           */
/* Comments: - Käyttää stability cut-offia suoraan lambdaan. Vanhassa        */
/*             serpentissä cut-off on exp(-lamda*t). Korvaa tää uudella      */
/*             hlfcut-parametrilla joka ottaa argumentiksi puoliintumis-     */
/*             ajan (vuosissa tai vapaissa yksiköissä)                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "radiotox.h"

#define FUNCTION_NAME "ProcessDecayData:"

/*****************************************************************************/

void ProcessDecayData(long nuc)
{
  long ace, ptr, rea, dec, spec, n, loc0, loc1, loc3;
  double T;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "nuc", DATA_ARRAY, nuc);

  /* Get pointer to ace data */

  if ((ace = (long)RDB[nuc + NUCLIDE_PTR_DECAY_ACE]) < 1)
    Die(FUNCTION_NAME, "Nuclide %s has no decay data", 
	GetText(nuc + NUCLIDE_PTR_NAME));

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "ace", ACE_ARRAY, ace);

  /* Check that nuclide and ACE ZAI match */

  if (RDB[nuc + NUCLIDE_ZAI] != ACE[ace + ACE_ZAI])
    Die(FUNCTION_NAME, "Mismatch in ZAI");

  /* Copy remaining parameters for decay type nuclide */

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_DECAY)
    {
      WDB[nuc + NUCLIDE_AWR] = ACE[ace + ACE_AWR];
      WDB[nuc + NUCLIDE_AW] = M_NEUTRON*ACE[ace + ACE_AWR];
    }
  else
    {
      /* Compare ACE and ENDF atomic weight ratios */
    }
  
  /* Reset reaction pointer */

  rea = -1;

  /* Loop over decay modes */

  if ((ptr = (long)ACE[ace + ACE_PTR_DECAY_LIST]) > 0)
    while ((dec = (long)ACE[ptr++]) > 0)
      {
	/* Get half-life for mode */

	T = log(2.0)/(ACE[ace + ACE_LAMBDA]*ACE[dec + DECAY_BR] + ZERO);
	
	/* Compare to cut-off */

	if (T < RDB[DATA_DEP_HALF_LIFE_CUTOFF])
	  {
	    /* Allocate memory for reaction data */

	    rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);

	    /* Put nuclide pointer */

	    WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

	    /* Types (use MT to store first value) */
	    
	    WDB[rea + REACTION_MT] = ACE[dec + DECAY_RTYP1];
	    WDB[rea + REACTION_RTYP2] = ACE[dec + DECAY_RTYP2];
	    WDB[rea + REACTION_RTYP3] = ACE[dec + DECAY_RTYP3];
	    WDB[rea + REACTION_RTYP4] = ACE[dec + DECAY_RTYP4];
	    WDB[rea + REACTION_RTYP5] = ACE[dec + DECAY_RTYP5];
	    
	    /* Target state (ground or isomeric) */
	    
	    WDB[rea + REACTION_RFS] = ACE[dec + DECAY_RFS];
	    
	    /* Q-value */
	    
	    WDB[rea + REACTION_Q] = ACE[dec + DECAY_Q];
	    
	    /* Branching */
	    
	    WDB[rea + REACTION_BR] = ACE[dec + DECAY_BR];
	    
	    /* Set type */

	    WDB[rea + REACTION_TYPE] = (double)REACTION_TYPE_DECAY;
	  }
      }

  /* If pointer is < 0, all decay modes are above cut-off */
  
  if (rea < 0)
    return;

  /* Loop over radiation spectra */

  if ((ptr = (long)ACE[ace + ACE_PTR_RAD_SPEC]) > 0)
    while ((spec = (long)ACE[ptr++]) > 0)
      {
	/* Include only gammas and xrays / annihilation radiation */

	if (((long)ACE[spec + RAD_SPEC_TYPE] == 0) ||
	    ((long)ACE[spec + RAD_SPEC_TYPE] == 9))
	  {
	    /* Pointer to lines */

	    loc0 = (long)ACE[spec + RAD_SPEC_PTR_DISC_E];
	    loc1 = (long)ACE[spec + RAD_SPEC_PTR_DISC_RI];

	    /* Loop over spectral lines */
	    
	    for (n = 0; n < (long)ACE[spec + RAD_SPEC_DISC_NE]; n++)
	      {
		/* Add line */

		loc3 = NewItem(nuc + NUCLIDE_PTR_PHOTON_LINE_SPEC, 
			       PHOTON_LINE_SPEC_BLOCK_SIZE);

		/* Store data */
		
		WDB[loc3 + PHOTON_LINE_SPEC_E] = ACE[loc0++];
		WDB[loc3 + PHOTON_LINE_SPEC_RI] = 
		  ACE[loc1++]*ACE[spec + RAD_SPEC_DISC_NORM];

		/* Add to specific intensity */
		
		WDB[nuc + NUCLIDE_SPEC_PHOTON_I] = 
		  RDB[nuc + NUCLIDE_SPEC_PHOTON_I] +
		  RDB[loc3 + PHOTON_LINE_SPEC_RI]; 
	      }	    
	  }
      }

  /* Check if photon line spectrum was given */

  if ((loc3 = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_LINE_SPEC]) > VALID_PTR)
    {
      /* Close list */

      CloseList(loc3);

      /* Sort by intensity */
      
      SortList(loc3, PHOTON_LINE_SPEC_RI, SORT_MODE_DESCEND);
    }

  /* Copy decay data */

  WDB[nuc + NUCLIDE_LAMBDA] = ACE[ace + ACE_LAMBDA];
  WDB[nuc + NUCLIDE_DECAY_E] = ACE[ace + ACE_DECAY_E];
  WDB[nuc + NUCLIDE_SFBR] = ACE[ace + ACE_SFBR];

  /* Set delayed neutron precursor flag */

  if ((long)ACE[ace + ACE_DELNU_PREC] == YES)
    SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_DELNU_PREC);

  /* Get specific radiotoxicities */

  n = 0;
  while ((long)toxdat[n][0] > 0)
    {
      /* Compare */
      
      if ((long)toxdat[n][0] == (long)ACE[ace + ACE_ZAI])
	break;
      else
	n++;
    }

  /* Check */
  
  if ((long)toxdat[n][0] > 0)
    {
      WDB[nuc + NUCLIDE_SPEC_ING_TOX] = toxdat[n][1];
      WDB[nuc + NUCLIDE_SPEC_INH_TOX] = toxdat[n][2];
    }
  else
    {
      WDB[nuc + NUCLIDE_SPEC_ING_TOX] = 0.0;
      WDB[nuc + NUCLIDE_SPEC_INH_TOX] = 0.0;
    }
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
