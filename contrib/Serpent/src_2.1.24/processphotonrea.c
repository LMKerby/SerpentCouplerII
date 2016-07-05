/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processphotonrea.c                             */
/*                                                                           */
/* Created:       2011/04/14 (JLe)                                           */
/* Last modified: 2014/11/25 (TKa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes photon reactions                                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessPhotonRea:"


/*****************************************************************************/

void ProcessPhotonRea(long rea) {
  long mt, nuc, loc0;

/*  fprintf(out, "Processing photon data...\n\n"); */

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get reaction mt */
  
  mt = (long)RDB[rea + REACTION_MT];

  /* Pointer to nuclide data */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check that data directory is set */
  if ((long)RDB[DATA_PHOTON_DATA_DIR] <= VALID_PTR)
    Die(FUNCTION_NAME, "Photon data directory not set");


  /* Allocate memory for distribution block */

  loc0 = NewItem(rea + REACTION_PTR_PHOTON_DIST, PHOTON_DIST_BLOCK_SIZE);


  /* Check mt */

  if (mt == 504 || mt == 516 || mt == 522) {
    /* TODO: Tässä luetaan Comptonille (504) , parinmuodostukselle (516) ja
     * valosähköiselle ilmiölle (522) relaksaatio- ja TTB-datat. Eli sama data
     * luetaan kolmeen kertaan. Tämä pitää muuttaa niin, että luetaan vain
     * yhteen kertaan. */


    /***** Read atomic relaxation data ***************************************/

    ProcessRelaxation(loc0, nuc);

    /*************************************************************************/


    /***** TTB data **********************************************************/

    ProcessTTB(loc0, nuc);

    /*************************************************************************/

    }

  if (mt == 502)
    {
    /***** Rayleigh scattering *********************************************/

    ProcessRayleigh(loc0, nuc);

    /***********************************************************************/

    }

  else if (mt == 504)
    {
    /***** Compton scattering **********************************************/

    ProcessCompton(loc0, nuc);

    /***********************************************************************/

    }

  else if (mt == 516)
    {
    /***** Pair production ***************************************************/

    ProcessPairProduction(loc0, nuc);

    /*************************************************************************/

    }

  else if (mt == 522)
    {

    /***** Photoelectric effect***********************************************/

    ProcessPhotoelectric(loc0, nuc);

    /*************************************************************************/

    }
  else
    Die(FUNCTION_NAME, "Invalid reaction mode");
}

/*****************************************************************************/
