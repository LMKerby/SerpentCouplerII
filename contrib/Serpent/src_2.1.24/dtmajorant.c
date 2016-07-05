/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dtmajorant.c                                   */
/*                                                                           */
/* Created:       2012/10/09 (JLe)                                           */
/* Last modified: 2014/04/04 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Returns delta-tracking majorant for neutrons and photons     */
/*                                                                           */
/* Comments: -Tää on aavistuksen verran sekava systeemi nyt, kun noita       */
/*            majoranttivaikutusalojan on sekä delta-trackingiä että         */
/*            0K dataa varten. Syy miksi tässä ei käytetä MacroXS():ää on    */
/*            se että optimointimoodissa 2 käytetään niitä moniryhmä-        */
/*            vaikutusaloja. Tämän funktion tarkoitus on lähinnä             */
/*            yksinkertaistaa noita pointterihässäköitä Trackin():ssä.       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DTMajorant:"

/*****************************************************************************/

double DTMajorant(long type, double E, long id)
{
  long rea, mat;
  double xs;

  /* Avoid compiler warning */

  rea = -1;
  
  /* Check type and get reaction pointer */

  if (type == PARTICLE_TYPE_NEUTRON)
    rea = (long)RDB[DATA_PTR_MAJORANT];
  else if (type == PARTICLE_TYPE_GAMMA)
    rea = (long)RDB[DATA_PTR_PHOTON_MAJORANT];
  else
    Die(FUNCTION_NAME, "Invalid particle type");

  /* Check Pointer (aikaisemmin on oletettu että st-moodissa   */
  /* majoranttia ei ole laskettu, joten toi voi mennä nulliksi */
  /* siitä syystä) */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get majorant cross section */

  xs = MajorantXS(rea, E, id);

  /* Add alpha cross section */

  if (type == PARTICLE_TYPE_NEUTRON)
    xs = xs + AlphaXS(E);

  /* Add poison xs (tässä voi tulla noi kahteen kertaan jos molemmat on */
  /* sama materiaali) */

  if ((mat = (long)RDB[DATA_MAX_XENON_PTR_MAT]) > VALID_PTR)
    xs = xs + PoisonXS(mat, E, MT_MACRO_TOTXS, id);
  if ((mat = (long)RDB[DATA_MAX_SAMARIUM_PTR_MAT]) > VALID_PTR)
    xs = xs + PoisonXS(mat, E, MT_MACRO_TOTXS, id);

  /* Return cross section */

  return xs;
}

/*****************************************************************************/
