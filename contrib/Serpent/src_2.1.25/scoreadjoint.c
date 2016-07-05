/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreadjoint.c                                 */
/*                                                                           */
/* Created:       2011/03/17 (JLe)                                           */
/* Last modified: 2014/11/04 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Adds collision point in history array etc.                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreAdjoint:"

/*****************************************************************************/

void ScoreAdjoint(long part, long mat, long rea, double flx, double x, 
		  double y, double z, double u, double v, double w, double E, 
		  double t, double wgt, double g, long id)
{
#ifdef OLD_HIST

  long hst, ptr;
  double xs;

  /* Pointer to history data */

  if ((hst = (long)RDB[part + PARTICLE_PTR_HIST]) < VALID_PTR)
    return;

  /* Check flux */

  CheckValue(FUNCTION_NAME, "flx", "", flx, ZERO, INFTY);

  if (mat < VALID_PTR)
    return;

  /***************************************************************************/
  
  /***** Score tally *********************************************************/

  /* NOTE: Miksi tässä ei käytetä tota funktion parametrina annettavaa    */
  /*       flx:ää? Noi vaikutusalat pitää myös kertoa g:llä. Mitä hittoa  */
  /*       toi flx edes tekee kun tuolla alempana sille haetaan taas uusi */
  /*       arvo? Ilmeisesti tää on vaan jotain testailua (14.12.2012)     */

  if ((rea = (long)RDB[mat + MATERIAL_PTR_TOTXS]) > VALID_PTR)
    flx = wgt/MacroXS(rea, E, id);

  if ((rea = (long)RDB[mat + MATERIAL_PTR_FISSXS]) > VALID_PTR)
    xs = MacroXS(rea, E, id);
  else
    xs = 0.0;

  /* Check if list is full (this loop is executed only after all positions */
  /* in the array are filled) */

  ptr = (long)RDB[part + PARTICLE_PTR_HIST];
  if ((long)RDB[ptr + HIST_WGT] > 0.0)
    {
      /* Loop over list */

      ptr = hst;
      while ((ptr = PrevItem(ptr)) != hst)
	{
	  /* Break at unassigned point */

	  if ((wgt = RDB[ptr + HIST_WGT]) < 0.0)
	    break;
	 
	  /* Get variables at previous collision point */

	  x = RDB[ptr + HIST_X];
	  y = RDB[ptr + HIST_Y];
	  z = RDB[ptr + HIST_Z];
	  u = RDB[ptr + HIST_U];
	  v = RDB[ptr + HIST_V];
	  w = RDB[ptr + HIST_W];
	  E = RDB[ptr + HIST_E];

	  /* Get flux */

	  flx = RDB[ptr + HIST_FLX];

	  /* Get pointer to reaction and material */

	  rea = (double)RDB[ptr + HIST_PTR_REA];
	  mat = (double)RDB[ptr + HIST_PTR_MAT];

	  ScoreMesh(part, mat, flx*xs, 0.0, x, y, z, E, t, wgt, g, id);
	}
    }

#endif
  
  /***************************************************************************/
}

/*****************************************************************************/
