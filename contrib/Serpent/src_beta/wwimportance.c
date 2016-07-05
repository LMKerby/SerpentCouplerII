/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : wwimportance.c                                 */
/*                                                                           */
/* Created:       2015/10/08 (JLe)                                           */
/* Last modified: 2015/10/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Return the importance for weight windows                     */
/*                                                                           */
/* Comments: - Tässä luupataan kaikkien yli ja palautetaan kumulatiivinen    */
/*             arvo.                                                         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WWImportance:"

/*****************************************************************************/

double WWImportance(double x, double y, double z, double E)
{
  long wwd, msh, erg, ptr;
  double f, p;

  /* Pointer to weight window structure */
  
  wwd = (long)RDB[DATA_PTR_WWD0];
  CheckPointer(FUNCTION_NAME, "(wwd)", DATA_ARRAY, wwd);

  /* Loop over weight window structures */

  while (wwd > VALID_PTR)
    {
      /* Pointer to energies */ 

      if ((erg = (long)RDB[wwd + WWD_PTR_ERG]) > VALID_PTR)
	{
	  /* Loop over energies and find interval */

	  while (erg > VALID_PTR)
	    {
	      /* Check interval */

	      if ((E > RDB[erg + WWD_ERG_EMIN]) && 
		  (E < RDB[erg + WWD_ERG_EMAX]))
		{
		  /* Get importance */
		  
		  f = RDB[erg + WWD_ERG_IMP];

		  /* Get exponential */

		  if ((p = RDB[wwd + WWD_POW]) > 0.0)
		    f = powl(f, p);
		  
		  /* Return value */

		  return f;
		}
	      
	      /* Next */

	      erg = NextItem(erg);
	    }
	}

      /* Pointer to mesh */
      
      if ((msh = (long)RDB[wwd + WWD_PTR_MESH]) < VALID_PTR)
	{
	  /* Pointer to next */

	  wwd = NextItem(wwd);

	  /* Cycle loop */

	  continue;
	}

      /* Pointer to data */
      
      if ((ptr = MeshPtr(msh, x, y, z)) > VALID_PTR)
	{
	  /* Pointer to structure */

	  ptr = (long)RDB[ptr];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Get importance */
	  
	  f = RDB[ptr + WWD_MESH_IMP];
	  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);

	  /* Get exponential */

	  if ((p = RDB[wwd + WWD_POW]) > 0.0)
	    f = powl(f, p);
      
	  /* Check normalization and scale */

	  if (RDB[wwd + WWD_NORM_FACT] > 0.0)
	    {
	      f = f*RDB[wwd + WWD_NORM_FACT];
	      CheckValue(FUNCTION_NAME, "f", "", f, 0.0, INFTY);
	    }

	  /* Check */
	  
	  if (f > 0.0)
	    return f;
	}

      /* Next */

      wwd = NextItem(wwd);
    }

  /* Return zero */

  return 0.0;
}

/*****************************************************************************/
