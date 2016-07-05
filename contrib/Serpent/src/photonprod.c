/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : photonprod.c                                   */
/*                                                                           */
/* Created:       2016/01/30 (JLe)                                           */
/* Last modified: 2016/04/13 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Produce photons from neutron collisions                      */
/*                                                                           */
/* Comments: - Ures sampling not accounted for                               */
/*                                                                           */
/*           - Probably does not work with implicit capture                  */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PhotonProd:"
 
/*****************************************************************************/

void PhotonProd(long nuc, double x, double y, double z, double u0, double v0, 
		double w0, double E0, double wgt, double t, long id)
{
  long mode, loc0, loc1, ptr, ntot, n, m, det, nb;
  double elaxs, totxs, inlxs, prodxs, yld, xs, mu, E;

  /* Check nuclied pointer */

  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check mode */

  if ((mode = (long)RDB[DATA_PHOTON_PRODUCTION]) == NO)
    return;

  /* Pointer to total */

  if ((loc0 = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_PROD]) < VALID_PTR)
    return;

  /* Get nuclide total and elastic cross sections */
  /* (use UresDiluMicroXS() instead?) */

  ptr = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  totxs = MicroXS(ptr, E0, id);

  /* Check total xs */

  if (totxs == 0.0)
    return;

  /* Get nuclide photon production xs */

  ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTPRODXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  prodxs = MicroXS(ptr, E0, id);

  /* Check */

  if (prodxs == 0.0)
    return;

  /* Avoid compiler warning */

  ntot = -1;

  /* Check mode */

  if (mode == PHOTON_PROD_ANA)
    {  
      /* Analog mode, get elastic cross section */
      
      ptr = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      elaxs = MicroXS(ptr, E0, id);

      /* Reject elastic reactions */
	  
      if (RandF(id) < elaxs/totxs)
	return;

      /* Calculate inelastic cross section */
      
      inlxs = totxs - elaxs;
      CheckValue(FUNCTION_NAME, "inlxs", "", inlxs, ZERO, INFTY);  
      
      /* Calculate yield */
      
      yld = prodxs/inlxs;
      
      /* Sample number of emitted photons */
      
      ntot = (long)yld;
      
      if (RandF(id) < yld - (double)ntot)
	ntot++;
      
      /* Check zero */
      
      if (ntot == 0)
	return;
    }
  else if (mode == PHOTON_PROD_IMP)
    {
      /* Implicit mode, get number of photons */

      ntot = 10;

      /* Calculate weight */

      wgt = wgt*prodxs/totxs/((double)ntot);
    }

  /* Get sum of partial photon production cross sections */
  /* (this is the first item in the list) */

  ptr = (long)RDB[loc0 + PHOTON_PROD_PTR_PRODXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  prodxs = MicroXS(ptr, E0, id);

  /* Loop over photons */

  for (n = 0; n < ntot; n++)
    {
      /* Re-sampling loop (12024.03c default-kirjastossa feilaa) */

      for (m = 0; m < 10; m++)
	{
	  /* Sample fraction of total production xs */
	  
	  xs = RandF(id)*prodxs;
	  
	  /* Loop over reactions */

	  loc1 = NextItem(loc0);
	  while (loc1 > VALID_PTR)
	    {
	      /* Subtract production xs */
	      
	      ptr = (long)RDB[loc1 + PHOTON_PROD_PTR_PRODXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      xs = xs - MicroXS(ptr, E0, id);
	      
	      /* Check */
	      
	      if (xs <= 0.0)
		break;
	      
	      /* Next reaction */
	      
	      loc1 = NextItem(loc1);
	    }
	  
	  /* Check sample */

	  if (loc1 > VALID_PTR)
	    break;	  
	}

      /* Check pointer */

      if (m == 10)
	Die(FUNCTION_NAME, "Failed to sample reaction (E = %E, xs = %E)", 
	    E0, xs);

      /* Sample direction */

      ptr = (long)RDB[loc1 + PHOTON_PROD_PTR_ANG];
      mu = SampleMu(-1, ptr, E0, id);

      /* Sample energy */

      ptr = (long)RDB[loc1 + PHOTON_PROD_PTR_ERG];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      SampleENDFLaw(-1, ptr, E0, &E, NULL, id);

      /***********************************************************************/

      /***** Store data to first detector for testing ************************/

      /* Pointer to detectors */
      
      det = (long)RDB[DATA_PTR_DET0];
      CheckPointer(FUNCTION_NAME, "(det)", DATA_ARRAY, det);

      /* Energy cut-off */

      if (E < 1E-3)
	continue;

      /* Get bin */

      nb = DetBin(det, -1, x, y, z, E, t, id);

      /* Get pointer to statistics */

      ptr = (long)RDB[det + DET_PTR_STAT];
      CheckPointer(FUNCTION_NAME, "(stat)", DATA_ARRAY, ptr);

      /* Add to data */

      if (nb > -1)
	AddBuf(1.0, wgt, ptr, id, -1, nb, 0);

      /* Pointer to second detector */

      ptr = NextItem(ptr);
      CheckPointer(FUNCTION_NAME, "(stat)", DATA_ARRAY, ptr);

      /* Bin by cosine */

      nb = (long)(200.0*(mu + 1.0)*0.5);
      AddBuf(1.0, wgt, ptr, id, -1, nb, 0);      

      /***********************************************************************/
    }
}

/*****************************************************************************/
