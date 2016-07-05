/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sabscattering.c                                */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2015/06/16 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Handles S(a,b) scattering laws                               */
/*                                                                           */
/* Comments: - Interpolation factorit ja interpoloinnin termit on eri päin   */
/*             kuin Serpent 1:ssä. Pitäisi tuottaa sama tulos (19.1.2012).   */
/*                                                                           */
/*           - Elastic -moodin tyyppiä 3 ei ole testattu (kaikki käyttää     */
/*             exact treatmentia?), mutta sen pitäisi mennä samalla tavalla  */
/*             kuin inelastic-moodin kulman arvonta.                         */
/*                                                                           */
/*           - Tästä on kaksi versiota, toi uudempi on koodattu 20.1.2011,   */
/*             vanhempi on lähempänä Serpent 1:n menetelmiä.                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SabScattering:"

/* Select which version to use */

#ifndef SABSCATTERING_USE_OLD_VERSION

/*****************************************************************************/

void SabScattering(long rea, double *E, double *u, double *v, double *w, 
		   long id)
{
  double mu, E0, r, d1, d2, a;
  long erg, law, l0, l1, l2, nc2, ctype, ne, ne2, n, i, k, l;

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get pointer to energy distribution */

  erg = (long)RDB[rea + REACTION_PTR_ERG];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Check initial energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Initial energy */

  E0 = *E;

  /* Avoid compiler warning */

  mu = 1.0;
  
  /* Get pointer to data */
  
  l0 = (long)RDB[erg + ERG_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

  /* Get scattering mode (elastic or inelastic) */

  law = (long)RDB[erg + ERG_LAW];

  /* Get cosine distribution type */

  ctype = (long)RDB[l0++];
  
  /* Get pointer to incident energy grid (recycle pointer) */
  
  erg = (long)RDB[l0++];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Sample energy and/or scattering cosine */

  if (ctype == 4)
    {
      /***********************************************************************/

      /***** Exact treatment of elastic scattering ***************************/

      /* Check scattering mode */

      if (law != 1002)
	Die(FUNCTION_NAME, "Invalid scattering mode");

      /* Get index to energy grid */

      if ((i = GridSearch(erg, E0)) < 0)
	{
	  /* Energy below or above grid, get grid size */

	  ne = (long)RDB[erg + ENERGY_GRID_NE];		
	  
	  /* Avoid compiler warning */
	  
	  i = -1;
	  
	  /* Check incident energy */
	  
	  if (E0 < RDB[erg + ENERGY_GRID_EMIN])
	    i = 0;
	  else if (E0 > RDB[erg + ENERGY_GRID_EMAX])
	    i = ne - 2;
	  else
	    Die(FUNCTION_NAME, "law 1002: i = %E", i);
	}
      
      /* Skip size */

      l0++;

      /* Pointers to data */
      
      l1 = l0;
      l2 = l0 + i + 1;
      
      /* Sample */
      
      d1 = RandF(id)*RDB[l2];

      /* Search */
      
      while (l2 != l1 + 1)
	{	    
	  n = (long)((l1 + l2)/2.0);
	  
	  if (d1 < RDB[n])
	    l2 = n;
	  else
	    l1 = n;		
	}
      
      /* Pointer to incident energy grid data (recycle pointer) */
      
      erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
      
      /* Pointer to value */

      erg = erg + l2 - l0 - 1;

      /* Calculate mu */
      
      mu = 1.0 - 2.0*RDB[erg]/E0;

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Other elastic or inelastic modes ********************************/
      
      /* Get interpolation factor */
      
      if ((r = GridFactor(erg, E0, id)) < 0.0)
	{
	  /* Get number of energies */
	  
	  ne = (long)RDB[erg + ENERGY_GRID_NE];		
	  i = -1;
	  
	  /* Check energy */
	  
	  if (E0 < RDB[erg + ENERGY_GRID_EMIN])
	    {
	      i = 0;
	      r = 0.0;
	    }
	  else if (E0 > RDB[erg + ENERGY_GRID_EMAX])
	    {
	      i = ne - 2;
	      r = 1.0;
	    }
	  else
	    Die(FUNCTION_NAME, "wtf?");
	}
      else
	{
	  /* Separate integer and decimal parts */
	  
	  i = (long)r;
	  r = r - (double)i;
	}
      
      /* Check interpolation factor */
      
      CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);
      
      /* Avoid compiler warning */

      nc2 = -1;
      l1 = -1;
      l2 = -1;

      /* Check scattering mode */

      if (law == 1004)
	{
	  /* Inelastic mode, get number of secondary energies and cosines */
      
	  ne2 = (long)RDB[l0++];
	  nc2 = (long)RDB[l0++];
      
	  /* Get secondary energy type */
	  
	  n = (long)RDB[l0++];
	  
	  /* Get index to energy bin (from sabcol-subroutine in MCNP source) */
	  
	  if (n == 0)
	    k = (long)(RandF(id)*ne2);
	  else if ((a = RandF(id)*(ne2 - 3.0)) > 1.0)
	    k = (long)a + 1;
	  else if (a > 0.6)
	    k = ne2 - 2;
	  else if (a > 0.5)
	    k = ne2 - 1;
	  else if (a > 0.1)
	    k = 1;
	  else
	    k = 0;
	  
	  /* Check value */
	  
	  CheckValue(FUNCTION_NAME, "k", " (law 1004)", k, 0, ne2 - 1);
	  
	  /* Pointers to distribution */
	  
	  l1 = (long)RDB[l0 + i] + k*(nc2 + 1);
	  l2 = (long)RDB[l0 + i + 1] + k*(nc2 + 1);
      
	  /* Get energy values */
      
	  d1 = RDB[l1];
	  d2 = RDB[l2];
	  
	  /* Interpolate energy value */
	  
	  *E = d1 + r*(d2 - d1);

	  /* Update pointers (cosines are listed after each energy) */

	  l1++;
	  l2++;
	}
      else if (law == 1002)
	{
	  /* NOTE: tää on kohtuullisen harvinainen sirontalaki, johon */
	  /* törmää esim. h/zr ja zr/h -kirjastoissa. Ei ole kunnolla */
	  /* testattu. */

	  /* Elastic, get number of cosines */
      
	  nc2 = (long)RDB[l0++];

	  /* Get pointers */
	  
	  l1 = l0 + i*nc2;
	  l2 = l0 + (i + 1)*nc2;
	}
      else
	Die(FUNCTION_NAME, "Invalid scattering mode");

      /* Check cosine distribution type */

      if (ctype == 3)
	{
	  /* Sample from discrete cosines */
	  
	  l = (long)(((double)nc2)*RandF(id));
	  
	  /* Get cosine values */
	  
	  d1 = RDB[l1 + l];
	  d2 = RDB[l2 + l];

	  /* Check values */
	  
	  CheckValue(FUNCTION_NAME, "d1 (law 1004)", "", d1, -1.0, 1.0);
	  CheckValue(FUNCTION_NAME, "d2 (law 1004)", "", d2, -1.0, 1.0);
	  
	  /* Interpolate mu value */
	  
	  mu = d1 + r*(d2 - d1);
	}
      else
	Die(FUNCTION_NAME, "Invalid or unsupported cosine distribution mode");

      /************************************************************************/
    }

  /* Check energy and cosines */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "r", "", *u**u+*v**v+*w**w - 1.0, -1E-5, 1E-5);

  /* Rotate direction cosines around a random azimuthal angle */
      
  AziRot(mu, u, v, w, id);      
}

/*****************************************************************************/

#else

/*****************************************************************************/

void SabScattering(long rea, double *E, double *u, double *v, double *w, 
		   long id)
{
  double mu, E0, u0, v0, w0, r, d1, d2, a;
  long erg, law, ptr, l0, l1, l2, nc2, mode, ne, ne2, n, i, k, l, type;

  /***************************************************************************/

  /***** Get initial values and check ****************************************/

  /* Check reaction pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Get pointer to energy distribution */

  erg = (long)RDB[rea + REACTION_PTR_ERG];
  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

  /* Check initial energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Get energy distribution type */

  law = (long)RDB[erg + ERG_LAW];

  /***************************************************************************/

  /***** Remember some values before the collision ***************************/

  /* Initial energy */

  E0 = *E;

  /* Initial direction cosines */

  u0 = *u;
  v0 = *v;
  w0 = *w;

  /* Avoid compiler warning */

  mu = 1.0;

  /* Get pointer to data */
  
  l0 = (long)RDB[erg + ERG_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(l0)", DATA_ARRAY, l0);

  /* Get scattering mode */

  mode = (long)RDB[l0++];
  
  /* Get pointer to incident energy grid */
  
  ptr = (long)RDB[l0++];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /***************************************************************************/

  /***** Sample energy and scattering angle **********************************/

  if (law == 1002)
    {
      /***********************************************************************/

      /***** S(a,b) elastic scattering ***************************************/

      /* Get number of cosines */
      
      nc2 = (long)RDB[l0++];

      /* Check mode */

      if (mode == 4)
	{
	  /* Exact treatment, get index to energy grid */

	  if ((i = GridSearch(ptr, E0)) < 0)
	    {
	      /* Energy below or above grid, get grid size */

	      ne = (long)RDB[ptr + ENERGY_GRID_NE];		

	      /* Avoid compiler warning */
	      
	      i = -1;
	  
	      /* Check incident energy */
	  
	      if (E0 < RDB[ptr + ENERGY_GRID_EMIN])
		i = 0;
	      else if (E0 > RDB[ptr + ENERGY_GRID_EMAX])
		i = ne - 2;
	      else
		Die(FUNCTION_NAME, "law 1002: i = %E", i);
	    }

	  /* Pointers to data */

	  l1 = l0;
	  l2 = l0 + i + 1;

	  /* Sample */

	  d1 = RandF(id)*RDB[l2];

	  /* Search */

	  while (l2 != l1 + 1)
	    {	    
	      n = (long)((l1 + l2)/2.0);
	      
	      if (d1 < RDB[n])
		l2 = n;
	      else
		l1 = n;		
	    }

	  /* Pointer to incident energy grid data */

	  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Pointer to value */

	  ptr = ptr + l2 - l0 - 1;

	  /* Calculate mu */

	  mu = 1.0 - 2.0*RDB[ptr]/E0;
	}
      else if (mode == 3)
	{
	  /* Discrete cosines */
	  
	  Die (FUNCTION_NAME, "Discrete cosines (not tested yet?)");
      
	  /* Find grid interval and calculate factor */

	  if ((i = GridSearch(ptr, E0)) < 0)
	    {
	      /* Energy below or above grid, get grid size */
	      
	      ne = (long)RDB[ptr + ENERGY_GRID_NE];		
	      
	      /* Avoid compiler warning */
	      
	      i = -1;
	      r = -1.0;
	      
	      /* Check incident energy */
	      
	      if (E0 < RDB[ptr + ENERGY_GRID_EMIN])
		{
		  i = 0;
		  r = 1.0;
		}
	      else if (E0 > RDB[ptr + ENERGY_GRID_EMAX])
		{
		  i = ne - 2;
		  r = 0.0;
		}
	      else
		Die(FUNCTION_NAME, "wtf?");
	    }
	  else
	    {
	      /* Interpolation factor */
	      
	      r = (RDB[l0 + i + 1] - E0)/(RDB[l0 + i + 1] - RDB[l0 + i]);
	      
	      /* Check value */
	      
	      CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);
	    }
	    
	  /* Get pointers */
	  
	  l1 = l0 + i*nc2;
	  l2 = l0 + (i + 1)*nc2;
	  
	  /* Sample bin */
	  
	  l = (long)(((double)nc2)*RandF(id));
	  
	  /* Get cosine values */
	  
	  d1 = RDB[l1 + l];
	  d2 = RDB[l2 + l];
	  
	  /* Check values */
	  
	  CheckValue(FUNCTION_NAME, "d1 (law 1002)", "", d1, -1.0, 1.0);
	  CheckValue(FUNCTION_NAME, "d2 (law 1002)", "", d2, -1.0, 1.0);
	  
	  /* Interpolate mu value */
	  
	  mu = d2 + r*(d1 - d2);
	}
      else
	Die(FUNCTION_NAME, "Invalid mode");
      
      /***********************************************************************/
    }
  else if (law == 1004)
    {
      /***** S(a,b) inelastic scattering *************************************/
      
      /* Get interpolation factor */
      
      if ((r = GridFactor(ptr, E0, id)) < 0.0)
	{
	  /* Get number of energies */
	    
	  ne = (long)RDB[ptr + ENERGY_GRID_NE];		
	  i = -1;
	  
	  /* Check energy */
	  
	  if (E0 < RDB[ptr + ENERGY_GRID_EMIN])
	    {
	      i = 0;
	      r = 0.0;
	    }
	  else if (E0 > RDB[ptr + ENERGY_GRID_EMAX])
	    {
	      i = ne - 2;
	      r = 1.0;
	    }
	  else
	    Die(FUNCTION_NAME, "wtf?");
	}
      else
	{
	  /* Separate integer and decimal parts */
	  
	  i = (long)r;
	  r = r - (double)i;
	}
      
      /* Check interpolation factor */
      
      CheckValue(FUNCTION_NAME, "r (law 1004)", "", r, 0.0, 1.0);
      
      /* Get number of secondary energies and cosines */
      
      ne2 = (long)RDB[l0++];
      nc2 = (long)RDB[l0++];
      
      /* Get secondary energy type and inelastic scattering mode */
      
      type = (long)RDB[l0++];
      
      /* Get index to energy bin (from sabcol-subroutine in MCNP source) */
      
      if (type == 0)
	k = (long)(RandF(id)*ne2);
      else if ((a = RandF(id)*(ne2 - 3.0)) > 1.0)
	k = (long)a + 1;
      else if (a > 0.6)
	k = ne2 - 2;
      else if (a > 0.5)
	k = ne2 - 1;
      else if (a > 0.1)
	k = 1;
      else
	k = 0;
      
      /* Check value */
      
      CheckValue(FUNCTION_NAME, "k", " (law 1004)", k, 0, ne2 - 1);
      
      /* Pointers to distribution */
      
      l1 = (long)RDB[l0 + i] + k*(nc2 + 1);
      l2 = (long)RDB[l0 + i + 1] + k*(nc2 + 1);
      
      /* Get energy values */
      
      d1 = RDB[l1];
      d2 = RDB[l2];
      
      /* Interpolate energy value */
      
      *E = d1 + r*(d2 - d1);
      
      /* Check cosine distribution mode */

      if (mode == 3)
	{
	  /* Sample from discrete cosines */
	  
	  l = (long)(((double)nc2)*RandF(id));
	  
	  /* Get cosine values */
	  
	  d1 = RDB[l1 + 1 + l];
	  d2 = RDB[l2 + 1 + l];
	  
	  /* Check values */
	  
	  CheckValue(FUNCTION_NAME, "d1 (law 1004)", "", d1, -1.0, 1.0);
	  CheckValue(FUNCTION_NAME, "d2 (law 1004)", "", d2, -1.0, 1.0);
	  
	  /* Interpolate mu value */
	  
	  mu = d1 + r*(d2 - d1);
	}
      else
	Die(FUNCTION_NAME, "Invalid inelastic scattering mode");
      
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid S(a,b) mode");

  /***************************************************************************/
  
  /* Check energy and cosines */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "r", "", *u**u+*v**v+*w**w - 1.0, -1E-5, 1E-5);

  /* Rotate direction cosines around a random azimuthal angle */
      
  AziRot(mu, u, v, w, id);      
}

/*****************************************************************************/

#endif
