/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : uresfactor.c                                   */
/*                                                                           */
/* Created:       2011/01/08 (JLe)                                           */
/* Last modified: 2014/05/22 (JLe)                                           */
/* Version:       2.1.21                                                     */
/*                                                                           */
/* Description: Samples unresolved resonance cross section from probability  */
/*              table.                                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UresFactor:"

/*****************************************************************************/

double UresFactor(long rea, double E, long id)
{
  long nuc, mat, urs, n, ptr, ptn, erg, ne, M, k0, k1, mt, type;
  double f, g, r, f0, f1, rnd, xs0, xs, adens, Emin, Emax;

  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  
  /* Get pointer to probability table data */

  if ((urs = (long)RDB[rea + REACTION_PTR_URES]) < VALID_PTR)
    return 1.0;

  /* Check boundaries */

  if ((E < RDB[rea + REACTION_URES_EMIN]) ||
      (E > RDB[rea + REACTION_URES_EMAX]))
    return 1.0;

  /* Check existing data */
 
  if ((f = TestValuePair(urs + URES_PTR_PREV_FACT, E, id)) > -INFTY)
    return f;
    
  /* Get reaction mt */

  mt = (long)RDB[rea + REACTION_MT];

  /* Check special */

  if ((mt == 101) || (mt == 1))
    {
      /***** Microscopic total or total absorption xs ************************/

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check ures sampling flag */

      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
	return 1.0;

      /* Get total cross section */
      
      xs0 = TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);
      xs = xs0;
      
      /* Check value (zero is not allowed) */

      CheckValue(FUNCTION_NAME, "xs0 (1)", "", xs0, 0.0, INFTY);

      /* Check mt */

      if (mt == 101)
	{
	  /* Total absorption, adjust for capture */

	  rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
	  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

	  xs = xs + MicroXS(rea, E, id);
	  xs = xs - TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);

	  /* Adjust for capture to isomeric state */

	  if ((rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS_ISO]) > VALID_PTR)
	    {
	      xs = xs + MicroXS(rea, E, id);
	      xs = xs - TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);
	    }
	}
      else
	{
	  /* Adjust for elastic */

	  rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
	  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

	  xs = xs + MicroXS(rea, E, id);
	  xs = xs - TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);
	  
	  /* Adjust for capture */
	  
	  rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
	  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

	  xs = xs + MicroXS(rea, E, id);
	  xs = xs - TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);

	  /* Adjust for capture to isomeric state */

	  if ((rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS_ISO]) > VALID_PTR)
	    {
	      xs = xs + MicroXS(rea, E, id);
	      xs = xs - TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);
	    }
	  
	  /* Adjust for fission */
	  
	  if ((rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS]) > VALID_PTR)
	    {
	      xs = xs + MicroXS(rea, E, id);
	      xs = xs - TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);
	    }
	}
      
      /* Calculate factor */

      if (xs0 > 0.0)
	f = xs/xs0;
      else
	f = 1.0;
            
      /***********************************************************************/
    }
  else if (mt < 0)
    {
      /***** Macroscopic reaction cross sections *****************************/

      /* Get total cross section */

      xs0 = TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);
      xs = xs0;

      /* Check value */
      
      CheckValue(FUNCTION_NAME, "xs0 (3)", "", xs0, 0.0, INFTY);

      /* Get pointer to material */

      mat = (long)RDB[rea + REACTION_PTR_MAT];
      CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

      /* Avoid compiler warning */

      ptr = -1;

      /* Get pointer to ures list  */
      
      if (mt == MT_MACRO_TOTXS)      
	ptr = (long)RDB[mat + MATERIAL_PTR_TOT_URES_LIST];
      else if (mt == MT_MACRO_ABSXS)      
	ptr = (long)RDB[mat + MATERIAL_PTR_ABS_URES_LIST];
      else if (mt == MT_MACRO_ELAXS)      
	ptr = (long)RDB[mat + MATERIAL_PTR_ELA_URES_LIST];
      else if (mt == MT_MACRO_FISSXS)      
	ptr = (long)RDB[mat + MATERIAL_PTR_FISS_URES_LIST];
      else if (mt == MT_MACRO_FISSE)    
	ptr = (long)RDB[mat + MATERIAL_PTR_FISS_URES_LIST];
      else if (mt == MT_MACRO_NSF)    
	ptr = (long)RDB[mat + MATERIAL_PTR_FISS_URES_LIST];
      else
	Die(FUNCTION_NAME, "Invalid mt %ld", mt);
      
      /* Check Pointer */

      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Reset reaction pointer (rewind list) */
	      
      rea = -1;

      /* Loop over reactions */

      while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
	{
	  /* Check reaction pointer */

	  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
	  
	  /* Check threshold */

	  if (E < Emin)
	    {
	      /* Break loop */

	      break;
	    }
	  else
	    {
	      /* Set multiplier */

	      if (mt == MT_MACRO_FISSE)    
		g = RDB[rea + REACTION_Q]*RDB[DATA_NORM_U235_FISSE]/U235_FISSQ; 
	      else if (mt == MT_MACRO_NSF)
		{
		  /* Get pointer to total nubar data */

		  ptn = (long)RDB[rea + REACTION_PTR_TNUBAR];
		  CheckPointer(FUNCTION_NAME, "(ptn)", DATA_ARRAY, ptr);

		  /* Get multiplier */
		  
		  g = Nubar(ptn, E, id);

		  /* Subtract delayed nubar data if delayed neutron emission */
		  /* is off */
	  
		  if (((long)RDB[DATA_USE_DELNU] == NO) && 
		      ((ptn = (long)RDB[rea + REACTION_PTR_DNUBAR]) 
		       > VALID_PTR))
		    g = g - Nubar(ptn, E, id);
		}
	      else
		g = RDB[rea + REACTION_WGT_F];
	      
	      /* Adjust cross section */
	  
	      xs = xs + g*adens*MicroXS(rea, E, id);
	      xs = xs - g*adens*TestValuePair(rea + REACTION_PTR_PREV_XS0, 
					      E, id);
	    }
	}

      /* Calculate factor */

      if (xs0 > 0.0)
	f = xs/xs0;
      else
	f = 0.0;
      
      /***********************************************************************/
    }
  else
    {
      /***** Microscopic reaction cross section ******************************/

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Check ures sampling flag */
      
      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
	return 1.0;

      /* Check if random number is already sampled */

      if ((rnd = TestValuePair(urs + URES_PTR_RND, E, id)) == -INFTY)
	{
	  /* Sample new value */
	  
	  rnd = RandF(id);
	  
	  /* Store value */
	  
	  StoreValuePair(urs + URES_PTR_RND, E, rnd, id);
	}

      /* Get pointer to energy grid */
      
      erg = (long)RDB[urs + URES_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
      
      /* Get number of grid points */
      
      ne = (long)RDB[erg + ENERGY_GRID_NE];
      
      /* Get interval */
      
      if ((n = GridSearch(erg, E)) < 0)
	f = 1.0;
      else
	{
	  /* Check index */
	  
	  CheckValue(FUNCTION_NAME, "n", "", n, 0.0, ne - 1);
	  
	  /* Check boundaries */
	  
	  if ((n < 0) || (n > ne - 2))
	    f = 1.0;
	  else
	    {
	      /* Get number of probabilities */
	      
	      M = (long)RDB[urs + URES_NP];
	      
	      /* Pointer to probability data */
	      
	      ptr = (long)RDB[urs + URES_PTR_PROB];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	      /* Find intervals */
	      
	      for (k0 = 0; k0 < M; k0++)
		if (rnd < RDB[ptr + n*M + k0])
		  break;
	      
	      for (k1 = 0; k1 < M; k1++)
		if (rnd < RDB[ptr + (n + 1)*M + k1])
		  break;
	      
	      /* Check values */
	      
	      CheckValue(FUNCTION_NAME, "k0", "", rnd, 0.0, 
			 RDB[ptr + n*M + k0]);
	      CheckValue(FUNCTION_NAME, "k1", "", rnd, 0.0, 
			 RDB[ptr + (n + 1)*M + k1]);

	      /* Get pointer to energy grid data */

	      ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA]; 
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Get interpolation type */

	      type = (long)RDB[rea + REACTION_URES_INT];

	      /* Avoid compiler warning */

	      r = -1.0;

	      /* Calculate interpolation factor */

	      if (type == 2)
		r = (E - RDB[ptr + n])/(RDB[ptr + n + 1] - RDB[ptr + n]);
	      else if (type == 5)
		r = log(E/RDB[ptr + n])/log(RDB[ptr + n + 1]/RDB[ptr + n]);
	      else
		Die(FUNCTION_NAME, "Invalid interpolation type");

	      /* Check factor */

	      CheckValue(FUNCTION_NAME, "r", "", r, 0.0, 1.0);

	      /* Pointer to factors */
	      
	      ptr = (long)RDB[urs + URES_PTR_FACT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Get values */
	      
	      f0 = RDB[ptr + n*M + k0];
	      f1 = RDB[ptr + (n + 1)*M + k1];
	      
	      /* Check interpolation scheme and interpolate */
	      
	      if (type == 2)
		f = r*(f1 - f0) + f0;
	      else if (type == 5)
		f = exp(r*log(f1/f0) + log(f0));
	      
	      /* Pointer to maximum factors */

	      ptr = (long)RDB[urs + URES_PTR_MAXF];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Check */

	      CheckValue(FUNCTION_NAME, "f", "", f, -50.0, RDB[ptr + n]);
	    }
	}
      
      /* Check mode */
      
      if ((long)RDB[urs + URES_IFF] == 0)
	{
	  /* Get un-adjusted cross section */

	  xs0 = TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);

	  /* Check value (may be zero for fission channel) */

	  CheckValue(FUNCTION_NAME, "xs0 (5)", "", xs0, 0.0, INFTY);
	  
	  /* Calculate factor */

	  if (xs0 > 0.0)
	    f = f/xs0;
	  else
	    f = 1.0;
	}

      /***********************************************************************/
    }

  /* Check factor */

  CheckValue(FUNCTION_NAME, "f", "", f, -50.0, INFTY);

  /* Remember value */

  StoreValuePair(urs + URES_PTR_PREV_FACT, E, f, id);

  /* Return factor */

  return f;
}

/*****************************************************************************/
