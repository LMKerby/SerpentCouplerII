/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : macroxs.c                                      */
/*                                                                           */
/* Created:       2011/01/02 (JLe)                                           */
/* Last modified: 2015/06/09 (TVi)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Interpolates macroscopic cross section                       */
/*                                                                           */
/* Comments: - Rutiinia muutettu radikaalisti 2.11.2011 (2.0.37)             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MacroXS:"

/*****************************************************************************/

double MacroXS(long rea0, double E, long id)
{
  long i, ptr, ptn, rea, erg, ne, mat, nuc, ncol, mt;
  double xs0, xs1, xs, adens, f, mult, Emin, Emax, Er, T;
  
  /* Check Pointer */

  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

  /* Get pointer to data */

  if ((ptr = (long)RDB[rea0 + REACTION_PTR_XS]) > VALID_PTR)
    {
      /***********************************************************************/

      /***** Interpolate pre-calculated data *********************************/

      /* Test existing data */

      if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
	return xs;
      
      /* Get mt */

      mt = (long)RDB[rea0 + REACTION_MT];
      
      /* Get pointer to material */
      
      mat = (long)RDB[rea0 + REACTION_PTR_MAT];
      CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

#ifdef DEBUG

      /* Sanity check for TMS */

      if (((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_MG) ||
	  (((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_CE) &&
	   (mt != MT_MACRO_TMP_MAJORANTXS)))
	Die(FUNCTION_NAME, "Pre-calculated data in %s %ld", 
	    GetText(mat + MATERIAL_PTR_NAME), mt);

#endif

      /* Reset cross sections */

      xs = 0.0;
      xs0 = 0.0;
   
      /* Get pointer to energy grid */

      erg = (long)RDB[rea0 + REACTION_PTR_EGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

      /* Get interpolation factor */

      if ((f = GridFactor(erg, E, id)) < 0)
	xs = 0.0;
      else
	{
	  /* Check interpolation factor */

	  CheckValue(FUNCTION_NAME, "f (1)", "", f, 0.0, MAX_EGRID_NE);

	  /* Separate integer and decimal parts of interpolation factor */
      
	  i = (long)f;
	  f = f - (double)i;
      
	  /* Get number of points */
      
	  ne = (long)RDB[rea0 + REACTION_XS_NE];
	  CheckValue(FUNCTION_NAME, "ne", "", ne, 2, MAX_EGRID_NE);
      
	  /* Check boundaries */
	  
	  if ((i < 0) || (i > ne - 1))
	    xs = 0.0;
	  else
	    {      
	      /* Get tabulated cross sections */
	  
	      xs0 = RDB[ptr + i];
	      xs1 = RDB[ptr + i + 1];

	      if (mt != MT_MACRO_TMP_MAJORANTXS) 
		{	      		
		  /* Interpolate in normal case */
		  
		  if (i == ne - 1)
		    xs = (1.0 - f)*xs0;
		  else
		    xs = f*(xs1 - xs0) + xs0;
		}
	      else 
		{     		
		  /* TMS-tapauksessa majoranttia ei interpoloida */
		  /* (histogrammimajorantti) */
		  
		  xs = xs0;
		}
      	    }
	}

      /* Add poison cross section */
      
      if (((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == YES) ||
	  ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == YES))
	if ((mt == MT_MACRO_TOTXS) || (mt == MT_MACRO_ABSXS))
	  xs = xs + PoisonXS(mat, E, mt, id);
      
      /* Remember non-adjusted cross section */

      xs0 = xs;

      /* Store value */

      StoreValuePair(rea0 + REACTION_PTR_PREV_XS0, E, xs0, id);
  
      /* Sample unresolved resonance probability table data */

      f = UresFactor(rea0, E, id);
      CheckValue(FUNCTION_NAME, "f (2)", "", f, -10.0, INFTY);

      /* Adjust cross section */
      
      xs = f*xs;

      /* Store cross section */

      if (xs > 0.0)
	StoreValuePair(rea0 + REACTION_PTR_PREV_XS, E, xs, id);

      /* Return interpolated value */
      
      return xs;

      /***********************************************************************/
    }

  /* Get mt */

  mt = (long)RDB[rea0 + REACTION_MT];
  
  /* Get pointer to material */
  
  mat = (long)RDB[rea0 + REACTION_PTR_MAT];
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  
  /* Check mode */

  if (((long)RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE) &&
      (mt != MT_MACRO_TMP_MAJORANTXS))
    {
      /***********************************************************************/
      
      /***** TMS-moodi ******************************************************/

#ifdef DEBUG

      /* Sanity check for TMS */

      if ((RDB[mat + MATERIAL_TMS_TMIN] == 0.0) || 
	  (RDB[mat + MATERIAL_TMS_TMIN] > RDB[mat + MATERIAL_TMS_TMAX]))
	Die(FUNCTION_NAME, "Error in temperature in %s", 
	    GetText(mat + MATERIAL_PTR_NAME));

#endif

      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      ncol = (long)GetPrivateData(ptr, id);

      /* Test existing data */
      
      if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, ncol, id)) > -INFTY)
	return xs;
      
      /* Reset cross sections */

      xs = 0.0;
      xs0 = 0.0;
   
      /* Get material temperature for on-the-fly temperature treatment */

      if ((T = GetTemp(mat, id)) > 0.0)
	{
	  /* Get pointer to partial list */

	  ptr = (long)RDB[rea0 + REACTION_PTR_PARTIAL_LIST];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Reset reaction pointer (rewind list) */
	  
	  rea = -1;

	  /* Loop over reactions */
	  
	  while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
	    {
	      /* Check reaction pointer */
	      
	      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
	      
	      /* Set multiplier */
	      
	      if (mt == MT_MACRO_FISSE)
		mult = RDB[rea + REACTION_Q]
		  *RDB[DATA_NORM_U235_FISSE]/U235_FISSQ; 
	      else if (mt == MT_MACRO_NSF)
		{
		  /* Get pointer to total nubar data */

		  ptn = (long)RDB[rea + REACTION_PTR_TNUBAR];
		  CheckPointer(FUNCTION_NAME, "(ptn)", DATA_ARRAY, ptr);

		  /* Get multiplier */
		  
		  mult = Nubar(ptn, E, id);

		  /* Subtract delayed nubar data if delayed neutron */
		  /* emission is off */
	  
		  if (((long)RDB[DATA_USE_DELNU] == NO) && 
		      ((ptn = (long)RDB[rea + REACTION_PTR_DNUBAR]) 
		       > VALID_PTR))
		    mult = mult - Nubar(ptn, E, id);
		}
	      else if (mt == MT_MACRO_INLPRODXS)
		mult = RDB[rea + REACTION_WGT_F] - 1.0;
	      else
		mult = RDB[rea + REACTION_WGT_F];

	      /* Pointer to nuclide */

	      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
	      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
		
	      /* Add to cross section */
	  
	      xs = xs + mult*adens*DopMicroXS(mat, rea, E, &Er, T, id);
	  
	      /* Add to non-adjusted cross section */
	      /* Muutettu TestValuePairiin E -> Er (TVi 2015-06-03) */
	      /* Koska DopMicroXS kutsuu MicroXSsaa nimenomaan suhteellisella
		 energialla */

	      xs0 = xs0 + adens*TestValuePair(rea + REACTION_PTR_PREV_XS0, 
					      Er, id);
	      
	      /* Check energy cut-off */
	      
	      if (E < Emin)
		break;
	    }
	  
	  /* Store non-adjusted cross section */
	  
	  StoreValuePair(rea0 + REACTION_PTR_PREV_XS0, E, xs0, id);
	  
	  /* Check if ures data exists */
	  
	  if ((ptr = (long)RDB[rea0 + REACTION_PTR_URES]) > VALID_PTR)
	    {
	      /* Calculate ures factor */
	      
	      if (xs0 > 0.0)
		f = xs/xs0;
	      else
		f = 1.0;
	      
	      /* Check value */
	      
	      CheckValue(FUNCTION_NAME, "f (3)", "", f, ZERO, INFTY);
	      
	      /* Store factor */
	      
	      StoreValuePair(ptr + URES_PTR_PREV_FACT, E, f, id);
	    }
	  	  
	  /* Store cross section */

	  if (xs > 0.0)
	    StoreValuePair(rea0 + REACTION_PTR_PREV_XS, ncol, xs, id);

	  /* Return interpolated value */
	  
	  return xs;
	}

      /***********************************************************************/
    }

  /***************************************************************************/
  
  /***** Calculate sum of partials *******************************************/

  /* Test existing data */

  if ((xs = TestValuePair(rea0 + REACTION_PTR_PREV_XS, E, id)) > -INFTY)
    return xs;
  
  /* Reset cross sections */
  
  xs = 0.0;
  xs0 = 0.0;
  
  /* Get pointer to partial list */
  
  ptr = (long)RDB[rea0 + REACTION_PTR_PARTIAL_LIST];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Reset reaction pointer (rewind list) */
  
  rea = -1;
  
  /* Loop over reactions */

  while (NextReaction(ptr, &rea, &adens, &Emin, &Emax, id) > VALID_PTR)
    {
      /* Check reaction pointer */

      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Set multiplier */

      if (mt == MT_MACRO_FISSE)
	mult = RDB[rea + REACTION_Q]*RDB[DATA_NORM_U235_FISSE]/U235_FISSQ; 
      else if (mt == MT_MACRO_NSF)
	{
	  /* Get pointer to total nubar data */

	  ptn = (long)RDB[rea + REACTION_PTR_TNUBAR];
	  CheckPointer(FUNCTION_NAME, "(ptn)", DATA_ARRAY, ptr);
	  
	  /* Get multiplier */
	  
	  mult = Nubar(ptn, E, id);

	  /* Subtract delayed nubar data if delayed neutron emission is off */
	  
	  if (((long)RDB[DATA_USE_DELNU] == NO) && 
	      ((ptn = (long)RDB[rea + REACTION_PTR_DNUBAR]) > VALID_PTR))
	    mult = mult - Nubar(ptn, E, id);
	}
      else if (mt == MT_MACRO_INLPRODXS)
	mult = RDB[rea + REACTION_WGT_F] - 1.0;
      else
	mult = RDB[rea + REACTION_WGT_F];

      /* Add to cross section */

      if(mt != MT_MACRO_TMP_MAJORANTXS)
	xs = xs + mult*adens*MicroXS(rea, E, id);

      /* In case of majorantxs, use MicroMajorantXS */
      else
	xs = xs + mult*adens*MicroMajorantXS(rea, E, id);
      
      /* Add to non-adjusted cross section */
      
      xs0 = xs0 + adens*TestValuePair(rea + REACTION_PTR_PREV_XS0, E, id);
      
      /* Check energy cut-off */
      
      if (E < Emin)
	break;
    }

  /* Store non-adjusted cross section */
  
  StoreValuePair(rea0 + REACTION_PTR_PREV_XS0, E, xs0, id);
  
  /* Check if ures data exists */
  
  if ((ptr = (long)RDB[rea0 + REACTION_PTR_URES]) > VALID_PTR)
    {
      /* Calculate ures factor */
      
      if (xs0 > 0.0)
	f = xs/xs0;
      else
	f = 1.0;
      
      /* Check value */
      
      CheckValue(FUNCTION_NAME, "f (3)", "", f, ZERO, INFTY);
      
      /* Store factor */
      
      StoreValuePair(ptr + URES_PTR_PREV_FACT, E, f, id);
    }

  /* Store cross section */

  if (xs > 0.0)
    StoreValuePair(rea0 + REACTION_PTR_PREV_XS, E, xs, id);

  /* Return interpolated value */

  return xs;
  
  /****************************************************************************/
}

/******************************************************************************/
