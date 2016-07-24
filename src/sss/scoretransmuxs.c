#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoretransmuxs.c                               */
/*                                                                           */
/* Created:       2011/04/21 (JLe)                                           */
/* Last modified: 2015/05/19 (TVi)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Scores transmutation cross sections and fission product      */
/*              yield weights for burnup calculation                         */
/*                                                                           */
/* Comments: TODO: add isomeric br calculation                               */
/*                                                                           */
/*           - Noi CheckValue() -funktiot pitää sit kattoa joka rutiinissa   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreTransmuXS:"

/*****************************************************************************/

void ScoreTransmuXS(double flx, long mat, double E, double wgt, long id)
{
  long rea, nuc, loc0, ptr, erg, i, dep, ncol, gcu;
  double val, E0, E1, E2, f, g, T, Er;

  /* Check burnup mode */

 if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO)
   return;

 /* Check if in coefficient calculation */

 if ((long)RDB[DATA_COEF_CALC_IDX] > 0)
   return;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Check material burn flag */

  if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
    return;  

  /* Check parent type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    Die(FUNCTION_NAME, "Parent type");

  /* Check flux */

  CheckValue(FUNCTION_NAME, "flx", "", flx, ZERO, INFTY);

  /* Check for B1 correction in flux */

  if ((long)RDB[DATA_B1_BURNUP_CORR] == YES)
    {
      /* Check mode */
      /*
      if ((long)RDB[DATA_MICRO_CALC_MODE] == MICRO_CALC_MODE_END)
        Error(0, "B1 burnup correction does not work in mode 1 right now");
      */

      Die(FUNCTION_NAME, "Tää ei kuitenkaan toimi");

      /* Get collision number */

      ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
      ncol = (long)GetPrivateData(ptr, id);

      /* Check for multiple levels */

      if ((long)RDB[DATA_MULTI_LEVEL_GCU] == YES)
	Die(FUNCTION_NAME, "Doesn't work with multiple levels");

      /* Get universe pointer and adjust flux */

      if ((gcu = TestValuePair(DATA_GCU_PTR_UNI, ncol, id)) > VALID_PTR)
        flx = flx*B1FluxCorr(gcu, E);
    }  

  /* Score burn flux */

  ptr = (long)RDB[mat + MATERIAL_PTR_BURN_FLUX];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(flx, wgt, ptr, id, 0);

  /***************************************************************************/

  /***** Score spectrum ******************************************************/

  /* Add to total flux */

  ptr = (long)RDB[mat + MATERIAL_PTR_FLUX_SPEC_SUM];
  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
  AddPrivateRes(ptr, flx*wgt, id);

  /* Check for spectrum-collapse mode */
  
  if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
    {
      /* Pointer to unionized grid */
      
      erg = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
      
      /* Get interpolation factor */
      
      f = GridFactor(erg, E, id);
      
      /* Separate integer and decimal parts */
      
      i = (long)f;
      f = f - (double)i;
      
      /* Check values */
      
      CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
      CheckValue(FUNCTION_NAME, "i", "", i, 0, 
		 (long)RDB[erg + ENERGY_GRID_NE] - 1);
      
      /* Score flux (NOTE: toi jälkimmäinen voi kosahtaa?) */
      
      ptr = (long)RDB[mat + MATERIAL_PTR_FLUX_SPEC];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      
      AddPrivateRes(ptr + i, (1.0 - f)*flx*wgt, id);
      AddPrivateRes(ptr + i + 1, f*flx*wgt, id);
      
      /* Exit subroutine of not in region */
      
      if ((E < RDB[DATA_BU_URES_EMIN]) || (E > RDB[DATA_BU_URES_EMAX]))
	return;
    }

  /***************************************************************************/

  /***** Score transmutation rates *******************************************/

  /* Score cross sections */
  
  loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

  /* Loop over reactions */

  i = 0;
  while ((dep = ListPtr(loc0, i++)) > VALID_PTR)
    {
      /* Check energy */
      
      if (E < RDB[dep + DEP_TRA_E0])
	break;
      
      /* Pointer to reaction data */
      
      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Check mt */

      if ((long)RDB[rea + REACTION_MT] < 16)
	Die(FUNCTION_NAME, "Error in mt");

      /* Pointer to nuclide */

      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
      
      /* Check mode and ures boundaries */

      if (((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO) ||
	  (((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES) &&
	   ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR) &&
	   (E > RDB[rea + REACTION_URES_EMIN]) && 
	   (E < RDB[rea + REACTION_URES_EMAX])))
	{
	  /* Get cross section */

	  if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_NONE)
	    val = MicroXS(rea, E, id);
	  else if ((T = GetTemp(mat, id)) < ZERO)
	    val = MicroXS(rea, E, id);
	  else	    
	    val = DopMicroXS(mat, rea, E, &Er, T, id);
      
	  /* Score */
      
	  ptr = (long)RDB[dep + DEP_TRA_PTR_RESU];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  AddPrivateRes(ptr, val*flx*wgt, id);
	}
    }

  /***************************************************************************/

  /***** Score fission yields ************************************************/

  /* Pointer to fission yields */

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST]) > VALID_PTR)
    {
      /* Loop over reactions */
      
      i = 0;
      while ((dep = ListPtr(loc0, i++)) > VALID_PTR)
	{
	  /* Check energy */
	  
	  if (E < RDB[dep + DEP_TRA_E0])
	    break;

	  /* Pointer to reaction data */

	  rea = (long)RDB[dep + DEP_TRA_PTR_REA];
	  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

	  /* Pointer to nuclide */

	  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
	  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	  /* Check mode and ures boundaries */
	  
	  if (((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == NO) ||
	      (((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES) &&
	       ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR) &&
	       (E > RDB[rea + REACTION_URES_EMIN]) && 
	       (E < RDB[rea + REACTION_URES_EMAX])))
	    {
	      /* Get interpolation energies */

	      E0 = RDB[rea + REACTION_FISSY_IE0];
	      E1 = RDB[rea + REACTION_FISSY_IE1];
	      E2 = RDB[rea + REACTION_FISSY_IE2];
	      
	      /* Check values */
	      
	      CheckValue(FUNCTION_NAME, "E1" ,"", E1, E0, E2);
	      CheckValue(FUNCTION_NAME, "E" ,"", E, E0, INFTY);
	      
	      /* Check interval */
	      
	      if (E < E2)
		{
		  /* Interpolation factor */
		  
		  if (E < E1)
		    {
		      /* Below interpolation energy */
		      
		      if (E0 < 0.0)
			g = 1.0;
		      else
			g = (E - E0)/(E1 - E0);
		    }
		  else
		    {
		      /* Above interpolation energy */
		      
		      if (E2 > 1000.0)
			g = 1.0;
		      else
			g = (E2 - E)/(E2 - E1);
		    }
		  
		  /* Check factor */
		  
		  CheckValue(FUNCTION_NAME, "g", "", g, 0.0, 1.0);
		  
		  /* Check parent */
		  
		  if ((ptr = (long)RDB[rea + REACTION_PTR_BRANCH_PARENT]) 
		      < VALID_PTR)
		    ptr = rea;
		  
		  /* Get cross section */

		  if ((long)RDB[mat + MATERIAL_TMS_MODE] == TMS_MODE_NONE)
		    val = MicroXS(ptr, E, id);
		  else if ((T = GetTemp(mat, id)) < ZERO)
		    val = MicroXS(ptr, E, id);
		  else		    
		    val = DopMicroXS(mat, ptr, E, &Er, T, id);

		  /* Pointer to data */
		  
		  ptr = (long)RDB[dep + DEP_TRA_PTR_RESU];
		  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
		  AddPrivateRes(ptr, g*val*flx*wgt, id);
		}
	    }
	}
    }
  
  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
