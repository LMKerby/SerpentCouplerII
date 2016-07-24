#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : calculatetransmuxs.c                           */
/*                                                                           */
/* Created:       2011/04/23 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Calculates one-group transmutation cross sections and        */
/*              fission yields                                               */
/*                                                                           */
/* Comments: - Values read from RES2 data block are truncated to 6 decimals  */
/*             to avoid round-off errors in reproducible MPI mode.           */
/*                                                                           */
/*           - B1-lasku ei toimi ures-datan kanssa                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CalculateTransmuXS:"

/*****************************************************************************/

void CalculateTransmuXS(long mat, long id)
{
  long loc0, i, loc1, dep, rea, ptr, i0, ne, erg, erg0, sz, i1, nuc, rea1;
  double g, sum, E, E0, E1, E2, Emin, Emax, *xs, *spec, flx;

  /* Check burnup mode and burn flag */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) ||
      (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)))
    return;
  
  /* Check decay only mode */

  if ((long)RDB[DATA_BURN_DECAY_CALC] == YES)
    return;

  /* Check divisor type */

  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
    Die(FUNCTION_NAME, "Divided parent material");

  /* Get sum of flux spectrum */
  
  ptr = (long)RDB[mat + MATERIAL_PTR_FLUX_SPEC_SUM];
  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
  flx = Truncate(GetPrivateRes(ptr), 6);

  /* Check spectrum-collapse method */
  
  if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
    {
      /* Pointer to unionized grid data */

      erg0 = (long)RDB[DATA_ERG_PTR_UNIONIZED_NGRID];
      CheckPointer(FUNCTION_NAME, "(erg0)", DATA_ARRAY, erg0);

      /* Number of grid points */

      sz = (long)RDB[erg0 + ENERGY_GRID_NE];

      /* Pointer to points */

      loc1 = (long)RDB[erg0 + ENERGY_GRID_PTR_DATA];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Allocate memory for spectrum */
      /*
      spec = Mem(MEM_ALLOC, sz, sizeof(double));  
      */
      spec = WorkArray(DATA_PTR_WORK_PRIVA_GRID1, PRIVA_ARRAY, sz, id);

      /* Get pointer to flux spectrum */

      loc0 = (long)RDB[mat + MATERIAL_PTR_FLUX_SPEC];
      CheckPointer(FUNCTION_NAME, "(loc0)", RES2_ARRAY, loc0);
      
      /* Reset sum */

      sum = 0.0;

      /* Read data */
  
      for (i = 0; i < sz; i++)
	{
	  spec[i] = Truncate(GetPrivateRes(loc0 + i), 6);
	  sum = sum + spec[i];
	}

      /* Check sum */

      if (flx > 0.0)
	if (fabs(sum/flx - 1.0) > 1E-5)
	  Die(FUNCTION_NAME, "Error in sum");

      /* Allocate memory for temporary array if microscopic data is not */
      /* reconstructed */

      if ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == NO)
	{
	  xs = WorkArray(DATA_PTR_WORK_PRIVA_GRID2, PRIVA_ARRAY, sz, id);
	  /*
	    xs = Mem(MEM_ALLOC, sz, sizeof(double));	
	  */
	}
      else
	xs = NULL;
    }
  else
    {
      /* Reset pointers and variables */

      erg0 = -1;
      sz = 0;
      loc0 = -1;
      loc1 = -1;
      spec = NULL;
      xs = NULL;
    }

  /***************************************************************************/

  /***** Transmutation cross sections ****************************************/

  /* Pointer to depletion list (onko toi lista aina olemassa?) */
  
  dep = (long)RDB[mat + MATERIAL_PTR_DEP_TRA_LIST];
  CheckPointer(FUNCTION_NAME, "(dep)", DATA_ARRAY, dep);

  /* Loop over reactions */
      
  while (dep > VALID_PTR)
    {
      /* Pointer to reaction */
      
      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
      
      /* Check mt */
      
      if ((long)RDB[rea + REACTION_MT] < 16)
	Die(FUNCTION_NAME, "Error in mt");

      /* Pointer to nuclide data */
      
      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set energy intervals */
      
      Emin = RDB[rea + REACTION_URES_EMIN];
      Emax = RDB[rea + REACTION_URES_EMAX];

      /* Check that boundaries are set */

      if (((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES) && (Emin >= Emax) &&
	  ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR))
	Die(FUNCTION_NAME, "Error in ures boundaries 1 %s %ld %E %E",
	    GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT],
	    Emin, Emax);

      /* Pointer to data */
     
      ptr = (long)RDB[dep + DEP_TRA_PTR_RESU];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      
      /* Get value and truncate */

      sum = Truncate(GetPrivateRes(ptr), 6);

      /* Check spectrum-collapse method */

      if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
	{
	  /* Check tallied value */
	  
	  if (sum > 0.0)
	    {
	      /* Check ures samling */
	      
	      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
		Die(FUNCTION_NAME, "Value should be zero");
	      else if ((long)RDB[rea + REACTION_PTR_URES] < VALID_PTR)
		Die(FUNCTION_NAME, "Value should be zero");
	    }

	  /* Get pointer to energy grid */
	  
	  erg = (long)RDB[rea + REACTION_PTR_EGRID];
	  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

	  /* Get pointer to grid data */

	  erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
	  CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

	  /* Get first energy point and number of points */
	  
	  i0 = (long)RDB[rea + REACTION_XS_I0];
	  ne = (long)RDB[rea + REACTION_XS_NE];

	  /* Pointer to cross section data */
	  
	  ptr = (long)RDB[rea + REACTION_PTR_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  /* Check reconstruction option */
	  
	  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == YES)
	    {
	      /* Copy pointer */

	      xs = &WDB[ptr];
	      
	      /* Set xs start index to zero */
	      
	      i1 = 0;
	    }
	  else
	    {
	      /* Check pointer */

	      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	      /* Reconstruct data on new array */
	      
	      InterpolateData(&RDB[loc1], xs, sz, &RDB[erg + i0], 
			      &RDB[ptr], ne, 0, &i1, &ne);
	      
	      /* One step back if not at the beginning */
	      
	      if (i1 > 0)
		i1--;
	      
	      /* Set energy grid start index equal to xs index */
	      
	      i0 = i1;
	    }

	  /* Avoid compiler warning */

	  E = -1.0;

	  /* Loop over data and add to sum */
	  
	  for (i = 0; i < ne; i++)
	    {
	      /* Get energy */

	      if (i0 + i > sz - 1)
		Die(FUNCTION_NAME, "Energy array dimension exceeded");
	      else
		E = RDB[loc1 + i0 + i];

	      /* Compare to spectrum boundaries (NOTE: Yhtäsuuruusmerkki */
	      /* tarvitaan) */

	      if (((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO) ||
		  (E <= Emin) || (E >= Emax))
		sum = sum + spec[i0 + i]*xs[i1 + i];
	    }
	}
   
      /* Divide sum */

      if (flx > 0.0)  
	sum = sum/flx;
      else if (sum > 0.0)
	Die(FUNCTION_NAME, "Error in sums");	    

      /* Store value */
      
      if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] < VALID_PTR)
	Die(FUNCTION_NAME, "Pointer error");
      else
	StoreValuePair(rea + REACTION_PTR_TRANSMUXS, mat, sum, id);

      /* Next reaction */

      dep = NextItem(dep);
    }
  
  /***************************************************************************/

  /***** Partial fission cross sections **************************************/

  /* Pointer to fission list (tätä ei välttämättä ole) */
  
  dep = (long)RDB[mat + MATERIAL_PTR_DEP_FISS_LIST];

  /* Loop over reactions */
  
  while (dep > VALID_PTR)
    {
      /* Pointer to reaction */
      
      rea = (long)RDB[dep + DEP_TRA_PTR_REA];
      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

      /* Pointer to nuclide data */
      
      nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Set energy intervals */
      
      Emin = RDB[rea + REACTION_URES_EMIN];
      Emax = RDB[rea + REACTION_URES_EMAX];

      /* Check that boundaries are set */

      if (((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == YES) && (Emin >= Emax) &&
	  ((long)RDB[rea + REACTION_PTR_URES] > VALID_PTR))
	Die(FUNCTION_NAME, "Error in ures boundaries 2 %s %ld %E %E",
	    GetText(nuc + NUCLIDE_PTR_NAME), (long)RDB[rea + REACTION_MT],
	    Emin, Emax);
      
      /* Pointer to total fission */
      
      rea1 = (long)RDB[nuc + NUCLIDE_PTR_TOTFISS_REA];
      CheckPointer(FUNCTION_NAME, "(rea1)", DATA_ARRAY, rea1);
      
      /* Get interpolation energies */
      
      E0 = RDB[rea + REACTION_FISSY_IE0];
      E1 = RDB[rea + REACTION_FISSY_IE1];
      E2 = RDB[rea + REACTION_FISSY_IE2];
      
      /* Check values */
      
      CheckValue(FUNCTION_NAME, "E1" ,"", E1, E0, E2);
      
      /* Pointer to data */
      
      ptr = (long)RDB[dep + DEP_TRA_PTR_RESU];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      
      /* Get value and truncate */
	  
      sum = Truncate(GetPrivateRes(ptr), 6);
	 
      /* Check spectrum-collapse method */
      
      if ((long)RDB[DATA_BU_SPECTRUM_COLLAPSE] == YES)
	{
	  /* Check tallied value */
	  
	  if (sum > 0.0)
	    {
	      /* Check ures samling */
	      
	      if ((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO)
		Die(FUNCTION_NAME, "Value should be zero");
	      else if ((long)RDB[rea + REACTION_PTR_URES] < VALID_PTR)
		Die(FUNCTION_NAME, "Value should be zero");
	    }

	  /* Check pointer to parent reaction */
	  
	  if ((ptr = (long)RDB[rea + REACTION_PTR_BRANCH_PARENT]) > VALID_PTR)
	    {
	      /* Get first energy point and number of points */
	      
	      i0 = (long)RDB[ptr + REACTION_XS_I0];
	      ne = (long)RDB[ptr + REACTION_XS_NE];
	      
	      /* Get pointer to energy grid */
	      
	      erg = (long)RDB[ptr + REACTION_PTR_EGRID];
	      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
	      
	      /* Get pointer to grid data */
	      
	      erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
	      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
	      
	      /* Pointer to cross section data */
	      
	      ptr = (long)RDB[ptr + REACTION_PTR_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    }
	  else
	    {
	      /* Get first energy point and number of points */
	      
	      i0 = (long)RDB[rea + REACTION_XS_I0];
	      ne = (long)RDB[rea + REACTION_XS_NE];
	      
	      /* Get pointer to energy grid */
	      
	      erg = (long)RDB[rea + REACTION_PTR_EGRID];
	      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
	      
	      /* Get pointer to grid data */
	      
	      erg = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
	      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);
	      
	      /* Pointer to cross section data */
	      
	      ptr = (long)RDB[rea + REACTION_PTR_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    }
	  
	  /* Check reconstruction option */
	  
	  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == YES)
	    {
	      /* Copy pointer */
	      
	      xs = &WDB[ptr];
	      
	      /* Set xs start index to zero */
	      
	      i1 = 0;
	    }
	  else
	    {
	      /* Reconstruct data on new array */
	      
	      InterpolateData(&RDB[loc1], xs, sz, &RDB[erg + i0], 
			      &RDB[ptr], ne, 0, &i1, &ne);
	      
	      /* One step back if not at the beginning */
	      
	      if (i1 > 0)
		i1--;
	      
	      /* Set energy grid start index equal to xs index */
	      
	      i0 = i1;
	    }

	  /* Avoid compiler warning */

	  E = -1.0;
	  
	  /* Loop over data and calculate sum */
	  
	  for (i = 0; i < ne; i++)
	    {
	      /* Get energy */

	      if (i0 + i > sz - 1)
		Die(FUNCTION_NAME, "Energy array dimension exceeded");
	      else
		E = RDB[loc1 + i0 + i];
	      
	      /* Compare to spectrum boundaries (NOTE: Yhtäsuuruusmerkki */
	      /* tarvitaan) */

	      if (((long)RDB[nuc + NUCLIDE_URES_SAMPLING] == NO) ||
		  (E <= Emin) || (E >= Emax))
		{
		  /* Compare to interval boundaries */
		  
		  if ((E >= E0) && (E <= E2))
		    {
		      /* Calculate interpolation factor */
		      
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
		      
		      CheckValue(FUNCTION_NAME, "f", "", g, 0.0, 1.0);
		      
		      /* Add to sum */
		      
		      sum = sum + g*spec[i0 + i]*xs[i1 + i];
		    }
		}
	    }
	}

      /* Divide sum */
      
      if (flx > 0.0)	  
	sum = sum/flx;
      else if (sum > 0.0)
	Die(FUNCTION_NAME, "Error in sums");	    

      /* Check value */
      
      CheckValue(FUNCTION_NAME, "sum", "", sum, 0.0, 1E+6);

      /* Store value */
      
      if ((long)RDB[rea + REACTION_PTR_TRANSMUXS] < VALID_PTR)
	Die(FUNCTION_NAME, "Pointer error");
      else
	StoreValuePair(rea + REACTION_PTR_TRANSMUXS, mat, sum, id);
      
      /* Add to total fission */

      if ((long)RDB[rea1 + REACTION_PTR_TRANSMUXS] < VALID_PTR)
	Die(FUNCTION_NAME, "Pointer error");
      else
	AddValuePair(rea1 + REACTION_PTR_TRANSMUXS, mat, sum, id);

      /* Next reaction */

      dep = NextItem(dep);
    }

  /***************************************************************************/

  /*
  if (((long)RDB[DATA_OPTI_RECONSTRUCT_MICROXS] == NO) && (xs != NULL))
    Mem(MEM_FREE, xs);

  if (spec != NULL)
    Mem(MEM_FREE, spec);
  */

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
