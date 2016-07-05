/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : sumtotxs.c                                     */
/*                                                                           */
/* Created:       2010/12/29 (JLe)                                           */
/* Last modified: 2014/06/04 (JLe)                                           */
/* Version:       2.1.21                                                     */
/*                                                                           */
/* Description: Sets reaction list for partials and calculates nuclide-wise  */
/*              total and total absorption cross section                     */
/*                                                                           */
/* Comments: - Total absorption is an array of zeros if nuclide has no       */
/*             absorption channels                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "SumTotXS:"

/*****************************************************************************/

void SumTotXS(long nuc)
{
  long tot, abs, erg, sz, loc0, loc1, ptr, pte, rea, i0, ne, i, type, ty;
  long ptp;
  double max, *E, *xs, tmp[2];
  
  /* Get type */

  type = (long)RDB[nuc + NUCLIDE_TYPE];

  /* Create reaction structures */

  tot = NewItem(nuc + NUCLIDE_PTR_TOTXS, REACTION_BLOCK_SIZE);

  if (type != NUCLIDE_TYPE_PHOTON)
    abs = NewItem(nuc + NUCLIDE_PTR_SUM_ABSXS, REACTION_BLOCK_SIZE);
  else
    abs = -1;

  /* Put nuclide pointers */

  WDB[tot + REACTION_PTR_NUCLIDE] = (double)nuc;

  if (type != NUCLIDE_TYPE_PHOTON)
    WDB[abs + REACTION_PTR_NUCLIDE] = (double)nuc;

  /* Put type */

  WDB[tot + REACTION_TYPE] = (double)REACTION_TYPE_SUM;

  /* Put mt */

  if (type != NUCLIDE_TYPE_PHOTON)
    WDB[tot + REACTION_MT] = 1.0;
  else
    WDB[tot + REACTION_MT] = 501.0;

  /* Put multiplier */

  WDB[tot + REACTION_WGT_F] = 1.0;

  if (type != NUCLIDE_TYPE_PHOTON)
    WDB[abs + REACTION_WGT_F] = 1.0;

  if (type != NUCLIDE_TYPE_PHOTON)
    {
      WDB[abs + REACTION_TYPE] = (double)REACTION_TYPE_SUM;
      WDB[abs + REACTION_MT] = 101.0;
    }

  /* Allocate memory for previous values */

  AllocValuePair(tot + REACTION_PTR_PREV_XS);
  AllocValuePair(tot + REACTION_PTR_PREV_XS0);

  if (type != NUCLIDE_TYPE_PHOTON)
    {
      AllocValuePair(abs + REACTION_PTR_PREV_XS);
      AllocValuePair(abs + REACTION_PTR_PREV_XS0);
    }

  /* Put awr */

  WDB[tot + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

  if (type != NUCLIDE_TYPE_PHOTON)
    WDB[abs + REACTION_AWR] = RDB[nuc + NUCLIDE_AWR];

  /* Reset maximum */

  max = 0.0;

  /* Get pointer to energy grid data */

  if ((erg = (long)RDB[nuc + NUCLIDE_PTR_EGRID]) < VALID_PTR)
    {
      /***********************************************************************/

      /***** No common grid in use *******************************************/

      /* Reset number of grid points and poiter */
      
      sz = 0;
      E = NULL;

      /* Loop over reaction modes and unionize grid */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
	{
	  /* Check type */

	  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
	    {
	      /* Pointer to nuclide-wise grid */

	      erg = (long)RDB[rea + REACTION_PTR_EGRID];
	      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

	      /* Number of points and pointer to data */

	      ne = (long)RDB[erg + ENERGY_GRID_NE];

	      ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Add points to common grid */

	      E = AddPts(E, &sz, &RDB[ptr], ne);

	      /* Add extra point just below and above boundaries (this is */
	      /* needed to account for the cut-off of S(a,b) data) */

	      tmp[0] = 0.999999999*RDB[rea + REACTION_EMIN];
	      tmp[1] = 1.000000001*RDB[rea + REACTION_EMAX];

	      /* Cut values to nuclide limits (this is needed to prevent */
	      /* non-zero at limiting points) */

	      if (tmp[0] < RDB[nuc + NUCLIDE_EMIN])
		tmp[0] = RDB[nuc + NUCLIDE_EMIN];

	      if (tmp[1] > RDB[nuc + NUCLIDE_EMAX])
		tmp[1] = RDB[nuc + NUCLIDE_EMAX];

	      /* Add points to array */

	      E = AddPts(E, &sz, tmp, 2);
	    }

	  /* Next reaction */

	  rea = NextItem(rea);
	}

      /* Generate grid */

      erg = MakeEnergyGrid(sz, 0, 0, -1, E, EG_INTERP_MODE_LIN);

      /* Put minimum and maximum energy */

      WDB[tot + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
      WDB[tot + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];

      if (type != NUCLIDE_TYPE_PHOTON)
	{
	  WDB[abs + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
	  WDB[abs + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];
	}

      /* Reset ures energy boundaries */

      WDB[tot + REACTION_URES_EMIN] = INFTY;
      WDB[tot + REACTION_URES_EMAX] = -INFTY;

      if (type != NUCLIDE_TYPE_PHOTON)
	{
	  WDB[abs + REACTION_URES_EMIN] = INFTY;
	  WDB[abs + REACTION_URES_EMAX] = -INFTY;
	}

      /* Put pointer to grid */

      WDB[tot + REACTION_PTR_EGRID] = (double)erg;

      if (type != NUCLIDE_TYPE_PHOTON)
	WDB[abs + REACTION_PTR_EGRID] = (double)erg;

      /* Put first point and number of points */

      WDB[tot + REACTION_XS_I0] = 0.0;
      WDB[tot + REACTION_XS_NE] = sz;

      if (type != NUCLIDE_TYPE_PHOTON)
	{  
	  WDB[abs + REACTION_XS_I0] = 0.0;
	  WDB[abs + REACTION_XS_NE] = sz;
	}

      /* Allocate memory for data */

      loc0 = ReallocMem(DATA_ARRAY, sz);

      if (type != NUCLIDE_TYPE_PHOTON)
	loc1 = ReallocMem(DATA_ARRAY, sz);
      else
	loc1 = -1;

      /* Put pointer */

      WDB[tot + REACTION_PTR_XS] = (double)loc0;
      
      if (type != NUCLIDE_TYPE_PHOTON)
	WDB[abs + REACTION_PTR_XS] = (double)loc1;

      /* Allocate memory for temporary cross section array */

      xs = (double *)Mem(MEM_ALLOC, sz, sizeof(double));

      /* Loop over reaction modes and reconstruct cross sections */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];

      while (rea > VALID_PTR)
	{
	  /* Check type */

	  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
	    {
	      /* Get ty */

	      ty = (long)RDB[rea + REACTION_TY];

	      /* Pointer to nuclide-wise grid */

	      erg = (long)RDB[rea + REACTION_PTR_EGRID];
	      CheckPointer(FUNCTION_NAME, "(erg)", DATA_ARRAY, erg);

	      /* Get number of points */
	      
	      ne = (long)RDB[erg + ENERGY_GRID_NE];

	      /* Pointer to energy grid data */

	      pte = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
	      CheckPointer(FUNCTION_NAME, "(pte)", DATA_ARRAY, pte);
	      
	      /* Get pointer to cross section data */

	      ptr = (long)RDB[rea + REACTION_PTR_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Reconstruct cross section */

	      InterpolateData(E, xs, sz, &RDB[pte], &RDB[ptr], ne, 0, NULL,
			      NULL);
	      
	      /* Loop over points and add to totals */

	      for (i = 0; i < sz; i++)
		{
		  /* Add value to total */

		  WDB[loc0 + i] = RDB[loc0 + i] + xs[i];

		  /* Compare to maximum */

		  if (RDB[loc0 + i] > max)
		    max = RDB[loc0 + i];

		  /* Add value to absorption */

		  if ((ty == 0) && (type != NUCLIDE_TYPE_PHOTON))
		    WDB[loc1 + i] = RDB[loc1 + i] + xs[i];
		}	      

	      /* Add reaction to total list */

	      ptr = NewItem(tot + REACTION_PTR_PARTIAL_LIST, 
			    RLS_DATA_BLOCK_SIZE);

	      /* Allocate memory for reaction counter */
		
	      ptp = AllocPrivateData(1, PRIVA_ARRAY);
	      WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

	      WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
	      WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;
	      WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
	      WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

	      /* Add reaction to absorption list */

	      if ((ty == 0) && (type != NUCLIDE_TYPE_PHOTON))
		{
		  ptr = NewItem(abs + REACTION_PTR_PARTIAL_LIST, 
				RLS_DATA_BLOCK_SIZE);

		  /* Allocate memory for reaction counter */
		
		  ptp = AllocPrivateData(1, PRIVA_ARRAY);
		  WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

		  WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
		  WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;
		  WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
		  WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
		}
	    }
	  
	  /* Next reaction */

	  rea = NextItem(rea);
	}

      /* Free temporary arrays */

      Mem(MEM_FREE, E);
      Mem(MEM_FREE, xs);

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Common energy grid available ************************************/

      /* Put minimum and maximum energy */

      WDB[tot + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
      WDB[tot + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];

      if (type != NUCLIDE_TYPE_PHOTON)
	{
	  WDB[abs + REACTION_EMIN] = RDB[erg + ENERGY_GRID_EMIN];
	  WDB[abs + REACTION_EMAX] = RDB[erg + ENERGY_GRID_EMAX];
	}

      /* Reset ures energy boundaries */

      WDB[tot + REACTION_URES_EMIN] = INFTY;
      WDB[tot + REACTION_URES_EMAX] = -INFTY;

      if (type != NUCLIDE_TYPE_PHOTON)
	{
	  WDB[abs + REACTION_URES_EMIN] = INFTY;
	  WDB[abs + REACTION_URES_EMAX] = -INFTY;
	}

      /* Put pointer to grid */

      WDB[tot + REACTION_PTR_EGRID] = (double)erg;

      if (type != NUCLIDE_TYPE_PHOTON)
	WDB[abs + REACTION_PTR_EGRID] = (double)erg;

      /* Get number of points */

      sz = (long)RDB[erg + ENERGY_GRID_NE];

      /* Put first point and number of points */

      WDB[tot + REACTION_XS_I0] = 0.0;
      WDB[tot + REACTION_XS_NE] = sz;

      if (type != NUCLIDE_TYPE_PHOTON)
	{
	  WDB[abs + REACTION_XS_I0] = 0.0;
	  WDB[abs + REACTION_XS_NE] = sz;
	}

      /* Allocate memory for data */

      loc0 = ReallocMem(DATA_ARRAY, sz);

      if (type != NUCLIDE_TYPE_PHOTON)
	loc1 = ReallocMem(DATA_ARRAY, sz);
      else
	loc1 = -1;

      /* Put pointer */

      WDB[tot + REACTION_PTR_XS] = (double)loc0;

      if (type != NUCLIDE_TYPE_PHOTON)
	WDB[abs + REACTION_PTR_XS] = (double)loc1;

      /* Loop over reactions */

      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
	{
	  /* Check type */

	  if ((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL)
	    {
	      /* Get ty */

	      ty = (long)RDB[rea + REACTION_TY];

	      /* Get index to first point and number of points */
	      
	      i0 = (long)RDB[rea + REACTION_XS_I0];
	      ne = (long)RDB[rea + REACTION_XS_NE];

	      /* Check values */

	      CheckValue(FUNCTION_NAME, "ne", "", i0 + ne, 2, sz);

	      /* Get pointer to data */

	      ptr = (long)RDB[rea + REACTION_PTR_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Loop over points and add to total */

	      for (i = 0; i < ne; i++)
		{
		  /* Add value to total */
		  
		  WDB[loc0 + i0 + i] = RDB[loc0 + i0 + i] + RDB[ptr + i];
		  
		  /* Compare to maximum */
		  
		  if (RDB[loc0 + i0 + i] > max)
		    max = RDB[loc0 + i0 + i];
		  
		  /* Add value to absorption */
		  
		  if ((ty == 0) && (type != NUCLIDE_TYPE_PHOTON))
		    WDB[loc1 + i0 + i] = RDB[loc1 + i0 + i] + RDB[ptr + i];
		}	      

	      /* Add reaction to total list */

	      ptr = NewItem(tot + REACTION_PTR_PARTIAL_LIST, 
			    RLS_DATA_BLOCK_SIZE);
	      
	      /* Allocate memory for reaction counter */
		
	      ptp = AllocPrivateData(1, PRIVA_ARRAY);
	      WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

	      WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
	      WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;
	      WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
	      WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];

	      /* Add reaction to absorption list */

	      if ((ty == 0) && (type != NUCLIDE_TYPE_PHOTON))
		{
		  ptr = NewItem(abs + REACTION_PTR_PARTIAL_LIST, 
				RLS_DATA_BLOCK_SIZE);
		  
		  /* Allocate memory for reaction counter */
		  
		  ptp = AllocPrivateData(1, PRIVA_ARRAY);
		  WDB[ptr + RLS_DATA_PTR_COUNT] = (double)ptp;

		  WDB[ptr + RLS_DATA_PTR_REA] = (double)rea;
		  WDB[ptr + RLS_DATA_COMP_IDX] = -1.0;
		  WDB[ptr + RLS_DATA_EMIN] = RDB[rea + REACTION_EMIN];
		  WDB[ptr + RLS_DATA_EMAX] = RDB[rea + REACTION_EMAX];
		}
	    }

	  /* Next reaction */

	  rea = NextItem(rea);
	}

      /***********************************************************************/
    }

  /* Put maximum cross section */

  WDB[nuc + NUCLIDE_MAX_TOTXS] = max;

  /* Close total list */

  ptr = (long)RDB[tot + REACTION_PTR_PARTIAL_LIST];
  CloseList(ptr);

  /* Sort list by minimum energy */
  
  SortList(ptr, RLS_DATA_EMIN, SORT_MODE_ASCEND);

  /* Check if absorption list was created */

  if (abs > VALID_PTR)
    if ((ptr = (long)RDB[abs + REACTION_PTR_PARTIAL_LIST]) > VALID_PTR)
      {
	/* Close list */
	
	CloseList(ptr);
	
	/* Sort list by minimum energy */
	
	SortList(ptr, RLS_DATA_EMIN, SORT_MODE_ASCEND);
      }
}

/*****************************************************************************/
