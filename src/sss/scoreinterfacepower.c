#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : scoreinterfacepower.c                          */
/*                                                                           */
/* Created:       2012/02/16 (JLe)                                           */
/* Last modified: 2016/02/17 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Scores power for multi-physics interface                     */
/*                                                                           */
/* Comments: - Polttoaineinterfacen aksiaalijako lisätty 3.4.2013            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ScoreInterfacePower:"

/*****************************************************************************/

void ScoreInterfacePower(double fissE, double wgt, double x, double y, 
			 double z, double t, long id)
{
  long loc0, loc1, ptr, ncol, idx, nz, nr, i, j, k, uni, na, nt, l;
  double zmin, zmax, x0, y0, r0, r, f, r2, phi;

  /* Check that interfaces are defined */

  if ((loc0 = (long)RDB[DATA_PTR_IFC0]) < VALID_PTR)
    return;

  /***************************************************************************/

  /***** Interface to thermal hydraulics *************************************/

  /* Get zone index */

  ptr = (long)RDB[DATA_PTR_ZONE_IDX];
  idx = (long)GetPrivateData(ptr, id);

  /* Loop over interfaces */

  while (loc0 > VALID_PTR)
    {
      /* Check flag and type */
      
      if (((long)RDB[loc0 + IFC_CALC_OUTPUT] == YES) &&
	  (((long)RDB[loc0 + IFC_TYPE] < IFC_TYPE_FUEP) || 
	   ((long)RDB[loc0 + IFC_TYPE] > IFC_TYPE_FPIP)))
	{
	  /* Get pointer to statistics */
		    
	  ptr = (long)RDB[loc0 + IFC_PTR_STAT];
	  CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

	  /* Check type */

	  if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH)
	    {
	      /* Get collision number */
	      
	      ncol = (long)RDB[DATA_PTR_COLLISION_COUNT];
	      ncol = (long)GetPrivateData(ncol, id);
	  
	      /* Get pointer to tet cell (set in ifcpoint.c) */

	      if ((loc1 = TestValuePair(loc0 + IFC_PTR_PREV_COL_CELL, ncol, id))
		  > VALID_PTR)
		{
		  /* Get index to statistics and score */
		  
		  if ((i = (long)RDB[loc1 + IFC_TET_MSH_STAT_IDX]) > -1)
		    AddBuf1D(fissE, wgt, ptr, id, i);
		}
	    }
	  else
	    {
	      /* Other types, get axial boundaries and number of bins */
	      
	      zmin = RDB[loc0 + IFC_ZMIN];
	      zmax = RDB[loc0 + IFC_ZMAX];
	      nz = (long)RDB[loc0 + IFC_NZ];
	      nr = (long)RDB[loc0 + IFC_NR];
	      
	      /* Find region (toi z-tarkistus tarvitaan että indeksi menee */
	      /* alarajalla oikein) */

	      if ((z >= zmin) && (z < zmax))
		if ((loc1 = (long)RDB[loc0 + IFC_PTR_SCORE]) > VALID_PTR)
		  if ((loc1 = SeekList(loc1, IFC_SCORE_REG_IDX, idx, 
				       SORT_MODE_ASCEND)) > VALID_PTR)
		    {
		      /* Get stat index */
		      
		      i = (long)RDB[loc1 + IFC_SCORE_STAT_IDX];
		      CheckValue(FUNCTION_NAME, "i", "", i, 0, 
				 (long)RDB[loc0 + IFC_STAT_NREG]);
		      
		      /* Get pointer to output data */
		      
		      loc1 = (long)RDB[loc1 + IFC_SCORE_PTR_OUT];
		      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
		      
		      /* Calculate axial bin */
		      
		      j = (long)((z - zmin)/(zmax - zmin)*((double)nz));
		      CheckValue(FUNCTION_NAME, "j", "", j, 0, nz - 1);
		      
		      /* Get center coordinates and radius */
		      
		      x0 = RDB[loc1 + IFC_OUT_X0];
		      y0 = RDB[loc1 + IFC_OUT_Y0];
		      r0 = RDB[loc1 + IFC_OUT_R];
		      
		      /* Check for zero radius */
		      
		      if (r0 < ZERO)
			Die(FUNCTION_NAME, "Zero radius");		      
		      
		      /* Get square radius */
		      
		      r = (x - x0)*(x - x0) + (y - y0)*(y - y0);
		      
		      /* Calculate radial zone (typerä toi jälkimmäinen) */
		      
		      if ((f = r/r0/r0) > 1.0)
			Die(FUNCTION_NAME, 
			    "Point is not inside %E %E : %E %E : %E",
			    x, y, x0, y0, r0);
		      else if (f == 1.0)
			f = 0.999999;
		      
		      /* Index */
		      
		      k = (long)(f*nr);		      
		      CheckValue(FUNCTION_NAME, "k", "", k, 0, nr - 1);
		      
		      /* Score */
		      
		      AddBuf(fissE, wgt, ptr, id, -1, i, j, k);
		    }
	    }
	}
      
      /* Next interface */
      
      loc0 = NextItem(loc0);
    }

  /***************************************************************************/

  /***** Interface to fuel performance codes *********************************/

  /* Get collision universe */
  
  ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
  uni = GetPrivateData(ptr, id);
  
  /* Check pointer */
  
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
  
  /* Check pointer to interface */
  
  if ((loc1 = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) > VALID_PTR)
    {
      /* Get coordinates */
      
      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_X];
      CheckPointer(FUNCTION_NAME, "(xptr)", PRIVA_ARRAY, ptr);
      x = GetPrivateData(ptr, id);
      
      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Y];
      CheckPointer(FUNCTION_NAME, "(yptr)", PRIVA_ARRAY, ptr);
      y = GetPrivateData(ptr, id);

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Z];
      CheckPointer(FUNCTION_NAME, "(zptr)", PRIVA_ARRAY, ptr);
      z = GetPrivateData(ptr, id);

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_T];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      t = GetPrivateData(ptr, id);

      /* Coordinate transformation to cold system */

      CoordExpans(loc1, &x, &y, &z, t, 1);

      /* Get the size of the statistics */
      
      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
      CheckPointer(FUNCTION_NAME, "(limptr)", DATA_ARRAY, ptr);

      nz = (long)RDB[ptr + FUEP_NZ];
      na = (long)RDB[ptr + FUEP_NA];
      nr = (long)RDB[ptr + FUEP_NR];
      nt = (long)RDB[ptr + FUEP_NT];

      /* Get pointer to axial output zones */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_Z];
      CheckPointer(FUNCTION_NAME, "(zptr)", DATA_ARRAY, ptr);

      /* Find interval */	  

      i = SearchArray(&RDB[ptr], z, nz + 1);

      /* Check */

      if (i < 0)
	{
	  Warn(FUNCTION_NAME,"Outside of axial zone");
	  return;
	}

      /* Calculate angle */

      phi = PolarAngle(x, y);
	  
      /* Check phi */

      CheckValue(FUNCTION_NAME, "phi", "", phi, 0.0, 2.0*PI);

      /* Get pointer to angular output zones */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_PHI];
      CheckPointer(FUNCTION_NAME, "(aptr)", DATA_ARRAY, ptr);

      /* Rotate if needed */

      if (phi > 2.0*PI+RDB[ptr])
	phi -= 2.0*PI;
      
      /* Find interval */

      j = SearchArray(&RDB[ptr], phi, na + 1);

      /* Check */

      if (j < 0)
	{
	  Warn(FUNCTION_NAME,"Outside of angular zone");
	  return;
	}

      /* Calculate square radius */

      r2 = x*x + y*y;

      /* Get pointer to radial output zones */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_R2];
      CheckPointer(FUNCTION_NAME, "(rptr)", DATA_ARRAY, ptr);

      /* Find interval */
	  
      k = SearchArray(&RDB[ptr], r2, nr + 1);

      /* Check */

      if (k < 0)
	{
	  Warn(FUNCTION_NAME,"Outside of radial zone");
	  return;
	}

      /* Find time bin */

      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_TB];
      CheckPointer(FUNCTION_NAME, "(TBptr)", DATA_ARRAY, ptr);

      /* Find bin*/

      l = SearchArray(&RDB[ptr], t, nt + 1);

      /* Check */

      if (l < 0)
	return;

      /* Score */

      ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
      CheckPointer(FUNCTION_NAME, "(Pptr)", DATA_ARRAY, ptr);
      
      AddBuf(fissE, wgt, ptr, id, -1, i, j, k, l);
    }

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
