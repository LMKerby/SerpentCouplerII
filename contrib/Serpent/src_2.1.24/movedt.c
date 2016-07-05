/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : movedt.c                                       */
/*                                                                           */
/* Created:       2012/10/05 (JLe)                                           */
/* Last modified: 2014/08/22 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Moves particle forward using delta-tracking                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MoveDT:"

/*****************************************************************************/

long MoveDT(long part, double majorant, double minxs, long *cell, double *xs0, 
	    double *x, double *y, double *z, double *l0, double *u, double *v, 
	    double *w, double E, long id)
{
  long ptr, type, mat, mat0, bc;
  double totxs, l, wgt;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Get particle type */

  type = (long)RDB[part + PARTICLE_TYPE];  

  /* Add to DT fraction counter */
		
  ptr = (long)RDB[RES_DT_TRACK_FRAC];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  AddBuf1D(1.0, 1.0, ptr, id, 2 - type);
  
  /* Check cross section and sample path length */
  
  CheckValue(FUNCTION_NAME, "majorant", "", majorant, ZERO, INFTY);
  l = -log(RandF(id))/majorant;

  /* Move particle to collision site */
		  
  *x = *x + *u*l;
  *y = *y + *v*l;
  *z = *z + *w*l;

  /* Set path length */
      
  *l0 = l;

  /* Reset previous material pointer and total xs */

  mat0 = -1;
  totxs = 0.0;

  /* Find location */

  *cell = WhereAmI(*x, *y, *z, *u, *v, *w, id);
  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

  /* Check if point is outside the geometry */

  if ((long)RDB[*cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
    {
      /* Check if track is stopped at outer boundary */
		  
      if ((long)RDB[DATA_STOP_AT_BOUNDARY] == NO)
	{
	  /* Set weight to test that albedo boundary conditions were not */
	  /* applied */

	  wgt = 1.0;

	  /* Apply repeated boundary conditions */

	  bc = BoundaryConditions(cell, x, y, z, u, v, w, &wgt, id);

	  /* Check that condition was applied */
	  
	  if (bc != YES)
	    Die(FUNCTION_NAME, "Repeated bc not apllied");
 
	  /* Check change in weight */

	  if (wgt != 1.0)
	    Die(FUNCTION_NAME, "Change in weight (albedo)");

	  /* Find location */

	  *cell = WhereAmI(*x, *y, *z, *u, *v, *w, id);
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);
	}
      else 
	{
	  /* Stop track at boundary */

	  StopAtBoundary(cell, x, y, z, l0, *u, *v, *w, id);
      
	  /* Set cross section */

	  *xs0 = -1.0;

	  /* Return surface crossing */
      
	  return TRACK_END_SURF;
	}
    }

  /* Get material pointer */

  mat = (long)RDB[*cell + CELL_PTR_MAT];
  mat = MatPtr(mat, id);
  
  /* Check pointer */
  
  if (mat != mat0)
    {
      /* Get total cross section */
      
      totxs = TotXS(mat, type, E, id);
      
      /* Remember material */
      
      mat0 = mat;
    }
  
  /* Check total xs */
  
  CheckValue(FUNCTION_NAME, "totxs", "", totxs, 0.0, INFTY);
  
  /* Compare to minimum value */
  
  if (totxs > minxs)
    {
      /* Sample between real and virtual collision */
      
      if (RandF(id) < totxs/majorant)
	{
	  /* Set cross section */
	  
	  *xs0 = totxs;
	  
	  /* Add to success counter */
	  
	  ptr = (long)RDB[RES_DT_TRACK_EFF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
	  AddBuf1D(1.0, 1.0, ptr, id, 2 - type);
	  
	  /* Return collision */
	  
	  return TRACK_END_COLL;
	}
      else
	{
	  /* Set cross section */
	  
	  *xs0 = -1.0;
	  
	  /* Add to failure counter */
	  
	  ptr = (long)RDB[RES_DT_TRACK_EFF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
	  AddBuf1D(1.0, 1.0, ptr, id, 4 - type);
	  
	  /* Return virtual */
	  
	  return TRACK_END_VIRT;
	}
    }
  else
    {
      /* Sample scoring */
      
      if (RandF(id) < minxs/majorant)
	{
	  /* Set cross section */
	  
	  *xs0 = minxs;
	  
	  /* Sample between real and virtual collision */
	  
	  if (RandF(id) < totxs/minxs)
	    {
	      /* Add to success counter */
	      
	      ptr = (long)RDB[RES_DT_TRACK_EFF];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
	      AddBuf1D(1.0, 1.0, ptr, id, 2 - type);
	      
	      /* Return collision */
	      
	      return TRACK_END_COLL;
	    }
	  else
	    {
	      /* Add to failure counter */
	      
	      ptr = (long)RDB[RES_DT_TRACK_EFF];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
	      AddBuf1D(1.0, 1.0, ptr, id, 4 - type);
	      
	      /* Return virtual */
	      
	      return TRACK_END_VIRT;
	    }
	}
      else
	{
	  /* Set cross section */
	  
	  *xs0 = -1.0;
	  
	  /* Add to failure counter */
	  
	  ptr = (long)RDB[RES_DT_TRACK_EFF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY,ptr);
	  AddBuf1D(1.0, 1.0, ptr, id, 4 - type);
	  
	  /* Return virtual */
	  
	  return TRACK_END_VIRT;
	}
    }
}

/*****************************************************************************/
