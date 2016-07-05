/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processtfb.c                                   */
/*                                                                           */
/* Created:       2012/01/08 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Sets up processes temperature feedback data                  */
/*                                                                           */
/* Comments: - This must be called before calculating majorants              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessTFB:"

/*****************************************************************************/

void ProcessTFB()
{
  long tfb, reg, nst, mat, mat0, loc0, surf, ptr, n, m, match, cell;
  double r;
  char tmpstr[MAX_STR];

  /* Loop over feedbacks */

  tfb = (long)RDB[DATA_PTR_TFB0];
  while (tfb > VALID_PTR)
    {
      /* Find nest */

      nst = RDB[DATA_PTR_NST0];
      if ((nst = SeekListStr(nst, NEST_PTR_NAME, GetText(tfb + TFB_PTR_NST)))
	  < VALID_PTR)
	Error(tfb, "Nest %s is not defined", GetText(tfb + TFB_PTR_NST));
      
      /* Initialize the power of the nest */
  
      WDB[tfb + TFB_ITER_POW] = 0.0;

      /* Put pointer to nest */

      WDB[tfb + TFB_PTR_NST] = (long)nst;

      /* Loop over nest regions */
      
      loc0 = (long)RDB[nst + NEST_PTR_REGIONS];
      while (loc0 > VALID_PTR)
	{
	  /* Pointer to cell */
	  
	  cell = (long)RDB[loc0 + NEST_REG_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
	  
	  /* Pointer to material */
	  
	  if ((mat = (long)RDB[cell + CELL_PTR_MAT]) < VALID_PTR)
	    {
	      /* No material, get pointer to next */
	      
	      loc0 = NextItem(loc0);
	      
	      /* Cycle loop */
	      
	      continue;
	    }

	  /* Use parent material for name comparison if divided for */
	  /* burnup calculation */

	  if ((mat0 = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) < VALID_PTR)
	    mat0 = mat;
	  
	  /* Loop over region list and find match */
	  
	  reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
	  while (reg > VALID_PTR)
	    {
	      /* Initialize some values */
      
	      WDB[reg + TFB_REG_ITER_C0] = 0.0;
	      WDB[reg + TFB_REG_ITER_C1] = 0.0;
	      WDB[reg + TFB_REG_ITER_C2] = RDB[tfb + TFB_TLIM];
	      WDB[reg + TFB_REG_ITER_TEMP] = RDB[tfb + TFB_TLIM];
              WDB[reg + TFB_REG_ITER_POW] = 0.0;
              WDB[reg + TFB_REG_ITER_POWIN] = 0.0;

	      /* Reset match */

	      match = NO;

	      /* Check if region is already assigned with a material */

	      if ((long)RDB[reg + TFB_REG_IDX] > 0)
		{
		  /* Yes, get pointer to parent */

		  ptr = (long)RDB[reg + TFB_REG_PTR_MAT];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  ptr = (long)RDB[ptr + MATERIAL_DIV_PTR_PARENT];
		  
		  /* Compare pointers */
		  
		  if (((long)RDB[reg + TFB_REG_PTR_MAT] == mat) ||
		      ((ptr > VALID_PTR) && (ptr == mat0)))
		    {
		      /* Duplicate region */

		      reg = DuplicateItem(reg);

		      /* Set flag */

		      match = YES;
		    }
		}
	      else
		{
		  /* No, compare names and set flag */
		  
		  if (CompareStr(mat0 + MATERIAL_PTR_NAME, 
				 reg + TFB_REG_PTR_MAT))
		    match = YES;
		}
	      
	      /* Check match and break */

	      if (match == YES)
		break;

	      /* Next */

	      reg = NextItem(reg);
	    }

	  /* Check region pointer */

	  if (reg > VALID_PTR)
	    {
	      /* Add counter */
	      
	      WDB[reg + TFB_REG_IDX] = WDB[reg + TFB_REG_IDX] + 1.0;
	      
	      /* Put pointer to material */
	      
	      WDB[reg + TFB_REG_PTR_MAT] = (double)mat;
		  
	      /* Put pointer to feedback region */
	      
	      if ((long)RDB[reg + NEST_REG_PTR_TFB_REG] > VALID_PTR)
		Error(tfb, "Multiple feedbaks for nest %s region %s", 
		      GetText(nst + NEST_PTR_NAME), 
		      GetText(mat + MATERIAL_PTR_NAME)); 
	      else
		WDB[loc0 + NEST_REG_PTR_TFB_REG] = (double)reg;
	      
	      /* Pointer to inside surface */
	      
	      if ((surf = (long)RDB[loc0 + NEST_REG_PTR_SURF_IN]) > VALID_PTR)
		{
		  /* Check type */
		  
		  if (((long)RDB[surf + SURFACE_TYPE] != SURF_CYL) &&
		      ((long)RDB[surf + SURFACE_TYPE] != SURF_CYLZ))
		    Error(tfb, "Invalid surface type in nest %s",
			  GetText(nst + NEST_PTR_NAME));
		  
		  /* Pointer to parameter list */
		  
		  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  
		  /* Put radius */
		  
		  WDB[reg + TFB_REG_R1] = RDB[ptr + 2];
		  WDB[reg + TFB_REG_ITER_R1] = RDB[ptr + 2];

		  /* Store pointer to radius */

		  WDB[reg + TFB_REG_PTR_RAD_IN] = (double)(ptr + 2);
		}
	      
	      /* Pointer to outside surface */
	      
	      if ((surf = (long)RDB[loc0 + NEST_REG_PTR_SURF_OUT]) > VALID_PTR)
		{
		  /* Check type */
		  
		  if (((long)RDB[surf + SURFACE_TYPE] != SURF_CYL) &&
		      ((long)RDB[surf + SURFACE_TYPE] != SURF_CYLZ))
		    Error(tfb, "Invalid surface type in nest %s",
			  GetText(nst + NEST_PTR_NAME));
		  
		  /* Pointer to parameter list */
		  
		  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  
		  /* Put radius */
		  
		  WDB[reg + TFB_REG_R0] = RDB[ptr + 2];
		  WDB[reg + TFB_REG_ITER_R0] = RDB[ptr + 2];

		  /* Store pointer to radius */

		  WDB[reg + TFB_REG_PTR_RAD_OUT] = (double)(ptr + 2);
		}	  
	    }
	  
	  /* Next nest region */

	  loc0 = NextItem(loc0);
	}

      /* Check that materials are found and reset indexes */

      reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
      while (reg > VALID_PTR)
	{
	  if ((long)RDB[reg + TFB_REG_IDX] > 0)
	    WDB[reg + TFB_REG_IDX] = 0.0;
	  else
	    Error(tfb, "Material %s not found in nest %s", 
		  GetText(reg + TFB_REG_PTR_MAT), 
		  GetText(nst + NEST_PTR_NAME));
	  
	  reg = NextItem(reg);
	}

      /* Sort list */

      reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
      SortList(reg, TFB_REG_R0, SORT_MODE_ASCEND);

      /* Set size */

      WDB[tfb + TFB_N_REG] = (double)ListSize(reg);

      /* Check radii, put indexes and get boundary condition */

      r = -1.0;
      n = 0;

      reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
      while (reg > VALID_PTR)
	{
	  /* Put index */

	  WDB[reg + TFB_REG_IDX] = (double)(n++);

	  /* Pointer to material */

	  mat = (long)RDB[reg + TFB_REG_PTR_MAT];
	  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
	  
	  /* Put material minimum temperature (used for pre-broadeing) */

	  if (((long)RDB[mat + MATERIAL_TMS_MODE] == NO) || 
	      (RDB[reg + TFB_REG_TMS_TMIN] < RDB[mat + MATERIAL_TMS_TMIN]))
	    {
	      /* Put temperature to material */

	      WDB[mat + MATERIAL_TMS_TMIN] = RDB[reg + TFB_REG_TMS_TMIN];

	      /* Put temperature to parent */

	      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
		WDB[ptr + MATERIAL_TMS_TMIN] = RDB[reg + TFB_REG_TMS_TMIN];
	    }

	  /* Put material maximum temperature (used for majorant) */

	  if (((long)RDB[mat + MATERIAL_TMS_MODE] == NO) || 
	      (RDB[reg + TFB_REG_TMS_TMAX] > RDB[mat + MATERIAL_TMS_TMAX]))
	    {
	      /* Put temperature to material */

	      WDB[mat + MATERIAL_TMS_TMAX] = RDB[reg + TFB_REG_TMS_TMAX];

	      /* Put temperature to parent */

	      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
		WDB[ptr + MATERIAL_TMS_TMAX] = RDB[reg + TFB_REG_TMS_TMAX];
	    }
	  
	  /* Put TMS flag */

	  WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

	  if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
	    WDB[ptr + MATERIAL_TMS_MODE] = (double)YES;

	  /* Check radii */

	  if ((r > 0.0) && (r != RDB[reg + TFB_REG_R0]))
	    Error(tfb, "Discontinuity before material %s",
		  GetText(mat + MATERIAL_PTR_NAME));
	  else
	    r = RDB[reg + TFB_REG_R1];

	  /* Next region */

	  reg = NextItem(reg);
	}

      /* Get number of regions and time bins */

      n = (long)RDB[tfb + TFB_N_REG];
      m = (long)RDB[DATA_DYN_NB];

      /* Allocate memory for mean power */
      
      sprintf(tmpstr, "TFB_POWER_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, n, m);
      WDB[tfb + TFB_PTR_MEAN_POW] = (double)ptr;  

      /* Allocate memory for volume-averaged temperature */
      
      sprintf(tmpstr, "TFB_VTEMP_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, n, m);
      WDB[tfb + TFB_PTR_MEAN_VTEMP] = (double)ptr;  

      /* Allocate memory for maximum temperature */
      
      sprintf(tmpstr, "TFB_MAX_TEMP_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, n, m);
      WDB[tfb + TFB_PTR_MAX_TEMP] = (double)ptr;  

      /* Allocate memory for minimum temperature */
      
      sprintf(tmpstr, "TFB_MIN_TEMP_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, n, m);
      WDB[tfb + TFB_PTR_MIN_TEMP] = (double)ptr;  

      /* Allocate memory for flux-averaged temperature */
      
      sprintf(tmpstr, "TFB_FTEMP_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, 2*n, m);
      WDB[tfb + TFB_PTR_MEAN_FTEMP] = (double)ptr;  

      /* Allocate memory for densities */
      
      sprintf(tmpstr, "TFB_MDENS_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, n, m);
      WDB[tfb + TFB_PTR_MEAN_MDENS] = (double)ptr;  

      /* Allocate memory for radii */
      
      sprintf(tmpstr, "TFB_RAD_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, n + 1, m);
      WDB[tfb + TFB_PTR_MEAN_RAD] = (double)ptr;  

      /* Allocate memory for heat conductivities (onko toi n + 1 typo?) */
      
      sprintf(tmpstr, "TFB_HC_%s", GetText(nst + NEST_PTR_NAME));

      ptr = NewStat(tmpstr, 2, n + 1, m);
      WDB[tfb + TFB_PTR_MEAN_HC] = (double)ptr;  

      /* Next feedback */
      
      tfb = NextItem(tfb);
    }
}

/*****************************************************************************/
