/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findmaterialpointers.c                         */
/*                                                                           */
/* Created:       2010/10/06 (JLe)                                           */
/* Last modified: 2015/05/27 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Finds material pointers for cells, nests and detectors, etc. */
/*                                                                           */
/* Comments: - Toi repro-flagi resetoidaan että fyysiset repro-materiaalit   */
/*             tulee majoranttiin mukaan.                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindMaterialPointers:"

/*****************************************************************************/

void FindMaterialPointers()
{
  long cell, mix, mat, ptr, uni, umsh, loc0;

  /* Check that materials are defined */

  if ((long)RDB[DATA_PTR_M0] < VALID_PTR)
    Error(0, "No material definitions in geometry");

  fprintf(out, "Linking materials to geometry...\n");

  /**************************************************************************/

  /***** Cells **************************************************************/

  /* Loop over cells */

  cell = RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {
      /* Check if material is defined */
      
      if ((long)RDB[cell + CELL_PTR_MAT] < VALID_PTR)
	{
	  /* Check fill pointer */

	  if ((long)RDB[cell + CELL_PTR_FILL] < VALID_PTR)
	    Die(FUNCTION_NAME, "Pointer error");

	  /* Set type */
	  
	  WDB[cell + CELL_TYPE] = (double)CELL_TYPE_FILL;

	  /* Put null pointer */
	  
	  WDB[cell + CELL_PTR_MAT] = NULLPTR;
	}
      else if (!strcmp(GetText(cell + CELL_PTR_MAT), "void"))
	{
	  /* Void cell, put null pointer */
	  
	  WDB[cell + CELL_PTR_MAT] = NULLPTR;

	  /* Set type */
	  
	  WDB[cell + CELL_TYPE] = CELL_TYPE_VOID;
	}
      else if (!strcmp(GetText(cell + CELL_PTR_MAT), "outside"))
	{
	  /* Outside cell, put null pointer */
	  
	  WDB[cell + CELL_PTR_MAT] = NULLPTR;

	  /* Set type */
	  
	  WDB[cell + CELL_TYPE] = CELL_TYPE_OUTSIDE;

	  /* Check union type */

	  if ((long)RDB[cell + CELL_PTR_SURF_COMP] > VALID_PTR)
	    Error(cell, "Unions are not allowed for outside cell %s",
		  GetText(cell + CELL_PTR_NAME));
	}

      else if ((mat = (long)RDB[cell + CELL_PTR_REG_MAT]) > VALID_PTR)
	{
	  /* Pointer set when dividing a nest */

	  WDB[cell + CELL_PTR_MAT] = (double)mat;

	  /* Set type */

	  WDB[cell + CELL_TYPE] = CELL_TYPE_MAT;

	  /* Set used-, physical- and majorant-flags */
	  
	  SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
	  SetOption(mat + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);
	  SetOption(mat + MATERIAL_OPTIONS, OPT_INCLUDE_MAJORANT);
	}
      else
	{
	  /* Set type */
	  
	  WDB[cell + CELL_TYPE] = CELL_TYPE_MAT;

	  /* Find material */

	  mat = RDB[DATA_PTR_M0];
	  if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, 
				 GetText(cell + CELL_PTR_MAT))) > VALID_PTR)
	    {
	      /* Put pointer */
		  
	      WDB[cell + CELL_PTR_MAT] = (double)mat;
	      
	      /* Set used-, physical- and majorant-flags */
	      
	      SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
	      SetOption(mat + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);
	      SetOption(mat + MATERIAL_OPTIONS, OPT_INCLUDE_MAJORANT);
	    }
	  else
	    Error(cell, "Material %s is not defined", 
		  GetText(cell + CELL_PTR_MAT));
	}

      /* Next */

      cell = NextItem(cell);
    }

  /**************************************************************************/

  /***** Unstructured mesh based geometries *********************************/

  /* Loop over structures (pointers are set in readumshgeometry.c) */

  umsh = (long)RDB[DATA_PTR_UMSH0];
  while (umsh > VALID_PTR)
    {
      /* Pointer to interface structure */
      
      loc0 = (long)RDB[umsh + UMSH_PTR_IFC];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Pointer to tet cells */

      loc0 = (long)RDB[loc0 + IFC_PTR_TET_MSH];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Loop over cells */
      
      while (loc0 > VALID_PTR)
	{
	  /* Pointer to geometry cell */
	  
	  cell = (long)RDB[loc0 + IFC_TET_MSH_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	  /* Get pointer to material */

	  mat = (long)RDB[cell + CELL_PTR_MAT];
	  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
	  
	  /* Set used-, physical- and majorant-flags */
	      
	  SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
	  SetOption(mat + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);
	  SetOption(mat + MATERIAL_OPTIONS, OPT_INCLUDE_MAJORANT);
	  
	  /* Next cell */

	  loc0 = NextItem(loc0);
	}
      
      /* Next */

      umsh = NextItem(umsh);
    }

  /**************************************************************************/

  /***** Mixtures ***********************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check used-flag */
      
      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_USED)
	{
	  /* Loop over mixture */

	  mix = (long)RDB[mat + MATERIAL_PTR_MIX];
	  while (mix > VALID_PTR)
	    {
	      /* Find material */
	      
	      ptr = (long)RDB[DATA_PTR_M0];
	      if ((ptr = SeekListStr(ptr, MATERIAL_PTR_NAME, 
				     GetText(mix + MIXTURE_PTR_MAT))) 
		  > VALID_PTR)
		{
		  /* Put pointer */
		  
		  WDB[mix + MIXTURE_PTR_MAT] = (double)ptr;
		  
		  /* Set used-flags */
		  
		  SetOption(ptr + MATERIAL_OPTIONS, OPT_USED);
		}
	      else
		Error(mat, "Material %s is not defined", 
		      GetText(mix + MIXTURE_PTR_MAT));
	      
	      /* Next */
	      
	      mix = NextItem(mix);
	    }
	}
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  /**************************************************************************/

  /***** DT block and force lists *******************************************/

  /* Loop over blocks */

  if ((ptr = (long)RDB[DATA_DT_PTR_BLOCK_LIST]) > VALID_PTR)
    while ((long)RDB[ptr] > VALID_PTR)
      {
	/* Find match */

	mat = (long)RDB[DATA_PTR_M0];
	while (mat > VALID_PTR)
	  {
	    /* Compare name */

	    if (CompareStr(ptr, mat + MATERIAL_PTR_NAME))
	      {
		/* Put option */

		WDB[mat + MATERIAL_DT_MODE] = (double)DT_MAT_BLOCK;

		/* Break loop */

		break;
	      }

	    /* Next material */
	    
	    mat = NextItem(mat);
	  }

	/* Check match */

	if (mat < VALID_PTR)
	  Note(0, "Material %s not found for blocked DT mode", GetText(ptr));

	/* Next */

	ptr++;
      }

  /* Loop over forces */

  if ((ptr = (long)RDB[DATA_DT_PTR_FORCE_LIST]) > VALID_PTR)
    while ((long)RDB[ptr] > VALID_PTR)
      {
	/* Find match */

	mat = (long)RDB[DATA_PTR_M0];
	while (mat > VALID_PTR)
	  {
	    /* Compare name */

	    if (CompareStr(ptr, mat + MATERIAL_PTR_NAME))
	      {
		/* Put option */

		WDB[mat + MATERIAL_DT_MODE] = (double)DT_MAT_FORCE;

		/* Break loop */

		break;
	      }

	    /* Next material */
	    
	    mat = NextItem(mat);
	  }

	/* Check match */

	if (mat < VALID_PTR)
	  Note(0, "Material %s not found for forced DT mode", GetText(ptr));

	/* Next */

	ptr++;
      }
  
  /**************************************************************************/

  /***** Normalization ******************************************************/

  /* NOTE: This will check and set pointers to parent materials */
  /*       pointers to divided are set in setnormalization.c    */

  ptr = (long)RDB[DATA_PTR_NORM];
  while (ptr > VALID_PTR)
    {
      /* Check if material is given */

      if ((long)RDB[ptr + NORM_PTR_MAT] > VALID_PTR)
	{
	  /* Loop over materials to find match */

	  mat = (long)RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Compare */

	      if (CompareStr(mat + MATERIAL_PTR_NAME, ptr + NORM_PTR_MAT))
		break;

	      /* Next material */

	      mat = NextItem(mat);
	    }
	  
	  /* Check pointer */

	  if (mat < VALID_PTR)
	    Error(0, "Material %s used for normalization is not defined",
		  GetText(ptr + NORM_PTR_MAT));

	  /* Put pointers */

	  WDB[ptr + NORM_PTR_MAT] = (double)mat;
	  WDB[mat + MATERIAL_PTR_NORM] = (double)ptr;
	}

      /* Next */

      ptr = NextItem(ptr);
    }

  /**************************************************************************/

  /***** Set some options on divided materials ******************************/

  /* NOTE: These flags apply to Serpent 1 -type division. Flags */
  /* for the new type division are set in dividezone.c. */

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get pointer to parent and set flags */
	  
      if (((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR) &&
	  ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_PHYSICAL_MAT))	
	{
	  /* Used-flag */

	  SetOption(ptr + MATERIAL_OPTIONS, OPT_USED);

	  /* Reset physical and majorant flags (this is probably redundant) */

	  ResetOption(ptr + MATERIAL_OPTIONS, OPT_PHYSICAL_MAT);
	  ResetOption(ptr + MATERIAL_OPTIONS, OPT_INCLUDE_MAJORANT);
	}

      /* Next material */

      mat = NextItem(mat);
    }

  /**************************************************************************/

  /***** Check that branch-replaced materials are not used ******************/

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check used and replaced flags */

      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_USED) &&
	  ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_REPLACED_MAT))
	Error(mat, 
	      "Material %s used in coefficient calculation is part of geometry",
	      GetText(mat + MATERIAL_PTR_NAME));

      /* Next material */

      mat = NextItem(mat);
    }

  /**************************************************************************/

  /***** Check that outside world exists *************************************/

  /* Get pointer to root universe */

  uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Check universe type */
  
  if ((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_NEST)
    {
      /* Pointer to nest */
	
      ptr = (long)RDB[uni + UNIVERSE_PTR_NEST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	
      /* Get pointer to regions */
      
      ptr = (long)RDB[ptr + NEST_PTR_REGIONS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Loop over regions */

      while (ptr > VALID_PTR)
	{
	  /* Pointer to cell */
	  
	  cell = (long)RDB[ptr + NEST_REG_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	  /* Check type */

	  if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
	    break;

	  /* Next region */

	  ptr = NextItem(ptr);
	}

      /* Check pointer */

      if (ptr < VALID_PTR)
	Note(0, "Root universe %s does not contain the \"outside world\"", 
	      GetText(uni + UNIVERSE_PTR_NAME));
    }
  else if ((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_CELL)
    {
      /* Pointer to cell list */
	
      ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	
      /* Loop over cell list */
	
      while (ptr > VALID_PTR)
	{
	  /* Pointer to cell */
	  
	  cell = (long)RDB[ptr + CELL_LIST_PTR_CELL];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	  /* Check type */

	  if ((long)RDB[cell + CELL_TYPE] == CELL_TYPE_OUTSIDE)
	    break;

	  /* Next region */
	  
	  ptr = NextItem(ptr);
	}
      
      /* Check pointer */
      
      if (ptr < VALID_PTR)
	Note(0, "Root universe %s does not contain the \"outside world\"", 
	     GetText(uni + UNIVERSE_PTR_NAME));
    }

  fprintf(out, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
