/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : whereami.c                                     */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2015/07/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Finds neutron location in universes                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WhereAmI:"

/*****************************************************************************/

long WhereAmI(double x, double y, double z, double u, double v, double w, 
	      long id)
{
  long uni, lvl0, lvl, cell, mat, nst, reg, lat, ptr, pbd, umsh, pbl, stl;
  long ncol, zone, idx, type, cgns, tra;

  /* Get pointer to root universe */
  
  uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

  /* Pointer to first level */

  lvl0 = (long)RDB[DATA_PTR_LVL0];
  CheckPointer(FUNCTION_NAME, "(lvl0)", DATA_ARRAY, lvl0);

  /* Get collision number */

  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Reset zone index */

  zone = 0;

  /* Loop over levels */

  while (lvl0 > VALID_PTR)
    {
      /* Pointer to private data */

      lvl = (long)RDB[lvl0 + LVL_PTR_PRIVATE_DATA];
      CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

      /* Check universe pointer */
  
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Do coordinate transformation */
      
      if ((tra = (long)RDB[uni + UNIVERSE_PTR_TRANS]) > VALID_PTR)
	CoordTrans(tra, &x, &y, &z, &u, &v, &w, id);

      /* Universe symmetry */

      if ((ptr = (long)RDB[uni + UNIVERSE_PTR_SYM]) > VALID_PTR)
	UniSym(ptr, &x, &y, &z, &u, &v, &w);

      /* Reset pointers */

      PutPrivateData(lvl + LVL_PRIV_PTR_NEST_REG, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_PTR_LAT, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_PTR_CELL, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_PTR_PBED, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_PTR_PEBBLE, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_PTR_UMSH, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_PTR_STL, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_PTR_MAT, NULLPTR, id);
      PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, NULLPTR, id);

      /* Put coordinates and direction cosines */
      
      PutPrivateData(lvl + LVL_PRIV_X, x, id);
      PutPrivateData(lvl + LVL_PRIV_Y, y, id);
      PutPrivateData(lvl + LVL_PRIV_Z, z, id);
      PutPrivateData(lvl + LVL_PRIV_U, u, id);
      PutPrivateData(lvl + LVL_PRIV_V, v, id);
      PutPrivateData(lvl + LVL_PRIV_W, w, id);

      /* Put coordinates to universe structure */

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_X];
      PutPrivateData(ptr, x, id);

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Y];
      PutPrivateData(ptr, y, id);

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Z];
      PutPrivateData(ptr, z, id);

      /* Reset last flag */
      
      PutPrivateData(lvl + LVL_PRIV_LAST, NO, id);
      
      /* Put level type */

      type = (long)RDB[uni + UNIVERSE_TYPE];
      PutPrivateData(lvl + LVL_PRIV_TYPE, type, id);
      
      /* Put universe pointer */

      PutPrivateData(lvl + LVL_PRIV_PTR_UNIV, uni, id);

      /* Put collision flag */

      StoreValuePair(uni + UNIVERSE_COL_COUNT, ncol, 1.0, id);

      /* Put gcu pointer */

      if ((long)RDB[DATA_OPTI_GC_CALC] == YES)
      	if ((ptr = (long)RDB[uni + UNIVERSE_PTR_GCU]) > VALID_PTR)
	  StoreValuePair(DATA_GCU_PTR_UNI, ncol, ptr, id);

      /* Put collision universe */

      ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
      PutPrivateData(ptr, uni, id);      

      /* Avoid warning messages */
      
      cell = -1;

      /* Check universe type */
      
      switch (type)
	{
	case UNIVERSE_TYPE_NEST:
	  {
	    /*****************************************************************/
	    
	    /***** Nest universe *********************************************/

	    /* Pointer to nest */
	    
	    nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
	    CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

	    /* Find nest region */
	    
	    reg = FindNestRegion(uni, nst, x, y, z, id);
	    CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

	    /* Get region index and update zone */

	    idx = (long)RDB[reg + NEST_REG_IDX];
	    zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
	    PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);

	    /* Put region pointer */
	    
	    PutPrivateData(lvl + LVL_PRIV_PTR_NEST_REG, reg, id);
	    
	    /* Get cell pointer */
	    
	    cell = (long)RDB[reg + NEST_REG_PTR_CELL];
	    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	    /* Put cell pointer */
	    
	    PutPrivateData(lvl + LVL_PRIV_PTR_CELL, cell, id);
	    
	    /* Put collision flag */
	    
	    StoreValuePair(cell + CELL_COL_COUNT, ncol, 1.0, id);

	    /* Put collision region */

	    StoreValuePair(nst + NEST_PTR_COL_REG, ncol, reg, id);

	    /* Check fill pointer */
	    
	    if ((ptr = RDB[reg + NEST_REG_PTR_FILL]) > VALID_PTR)
	      {
		/* Filled region, update universe pointer */
		
		uni = ptr;
	      }
	    else
	      {
		/* Put zone index */

		ptr = (long)RDB[DATA_PTR_ZONE_IDX];
		PutPrivateData(ptr, zone, id);

		/* Get material pointer */
	    
		mat = (long)RDB[cell + CELL_PTR_MAT];
		mat = MatPtr(mat, id);
		
		/* Put private pointer */
		
		PutPrivateData(lvl + LVL_PRIV_PTR_MAT, mat, id);
	    
		/* Put last flag */
		
		PutPrivateData(lvl + LVL_PRIV_LAST, YES, id);
		
		/* Return cell pointer */
		
		return cell;
	      }
	    
	    /* Break case */
	    
	    break;
	    
	    /*****************************************************************/
	  }      
	  
	case UNIVERSE_TYPE_CELL:
	  {
	    /*****************************************************************/

	    /***** Cell universe *********************************************/
	    
	    /* Find cell */

	    if ((cell = FindUniverseCell(uni, x, y, z, &idx, id)) < VALID_PTR)
	      {
		/* Geometry error, put last flag */
		
		PutPrivateData(lvl + LVL_PRIV_LAST, YES, id);

		/* Check void and plotter modes */

		if ((long)RDB[DATA_IGNORE_VOID_CELLS] == YES)
		  return (long)RDB[DATA_PTR_VOID_CELL];
		else if ((long)RDB[DATA_PLOTTER_MODE] == YES)
		  return cell;
		else
		  TrackingError(TRACK_ERR_CELL_SEARCH, -1, -1, -1, id);
	      }

	    /* Check pointer */

	    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	    /* Update zone */

	    zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
	    PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);

	    /* Put cell pointer */

	    PutPrivateData(lvl + LVL_PRIV_PTR_CELL, cell, id);

	    /* Put collision flag */

	    StoreValuePair(cell + CELL_COL_COUNT, ncol, 1.0, id);
	    
	    /* Check fill pointer and call recursively */
	    
	    if ((ptr = RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
	      {
		/* Do coordinate transformation */
      
		if ((tra = (long)RDB[cell + CELL_PTR_TRANS]) > VALID_PTR)
		  CoordTrans(tra, &x, &y, &z, &u, &v, &w, id);

		/* Filled cell, update universe pointer */
		
		uni = ptr;
	      }
	    else
	      {
		/* Put zone index */

		ptr = (long)RDB[DATA_PTR_ZONE_IDX];
		PutPrivateData(ptr, zone, id);

		/* Get material pointer */
	    
		mat = (long)RDB[cell + CELL_PTR_MAT];
		mat = MatPtr(mat, id);
		
		/* Put private pointer */
		
		PutPrivateData(lvl + LVL_PRIV_PTR_MAT, mat, id);
	    
		/* Put last flag */
		
		PutPrivateData(lvl + LVL_PRIV_LAST, YES, id);

		/* Return cell pointer */
		
		return cell;
	      }

	    /* Break case */

	    break;

	    /*****************************************************************/
	  }      

	case UNIVERSE_TYPE_LATTICE:
	  {
	    /*****************************************************************/

	    /***** Lattice universe ******************************************/

	    /* Pointer to lattice */

	    lat = (long)RDB[uni + UNIVERSE_PTR_LAT];
	    CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);
	    
	    /* Put lattice pointer */
	    
	    PutPrivateData(lvl + LVL_PRIV_PTR_LAT, lat, id);
	    
	    /* Find lattice universe */
	    
	    if ((ptr = FindLatticeRegion(lat, lvl, &x, &y, &z, &idx, id)) 
		< VALID_PTR)
	      {
		/* Put last flag */
		
		PutPrivateData(lvl + LVL_PRIV_LAST, YES, id);

		/* Check plotter mode */
		
		if ((long)RDB[DATA_PLOTTER_MODE] == YES)
		  return GEOM_ERROR_NO_CELL;
		else
		  TrackingError(TRACK_ERR_LATTICE, -1, -1, -1, id);
	      }

	    /* Update zone */

	    zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
	    PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);

	    /* Update universe pointer */
		
	    uni = ptr;

	    /* Break case */
	    
	    break;
	    
	    /*****************************************************************/
	  }
	  
	case UNIVERSE_TYPE_PBED:
	  {
	    /*****************************************************************/
	    
	    /***** Explicit stochastic geometry ******************************/

	    /* Pointer to geometry */
	    
	    pbd = (long)RDB[uni + UNIVERSE_PTR_PBED];
	    CheckPointer(FUNCTION_NAME, "(pbd)", DATA_ARRAY, pbd);
	    
	    /* Put level pointer */
	    
	    PutPrivateData(lvl + LVL_PRIV_PTR_PBED, pbd, id);
	    
	    /* Find universe */

	    ptr = FindPBRegion(uni, pbd, &x, &y, &z, &pbl, &idx, id);

	    /* Update zone (JLe: tää korjattiin 4.6.2013 / 2.1.15) */

	    zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
	    PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);

	    /* Put direct pointer to pebble */
	    
	    PutPrivateData(lvl + LVL_PRIV_PTR_PEBBLE, pbl, id);

	    /* Update universe pointer */
		
	    uni = ptr;

	    /* Break case */

	    break;
	    
	    /*****************************************************************/
	  }
	case UNIVERSE_TYPE_UMSH:
	  {
	    /*****************************************************************/
	    
	    /***** Unstructured mesh based geometry **************************/

	    /* Pointer to geometry */
	    
	    umsh = (long)RDB[uni + UNIVERSE_PTR_UMSH];
	    CheckPointer(FUNCTION_NAME, "(umsh)", DATA_ARRAY, umsh);
	    
	    /* Put level pointer */
	    
	    PutPrivateData(lvl + LVL_PRIV_PTR_UMSH, umsh, id);
	    
	    /* Reset cell */

	    cell = -1;

	    /* Check if cell pointer was set by NearestBoundary() */

	    if ((cgns = TestValuePair(uni + UNIVERSE_PTR_NEXT_CELL, ncol, id))
		> VALID_PTR)
	      {
		/* For some reason this does not work with 100% certainty, */
		/* maybe because the neutron is actually moved forward by  */
		/* extrapolation distance, and may end up in a different   */
		/* cell. */

		/* Get pointer to interface */
		ptr = (long)RDB[umsh + UMSH_PTR_IFC];

		if (InTetCell(ptr, cgns, x, y, z, YES, id) == YES)
		  {
		    /* Get cell from tet */

		    cell = (long)RDB[cgns + IFC_TET_MSH_PTR_CELL];
		    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

		  }

	      }

	    /* Check cell pointer */

	    if (cell < VALID_PTR)
	      {
		/* Get pointer to interface structure */

		ptr = (long)RDB[umsh + UMSH_PTR_IFC];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

		/* Find cell */

		if ((cell = FindTetCell(ptr, x, y, z, id)) > VALID_PTR)
		  {
		    /* Put collision flag */

		    StoreValuePair(ptr + IFC_PTR_PREV_COL_CELL, ncol, cell, 
				   id);

		    /* Store pointer to tet-cell */

		    cgns = cell;

		    /* Get pointer to geometry cell */

		    cell = (long)RDB[cell + IFC_TET_MSH_PTR_CELL];
		    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
		  }
	      }

	    /* Check cell pointer */

	    if (cell < VALID_PTR)
	      {
		/* Get pointer to background universe */

		uni = (long)RDB[umsh + UMSH_PTR_BG_UNIV];
		CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

		/* Index (nää pitää kattoa vielä myöhemmin) */

		idx = 0;

		/* Update zone */

		zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
		PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);
	      }
	    else
	      {
		/* Get index */

		idx = (long)RDB[cell + CELL_UMSH_IDX];

		/* Update zone */

		zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
		PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);

		/* Put cell pointer */

		PutPrivateData(lvl + LVL_PRIV_PTR_CELL, cell, id);

		/* Put collision flag (tätä ei tarvita mihinkään?) */

		StoreValuePair(cell + CELL_COL_COUNT, ncol, 1.0, id);

		/* Check fill pointer */
	    
		if ((ptr = RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
		  Die(FUNCTION_NAME, "Fill pointer is not null");
		
		/* Get material pointer */
		    
		mat = (long)RDB[cell + CELL_PTR_MAT];
		CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

		mat = MatPtr(mat, id);
		CheckPointer(FUNCTION_NAME, "(mat2)", DATA_ARRAY, mat);

		/* Put private pointer */
		
		PutPrivateData(lvl + LVL_PRIV_PTR_MAT, mat, id);
		
		/* Put last flag */
		
		PutPrivateData(lvl + LVL_PRIV_LAST, YES, id);
		
		/* Save tet, where collision happened */
		
		ptr = (long)RDB[cell + CELL_PTR_PREV_TET];
		CheckPointer(FUNCTION_NAME, "(ptr/tet)", PRIVA_ARRAY, ptr);

		/* Store tet cell */

		PutPrivateData(ptr, cgns, id);

		/* Return cell pointer */
		
		return cell;
	       }
	    
	    /* Break case */

	    break;
	    
	    /*****************************************************************/
	  }
	case UNIVERSE_TYPE_STL:
	  {
	    /*****************************************************************/
	    
	    /***** STL based geometry ****************************************/
	    
	    /* Pointer to geometry */
	    
	    stl = (long)RDB[uni + UNIVERSE_PTR_STL];
	    CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);
	    
	    /* Put level pointer */
	    
	    PutPrivateData(lvl + LVL_PRIV_PTR_STL, stl, id);
	    
	    /* Find solid */

	    if ((ptr = FindSTLSolid(stl, x, y, z, u, v, w, YES, id)) > 
		VALID_PTR)
	      {
		/* Get pointer to cell */

		cell = (long)RDB[ptr + STL_SOLID_PTR_CELL];
		CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

		/* Get index */

		idx = (long)RDB[ptr + STL_SOLID_REG_IDX];

		/* Update zone */

		zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
		PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);

		/* Put cell pointer */

		PutPrivateData(lvl + LVL_PRIV_PTR_CELL, cell, id);

		/* Put collision flag */

		StoreValuePair(cell + CELL_COL_COUNT, ncol, 1.0, id);

		/* Check fill pointer */
	    
		if ((ptr = RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
		  Die(FUNCTION_NAME, "Fill pointer is not null (muuta tää)");
		
		/* Get material pointer */
		    
		mat = (long)RDB[cell + CELL_PTR_MAT];
		mat = MatPtr(mat, id);

		/* Put private pointer */
		
		PutPrivateData(lvl + LVL_PRIV_PTR_MAT, mat, id);
		
		/* Put last flag */
		
		PutPrivateData(lvl + LVL_PRIV_LAST, YES, id);
		
		/* Return cell pointer */
		
		return cell;
	      }
	    else if (ptr == GEOM_ERROR_MULTIPLE_CELLS)
	      {
		/* Check plotter mode */

		if ((long)RDB[DATA_PLOTTER_MODE] == YES)
		  return GEOM_ERROR_MULTIPLE_CELLS;
		else
		  Error(stl, "Overlapping solids at [%1.2E, %1.2E, %1.2E]",
			x, y, z);
	      }
	    else
	      {
		/* Get pointer to background universe */

		uni = (long)RDB[stl + STL_PTR_BG_UNIV];
		CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

		/* Index (nää pitää kattoa vielä myöhemmin) */

		idx = 0;

		/* Update zone */

		zone = zone + ((long)RDB[lvl0 + LVL_ZONE_IDX_MULT])*idx;
		PutPrivateData(lvl + LVL_PRIV_ZONE_IDX, zone, id);
	      }
	    
	    /* Break case */

	    break;
	    
	    /*****************************************************************/
	  }

	default:
	  {
	    /* Invalid type */
	    
	    Die(FUNCTION_NAME, "Invalid universe type");
	  }
	}

      /* Next level */

      lvl0 = NextItem(lvl0);
    }

  /* Error */

  Die(FUNCTION_NAME, "No material cell");

  /* Avoid warning message */

  return 0;
}

/*****************************************************************************/
