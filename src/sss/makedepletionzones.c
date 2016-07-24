#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : makedepletionzones.c                           */
/*                                                                           */
/* Created:       2012/05/10 (JLe)                                           */
/* Last modified: 2016/01/20 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Recursive algorithm that divides materials into depletion    */
/*              zones                                                        */
/*                                                                           */
/* Comments: - Serpent 1 -type division is done at divideburnmat.c           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MakeDepletionZones:"

/*****************************************************************************/

void MakeDepletionZones(long loc0, long  uni0, long lvl, long recu, long idx0,
			long *idx1)
{
  long nst, ptr, reg, cell, lat, lst, pbd, pbl, stl, sld, i, j, mat, div, uni;
  long n, n0, n1, n2, n3, n4, n5, tot;

  /* Check pointer */

  if ((long)RDB[DATA_PTR_DIV0] < VALID_PTR)
    return;

  /* Update level pointer or get pointer to first level and universe */

  if (lvl < VALID_PTR)
    {
      lvl = (long)RDB[DATA_PTR_LVL0];
      uni0 = (long)RDB[DATA_PTR_ROOT_UNIVERSE];

      /* Print */

      fprintf(out, "Dividing materials into depletion zones...\n\n");
    }
  else
    lvl = NextItem(lvl);

  /* Get level number */

  j = (long)RDB[lvl + LVL_NUMBER];

  /* Check level and universe pointers */

  CheckPointer(FUNCTION_NAME, "(lvl)", DATA_ARRAY, lvl);
  CheckPointer(FUNCTION_NAME, "(uni0)", DATA_ARRAY, uni0);

  /* Check infinite loop */

  if (recu++ > 1000)
    Die(FUNCTION_NAME, "Infinite geometry loop involving universe %s",
	GetText(uni0 + UNIVERSE_PTR_NAME));

  /***************************************************************************/

  /***** Divide materials ****************************************************/
 
  /* Check universe type */
  
  switch((long)RDB[uni0 + UNIVERSE_TYPE])
    {
    case UNIVERSE_TYPE_NEST:
      {
	/*********************************************************************/

	/***** Nest universe *************************************************/

	/* Pointer to nest */
	
	nst = (long)RDB[uni0 + UNIVERSE_PTR_NEST];
	CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

	/* Get pointer to regions */
	
	reg = (long)RDB[nst + NEST_PTR_REGIONS];
	CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

	/* Loop over regions */ 
	
	while (reg > VALID_PTR)
	  {
	    /* Get region index */
	    
	    i = (long)RDB[reg + NEST_REG_IDX];

	    /* Check index */

	    if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
	      Die(FUNCTION_NAME, "Invalid region index (nest)");

	    /* Calculate zone index */

	    idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);

	    /* Check maximum */

	    if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
	      Die(FUNCTION_NAME, "Indexing error (nest)");

	    /* Pointer to cell */

	    cell = (long)RDB[reg + NEST_REG_PTR_CELL];
	    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	    /* Check fill and material pointers */
	    
	    if ((uni = (long)RDB[reg + NEST_REG_PTR_FILL]) > VALID_PTR)
	      {
		/* Filled region, call recursively */
		
		MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1);
	      }
	    else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
	      {
		/* Material region, divide */

		DivideZone(mat, idx1, uni0, j);
	      }

	    /* Next region */

	    reg = NextItem(reg);
	  }
	
	/* Break case */
	    
	break;
	
	/*********************************************************************/
      }
    case UNIVERSE_TYPE_CELL:
      {
	/*********************************************************************/

	/***** Cell universe *************************************************/
	
	/* Pointer to cell list */
	
	lst = (long)RDB[uni0 + UNIVERSE_PTR_CELL_LIST];
	CheckPointer(FUNCTION_NAME, "(lst)", DATA_ARRAY, lst);  
	
	/* Loop over cell list */
	
	while (lst > VALID_PTR)
	  {
	    /* Get region index */
	    
	    i = (long)RDB[lst + CELL_LIST_REG_IDX];

	    /* Check index */

	    if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
	      Die(FUNCTION_NAME, "Invalid region index (cell)");

	    /* Calculate zone index */

	    idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);

	    /* Check maximum */

	    if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
	      Die(FUNCTION_NAME, "Indexing error (cell)");

	    /* Pointer to cell */
	    
	    cell = (long)RDB[lst + CELL_LIST_PTR_CELL];
	    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	    /* Check fill and material pointers */
	    
	    if ((uni = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
	      {
		/* Filled region, call recursively */

		MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1);
	      }
	    else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
	      {
		/* Material region, divide */
		
		DivideZone(mat, idx1, uni0, j);
	      }

	    /* Next */

	    lst = NextItem(lst);
	  }
	
	/* Break case */
	
	break;
	
	/*********************************************************************/
      }
    case UNIVERSE_TYPE_LATTICE:
      {
	/*****************************************************************/

	/***** Lattice universe **********************************************/

	/* Pointer to lattice */
	
	lat = (long)RDB[uni0 + UNIVERSE_PTR_LAT];
	CheckPointer(FUNCTION_NAME, "(lat)", DATA_ARRAY, lat);
	
	/* Check type */
	
	if ((long)RDB[lat + LAT_TYPE] == LAT_TYPE_CLU)
	  {
	    /*****************************************************************/

	    /***** Circular array ********************************************/
	    
	    /* Get pointer to rings */
	    
	    reg = (long)RDB[lat + LAT_PTR_FILL];
	    CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

	    /* Reset index */

	    i = 0;
	    
	    /* Loop over rings */ 
	    
	    while (reg > VALID_PTR)
	      {
		/* Pointer to items */ 
		
		ptr = (long)RDB[reg + RING_PTR_FILL];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	
		/* Loop over items */
		
		while ((uni = (long)RDB[ptr]) > VALID_PTR)
		  {
		    /* Check index */

		    if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
		      Die(FUNCTION_NAME, "Invalid region index (lat1)");

		    /* Calculate zone index */

		    idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);

		    /* Check maximum */

		    if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
		      Die(FUNCTION_NAME, "Indexing error (lat1)");

		    /* Update region index */
		    
		    i++; 

		    /* Call recursively */

		    MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1); 
		    
		    /* Next */
		    
		    ptr++;
		  }

		/* Next ring */

		reg = NextItem(reg);
	      }
	    
	    /*****************************************************************/
	  }
	else
	  {
	    /*****************************************************************/

	    /***** Simple types **********************************************/

	    /* Pointer to items */ 
	    
	    ptr = (long)RDB[lat + LAT_PTR_FILL];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	    /* Reset index */

	    i = 0;

	    /* Loop over items (tuolla voi olla NULLPTR välissä) */
	    
	    for (n = 0; n < (long)RDB[lat + LAT_NTOT]; n++)
	      if ((uni = (long)RDB[ptr + n]) > VALID_PTR)
		{
		  /* Check index */
		  
		  if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
		    Die(FUNCTION_NAME, "Invalid region index (lat2)");

		  /* Calculate zone index */
		  
		  idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);
		  
		  /* Check maximum */
		  
		  if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
		    Die(FUNCTION_NAME, "Indexing error (lat2)");
		  
		  /* Update region index */
		  
		  i++; 
		  
		  /* Call recursively */

		  MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1); 
		}
	    
	    /*****************************************************************/
	  }

	/* Break case */
	
	break;
	
	/*********************************************************************/
      }
    case UNIVERSE_TYPE_PBED:
      {
	/*********************************************************************/

	/***** Explicit stochastic geometry **********************************/
	
	/* Pointer to geometry */
	
	pbd = (long)RDB[uni0 + UNIVERSE_PTR_PBED];
	CheckPointer(FUNCTION_NAME, "(pbd)", DATA_ARRAY, pbd);	
	
	/* Reset index */

	i = 0;

	/* Pointer to background universe */

	uni = (long)RDB[pbd + PBED_PTR_BG_UNIV];
	CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);	

	/* Check index */

	if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
	  Die(FUNCTION_NAME, "Invalid region index (pb1)");

	/* Calculate zone index */

	idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);

	/* Check maximum */
	
	if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
	  Die(FUNCTION_NAME, "Indexing error (pb1)");
	
	/* Update region index */
	
	i++; 
		
	/* Call recursively */
		
	MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1); 

	/* Check maximum level (Tässä on tarkoituksena että loopataan   */
	/* partikkelityyppien eikä partikkelien yli, jos on varmaa että */
	/* niitä ei enää jaeta alialueisiin. Muutos tehtiin GC:n koko-  */
	/* sydänlaskuissa ilmenneen hyytymisongelman ratkaisemiseksi    */
	/* 29.3.2014. Tämä on aika viritelmä, mutta ainakaan sen ei     */
	/* pitäisi rikkoa muuta kuin korkeintaan pbed-geometrioiden     */
	/* palama-aluejako.) */

	if (((long)RDB[DATA_MAX_DIV_SEP_LVL] < 0) ||
	    ((long)RDB[DATA_MAX_DIV_SEP_LVL] > j - 3))
	  {
	    /* Loop over pebble types */

	    pbl = (long)RDB[pbd + PBED_PTR_PEBBLE_TYPES];
	    while (pbl > VALID_PTR)
	      {
		/* Pointer to universe (huom pointteri eri kuin pebble- */
		/* rakenteessa) */

		uni = (long)RDB[pbl + PEBTYPE_PTR_UNIV];
		CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);	
		
		/* Check index */
		
		if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
		  Die(FUNCTION_NAME, "Invalid region index (pb2)");
		
		/* Calculate zone index */
		
		idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);
		
		/* Check maximum */
		
		if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
		  Die(FUNCTION_NAME, "Indexing error (pb2)");
		
		/* Update region index */
		
		i++; 
		
		/* Call recursively */
		
		MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1); 
		
		/* Next pebble type */
		
		pbl = NextItem(pbl);
	      }
	  }
	else
	  {
	    /* Loop over pebbles */

	    pbl = (long)RDB[pbd + PBED_PTR_PEBBLES];
	    while (pbl > VALID_PTR)
	      {
		/* Pointer to universe */
		
		uni = (long)RDB[pbl + PEBBLE_PTR_UNIV];
		CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);	
		
		/* Check index */
		
		if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
		  Die(FUNCTION_NAME, "Invalid region index (pb2)");
		
		/* Calculate zone index */
		
		idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);
		
		/* Check maximum */
		
		if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
		  Die(FUNCTION_NAME, "Indexing error (pb2)");
		
		/* Update region index */
		
		i++; 
		
		/* Call recursively */
		
		MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1); 
		
		/* Next pebble */
		
		pbl = NextItem(pbl);
	      }
	  }

	/* Break case */
	
	break;

	/*********************************************************************/
      }
    case UNIVERSE_TYPE_STL:
      {
	/*********************************************************************/

	/***** Unstructured surface based geometry (STL) *********************/
	
	/* Pointer to geometry */
	
	stl = (long)RDB[uni0 + UNIVERSE_PTR_STL];
	CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);	
	
	/* Reset index */

	i = 0;

	/* Pointer to background universe */

	uni = (long)RDB[stl + STL_PTR_BG_UNIV];
	CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);	

	/* Check index */

	if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
	  Die(FUNCTION_NAME, "Invalid region index (stl)");

	/* Calculate zone index */

	idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);

	/* Check maximum */
	
	if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
	  Die(FUNCTION_NAME, "Indexing error (stl)");
	
	/* Call recursively */
		
	MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1); 

	/* Loop over solids */

	sld = (long)RDB[stl + STL_PTR_SOLIDS];
	while (sld > VALID_PTR)
	  {
	    /* Get region index */

	    i = (long)RDB[sld + STL_SOLID_REG_IDX];

	    /* Check index */

	    if (i > (long)RDB[lvl + LVL_MAX_REGIONS] - 1)
	      Die(FUNCTION_NAME, "Invalid region index (STL solid)");

	    /* Check numerical limit */

	    if (idx0 > LONG_MAX - i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]))
	      Die(FUNCTION_NAME, "Multiplier exceed maximum of long int");

	    /* Calculate zone index */

	    idx1[j] = idx0 + i*((long)RDB[lvl + LVL_ZONE_IDX_MULT]);

	    /* Check maximum */

	    if (idx1[j] > (long)RDB[lvl + LVL_CUM_MAX_REGIONS] - 1)
	      Die(FUNCTION_NAME, "Indexing error (cell)");

	    /* Pointer to cell */
	    
	    cell = (long)RDB[sld + STL_SOLID_PTR_CELL];
	    CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);

	    /* Check fill and material pointers */
	    
	    if ((uni = (long)RDB[cell + CELL_PTR_FILL]) > VALID_PTR)
	      {
		/* Filled region, call recursively */

		MakeDepletionZones(loc0, uni, lvl, recu, idx1[j], idx1);
	      }
	    else if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
	      {
		/* Material region, divide */
		
		DivideZone(mat, idx1, uni0, j);
	      }

	    /* Pointer to next solid */
		
	    sld = NextItem(sld);
	  }

	/* Break case */
	
	break;

	/*********************************************************************/
      }
    default:
      {
	/* Invalid type */
	
	Die(FUNCTION_NAME, "Invalid universe type %ld",
	    (long)RDB[uni0 + UNIVERSE_TYPE]);
      }
    }

  /***************************************************************************/

  /***** Sort lists, etc. ****************************************************/

  /* Check if completed */

  if ((long)RDB[lvl + LVL_NUMBER] == 0)
    {
      /* Reset total count */

      tot = 0;

      /* Loop over materials and print */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
	{
	  /* Pointer to divisor */

	  if ((div = (long)RDB[mat + MATERIAL_PTR_DIV]) > VALID_PTR)
	    {
	      /* Get sizes */
	      
	      n0 = (long)RDB[mat + MATERIAL_DIV_N_ZONES];
	      n1 = (long)RDB[div + DIV_NX];
	      n2 = (long)RDB[div + DIV_NY];
	      n3 = (long)RDB[div + DIV_NZ];
	      n4 = (long)RDB[div + DIV_NRAD];
	      n5 = (long)RDB[div + DIV_NSEG];
	      
	      /* Check division */
	      
	      if ((n0 > 0) && (n0*n1*n2*n3*n4*n5 > 1))
		{
		  /* Print */
		  
		  fprintf(out, "Material %s:\n\n", 
			  GetText(mat + MATERIAL_PTR_NAME));
		  
		  if (n0 > 1)
		    fprintf(out, " - %ld cells\n", n0);
		  else if (n0 == 1)
		    fprintf(out, " - 1 cell\n");

		  if (n1 > 1)
		    fprintf(out, " - %ld sub-regions in x-direction\n", n1);

		  if (n2 > 1)
		    fprintf(out, " - %ld sub-regions in y-direction\n", n2);

		  if (n3 > 1)
		    fprintf(out, " - %ld sub-regions in z-direction\n", n3);

		  if (n4 > 1)
		    fprintf(out, " - %ld radial sub-regions\n", n4);

		  if (n5 > 1)
		    fprintf(out, " - %ld angular sub-regions\n", n5);
		  
		  if ((long)RDB[div + DIV_LIMS_CHECK] == NO)
		    fprintf(out, " - 1 extra zone\n");

		  if ((long)RDB[div + DIV_LIMS_CHECK] == YES)
		    fprintf(out, " - %ld depletion zones in total\n\n", 
			    n0*n1*n2*n3*n4*n5);
		  else
		    fprintf(out, " - %ld depletion zones in total\n\n", 
			    n0*n1*n2*n3*n4*n5 + 1);

		  /* Add to total */

		  tot = tot + n0*n1*n2*n3*n4*n5;
		}
	    }

	  /* Next material */
	  
	  mat = NextItem(mat);
	}

      /* Loop over divisors */

      if (tot > 1000)
	fprintf(out, "Sorting lists (this may take a while)...\n");

      div = (long)RDB[DATA_PTR_DIV0];
      while (div > VALID_PTR)
	{
	  /* Pointer to material list */

	  if ((ptr = (long)RDB[div + DIV_PTR_MAT_LIST]) > VALID_PTR)
	    {
	      /* Close list */

	      CloseList(ptr);

	      /* Sort */

	      SortList(ptr, DIV_MAT_LIST_ZONE_IDX, SORT_MODE_ASCEND);
	    }

	  /* Next divisor */

	  div = NextItem(div);
	}

      if (tot > 1000)
	fprintf(out, "OK.\n\n");
    }

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
