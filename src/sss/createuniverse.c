/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : createuniverse.c                               */
/*                                                                           */
/* Created:       2010/10/06 (JLe)                                           */
/* Last modified: 2015/06/25 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Creates universes                                            */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CreateUniverse:"

/*****************************************************************************/

long CreateUniverse(long loc0, char *name, long level)
{
  long loc1, loc2, lst, ptr, uni, cell, nst, reg, lat, pbd, umsh, stl, n;

  /* Check level count */

  if (level > MAX_GEOMETRY_LEVELS)
    Error(loc0, "Maximum number of geometry levels exceeded (infinite loop?)");

  /* Compare level to maximum */

  if (level + 1 > (long)RDB[DATA_GEOM_LEVELS])
    WDB[DATA_GEOM_LEVELS] = (double)(level + 1.0);

  /***************************************************************************/
 
  /***** Check if universe exists ********************************************/

  /* Loop over universes */

  uni = RDB[DATA_PTR_U0];
  while (uni > VALID_PTR)
    {
      /* Compare names */
      
      if (!strcmp(GetText(uni + UNIVERSE_PTR_NAME), name)) 
	return uni;
      
      /* Next universe */
      
      uni = NextItem(uni);
    }

  /* Create new universe */

  uni = NewItem(DATA_PTR_U0, UNIVERSE_BLOCK_SIZE);

  /* Put name */

  WDB[uni + UNIVERSE_PTR_NAME] = (double)PutText(name);

  /* Put level */
  
  WDB[uni + UNIVERSE_LEVEL] = (double)level;

  /* Reset pointers */

  WDB[uni + UNIVERSE_PTR_CELL_LIST] = NULLPTR;
  WDB[uni + UNIVERSE_PTR_NEST] = NULLPTR;
  WDB[uni + UNIVERSE_PTR_LAT] = NULLPTR;
  WDB[uni + UNIVERSE_PTR_PBED] = NULLPTR;
  WDB[uni + UNIVERSE_PTR_UMSH] = NULLPTR;
  WDB[uni + UNIVERSE_PTR_SYM] = NULLPTR;

  /* Allocate memory for collision counter */

  AllocValuePair(uni + UNIVERSE_COL_COUNT);

  /* Allocate memory for coordinates */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[uni + UNIVERSE_PTR_PRIVA_X] = (double)ptr;

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[uni + UNIVERSE_PTR_PRIVA_Y] = (double)ptr;

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[uni + UNIVERSE_PTR_PRIVA_Z] = (double)ptr;

  /* Onko tän ajan pakko olla universe-rakenteessa? */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[uni + UNIVERSE_PTR_PRIVA_T] = (double)ptr;

  /* Allocate memory for previous region */

  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[uni + UNIVERSE_PTR_PREV_REG] = (double)ptr;

  /***************************************************************************/

  /***** Cells ***************************************************************/

  /* Loop over cells */

  cell = RDB[DATA_PTR_C0];
  while (cell > VALID_PTR)
    {

      /* Compare names and check used-flag */

      if (!((long)RDB[cell + CELL_OPTIONS] & OPT_USED))
	if (CompareStr(cell + CELL_PTR_UNI, uni + UNIVERSE_PTR_NAME))
	  {
	    /* Set used-flag */

	    SetOption(cell + CELL_OPTIONS, OPT_USED);

	    /* Put pointer */
	    
	    WDB[cell + CELL_PTR_UNI] = (double)uni;
	    
	    /* Create new item in universe cell list */
	    
	    lst = NewItem(uni + UNIVERSE_PTR_CELL_LIST, CELL_LIST_BLOCK_SIZE);
	    
	    /* Put pointer */
	    
	    WDB[lst + CELL_LIST_PTR_CELL] = (double)cell;

	    /* Allocate memory from private array */

	    ptr = AllocPrivateData(1, PRIVA_ARRAY);

	    /* Put pointer */
	    
	    WDB[lst + CELL_LIST_PTR_COUNT] = (double)ptr;

	    /* Check if fill pointer is set */
	    
	    if (RDB[cell + CELL_PTR_FILL] > VALID_PTR)
	      {
		/* Call recursively */
		
		sprintf(name, "%s", GetText(cell + CELL_PTR_FILL));
		loc1 = CreateUniverse(cell, name, level + 1);
		
		/* Put pointer */
		
		WDB[cell + CELL_PTR_FILL] = (double)loc1;
	      }
	    
	    /* Put universe type */
	    
	    WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_CELL;
	  }
      
      /* Next cell */
      
      cell = NextItem(cell);
    }

  /* Check if cells are defined */
  
  if ((ptr = (long)RDB[uni + UNIVERSE_PTR_CELL_LIST]) > 0)
    {
      /* Close list (NOTE: Tää oli processinput.c:ssä 8.11.2010 asti) */

      CloseList(ptr);

      /* Return pointer to universe */

      return uni;
    }

  /***************************************************************************/

  /***** Nests ***************************************************************/

  /* Loop over nests */

  nst = RDB[DATA_PTR_NST0];
  while (nst > VALID_PTR)
    {
      /* Compare names and check used-flag */
      
      if ((CompareStr(nst + NEST_PTR_NAME, uni + UNIVERSE_PTR_NAME)) &
	  !((long)RDB[nst + NEST_OPTIONS] & OPT_USED))
	{
	  /* Set used-flag */

	  SetOption(nst + NEST_OPTIONS, OPT_USED);

	  /* Put pointers */
	  
	  WDB[nst + NEST_PTR_UNI] = (double)uni;
	  WDB[uni + UNIVERSE_PTR_NEST] = (double)nst;

	  /* Get pointer to regions */
	      
	  reg = (long)RDB[nst + NEST_PTR_REGIONS];
	  CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

	  /* Close list */

	  CloseList(reg);

	  /* Loop over regions */ 
	      
	  while (reg > VALID_PTR)
	    {
	      /* Check if fill pointer is set */

	      if (RDB[reg + NEST_REG_PTR_FILL] > VALID_PTR)
		{
		  /* Call recursively */

		  sprintf(name, "%s", GetText(reg + NEST_REG_PTR_FILL));
		  loc1 = CreateUniverse(nst, name, level + 1);

		  /* Put pointer */

		  WDB[reg + NEST_REG_PTR_FILL] = (double)loc1;
		}

	      /* Next region */

	      reg = NextItem(reg);
	    }

	  /* Put universe type */

	  WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_NEST;

	  /* Return pointer to universe */

	  return uni;
	}

      /* Next nest */

      nst = NextItem(nst);
    }

  /***************************************************************************/

  /***** Lattices ************************************************************/

  /* Loop over lattices */

  lat = RDB[DATA_PTR_L0];
  while (lat > VALID_PTR)
    {
      /* Compare names and check used-flag */

      if ((CompareStr(lat + LAT_PTR_NAME, uni + UNIVERSE_PTR_NAME)) &
	  !((long)RDB[lat + LAT_OPTIONS] & OPT_USED))
	{
	  /* Set used-flag */

	  SetOption(lat + LAT_OPTIONS, OPT_USED);

	  /* Put pointers */

	  WDB[lat + LAT_PTR_UNI] = (double)uni;
	  WDB[uni + UNIVERSE_PTR_LAT] = (double)lat;

	  /* Check type */

	  if ((long)RDB[lat + LAT_TYPE] == LAT_TYPE_CLU)
	    {
	      /***** Circular array ******************************************/
	      
	      /* Get pointer to rings */
	      
	      reg = (long)RDB[lat + LAT_PTR_FILL];
	      CheckPointer(FUNCTION_NAME, "(reg)", DATA_ARRAY, reg);

	      /* Loop over rings */ 
	      
	      while (reg > VALID_PTR)
		{
		  /* Loop over items */ 
		  
		  ptr = (long)RDB[reg + RING_PTR_FILL];
		  while ((long)RDB[ptr] > VALID_PTR)
		    {
		      /* Call recursively */
		      
		      sprintf(name, "%s", GetText(ptr));
		      loc1 = CreateUniverse(lat, name, level + 1);
		      
		      /* Put pointer */
		      
		      WDB[ptr++] = (double)loc1;
		    }

		  /* Next region */

		  reg = NextItem(reg);
		}

	      /***************************************************************/
	    }
	  else
	    {
	      /***** Simple types ********************************************/
	      
	      /* Loop over items */ 
	  
	      ptr = (long)RDB[lat + LAT_PTR_FILL];
	      while ((long)RDB[ptr] > VALID_PTR)
		{
		  /* Get universe  name */

		  sprintf(name, "%s", GetText(ptr));

		  /* Check if intentionally undefined (dots) */

		  for (n = 0; n < strlen(name); n++)
		    if (name[n] != '.')
		      break;
		  
		  /* Call recursively or put null pointer */
		  
		  if (n < strlen(name))
		    loc1 = CreateUniverse(lat, name, level + 1);
		  else
		    loc1 = NULLPTR;
		    
		  /* Put pointer */
		  
		  WDB[ptr++] = (double)loc1;
		}

	      /***************************************************************/
	    }

	  /* Put universe type */

	  WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_LATTICE;

	  /* Return pointer to universe */

	  return uni;
	}

      /* Next lattice */

      lat = NextItem(lat);
    }

  /***************************************************************************/

  /***** Pebble-bed geometries ***********************************************/

  /* Loop over geometries */

  pbd = RDB[DATA_PTR_PB0];
  while (pbd > VALID_PTR)
    {
      /* Compare names and check used-flag */
      
      if ((CompareStr(pbd + PBED_PTR_NAME, uni + UNIVERSE_PTR_NAME)) &
	  !((long)RDB[pbd + PBED_OPTIONS] & OPT_USED))
	{
	  /* Set used-flag */

	  SetOption(pbd + PBED_OPTIONS, OPT_USED);

	  /* Put pointers */

	  WDB[pbd + PBED_PTR_UNI] = (double)uni;
	  WDB[uni + UNIVERSE_PTR_PBED] = (double)pbd;

	  /* Call recursively for background universe */
	  
	  sprintf(name, "%s", GetText(pbd + PBED_PTR_BG_UNIV));
	  loc1 = CreateUniverse(pbd, name, level + 1);

	  /* Put pointer */

	  WDB[pbd + PBED_PTR_BG_UNIV] = (double)loc1;

	  /* Loop over pebbles */
	  
	  loc1 = (long)RDB[pbd + PBED_PTR_PEBBLES];
	  while (loc1 > VALID_PTR)
	    {
	      /* Call recursively for pebble */
	  
	      sprintf(name, "%s", GetText(loc1 + PEBBLE_PTR_UNIV));
	      ptr = CreateUniverse(pbd, name, level + 1);

	      /* Put pointer */

	      WDB[loc1 + PEBBLE_PTR_UNIV] = (double)ptr;

	      /* Loop over types */

	      loc2 = (long)RDB[pbd + PBED_PTR_PEBBLE_TYPES];
	      while(loc2 > VALID_PTR)
		{
		  /* Compare universe pointer */

		  if ((long)RDB[loc2 + PEBTYPE_PTR_UNIV] == ptr)
		    {
		      /* Add counter */

		      WDB[loc2 + PEBTYPE_COUNT] = 
			RDB[loc2 + PEBTYPE_COUNT] + 1.0;

		      /* Break loop */

		      break;
		    }

		  /* Next type */

		  loc2 = NextItem(loc2);
		}

	      /* Check pointer */

	      if (loc2 < VALID_PTR)
		{
		  /* No previous definition */

		  loc2 = NewItem(pbd + PBED_PTR_PEBBLE_TYPES, 
				 PEBTYPE_BLOCK_SIZE);

		  /* Put universe pointer */

		  WDB[loc2 + PEBTYPE_PTR_UNIV] = (double)ptr;

		  /* Init counter */

		  WDB[loc2 + PEBTYPE_COUNT] = 1.0;
		}

	      /* Next pebble */

	      loc1 = NextItem(loc1);
	    }

	  /* Put universe type */

	  WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_PBED;

	  /* Return pointer to universe */

	  return uni;
	}

      /* Next geometry */

      pbd = NextItem(pbd);
    }

  /***************************************************************************/

  /***** Unstructured mesh based geometries **********************************/

  /* Loop over geometries */

  umsh = RDB[DATA_PTR_UMSH0];
  while (umsh > VALID_PTR)
    {
      /* Compare names and check used-flag */
      
      if ((CompareStr(umsh + UMSH_PTR_NAME, uni + UNIVERSE_PTR_NAME)) &
	  !((long)RDB[umsh + UMSH_OPTIONS] & OPT_USED))
	{
	  /* Set used-flag */

	  SetOption(umsh + UMSH_OPTIONS, OPT_USED);

	  /* Put pointers */

	  WDB[umsh + UMSH_PTR_UNI] = (double)uni;
	  WDB[uni + UNIVERSE_PTR_UMSH] = (double)umsh;

	  /* Call recursively for background universe */
	  
	  sprintf(name, "%s", GetText(umsh + UMSH_PTR_BG_UNIV));
	  loc1 = CreateUniverse(umsh, name, level + 1);

	  /* Put pointer */

	  WDB[umsh + UMSH_PTR_BG_UNIV] = (double)loc1;

	  /* Put universe type */

	  WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_UMSH;

	  /* Allocate memory for next cell */

	  AllocValuePair(uni + UNIVERSE_PTR_NEXT_CELL);

	  /* Return pointer to universe */

	  return uni;
	}

      /* Next geometry */

      umsh = NextItem(umsh);
    }

  /***************************************************************************/

  /***** STL geometries ******************************************************/

  /* Loop over geometries */
 
  stl = RDB[DATA_PTR_STL0];
  while (stl > VALID_PTR)
    {
      /* Compare names and check used-flag */
      
      if ((CompareStr(stl + STL_PTR_NAME, uni + UNIVERSE_PTR_NAME)) &
	  !((long)RDB[stl + STL_OPTIONS] & OPT_USED))
	{
	  /* Set used-flag */

	  SetOption(stl + STL_OPTIONS, OPT_USED);

	  /* Put pointers */

	  WDB[stl + STL_PTR_UNI] = (double)uni;
	  WDB[uni + UNIVERSE_PTR_STL] = (double)stl;

	  /* Call recursively for background universe */
	  
	  sprintf(name, "%s", GetText(stl + STL_PTR_BG_UNIV));
	  loc1 = CreateUniverse(stl, name, level + 1);

	  /* Put pointer */

	  WDB[stl + STL_PTR_BG_UNIV] = (double)loc1;

	  /* Put universe type */

	  WDB[uni + UNIVERSE_TYPE] = (double)UNIVERSE_TYPE_STL;

	  /* Return pointer to universe */

	  return uni;
	}

      /* Next geometry */

      stl = NextItem(stl);
    }

  /***************************************************************************/

  /* Universe is not defined */
  
  if ((level == 0) && ((long)RDB[DATA_PTR_ROOT_UNIVERSE] < VALID_PTR))
    Error(loc0, "Root universe %s is not defined", name);
  else
    Error(loc0, "Universe %s is not defined", name);
  
  /* Avoid compiler warning */

  return NULLPTR;
}

/*****************************************************************************/
