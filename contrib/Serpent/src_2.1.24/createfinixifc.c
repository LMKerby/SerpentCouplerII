/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : createfinixifc.c                                */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Creates an IFC-template for the ReadInterface()              */
/*              for the initialization of the FINIX interface                */
/*                                                                           */
/* Comments:   -Tän vois periaatteessa kirjoittaa suoraan muistiinkin        */
/*              mutta lienee järkevintä luoda kaikki rajapinnat vasta        */
/*              ReadInterface():ssa                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "CreateFinixIFC:"

/*****************************************************************************/

void CreateFinixIFC()
{
  long fib, loc0, loc1, found, nu, ptr, i;
  char tmpstr[MAX_STR];
  double mem0, mem1;

  /* Get memory size before creation of the interface */

  mem0 = RDB[DATA_ALLOC_MAIN_SIZE];

  /* Create new interface */

  loc0 = NewItem(DATA_PTR_IFC0, IFC_BLOCK_SIZE);

  /* Print filename to string */

  sprintf(tmpstr,"./Finix.ifc");

  /* Put file name and some PARAM_N_COMMON things*/

  WDB[loc0 + IFC_PTR_INPUT_FNAME] = PutText(tmpstr);
  WDB[loc0 + PARAM_PTR_NAME] = PutText("ifc");
  WDB[loc0 + PARAM_PTR_FNAME] = PutText(tmpstr);
  WDB[loc0 + PARAM_LINE] = 0;
  WDB[loc0 + IFC_PTR_MAT] = NULLPTR;

  /* Create interface */

  ReadInterface(loc0, 0);

  /* Get memory size after creation of the interface */

  mem1 = RDB[DATA_ALLOC_MAIN_SIZE];

  /* Store ifc-memory size to interface (for MPI transfer) */

  WDB[loc0 + IFC_MEM_SIZE] = mem1 - mem0;

  /* Put pointer to interface to pins */

  fib = (long)RDB[DATA_PTR_FIN0];

  while(fib > VALID_PTR)
    {

      loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

      /* Locate correct interface pin */

      while(loc1 > VALID_PTR)
	{
	  /* Reset found flag */

	  found=0;

	  /* Get number of pins in this interface pin     */
	  /* might be larger than 1 due to axial segments */

	  nu = WDB[loc1 + IFC_FUEP_N_UNI];

	  /* Get pointer to segment list */

	  ptr = (long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST];

	  /* If one of the segments corresponds to this FINIX pin */
	  /* Activate found flag */

	  for(i=0; i < nu; i++)
	    {
	      if(CompareStr(ptr + i, fib + FINIX_PTR_UNI_NAME))
		found=1;
	    }

	  if(found==1)
	    break;
	  else
	    loc1 = NextItem(loc1);

	}

      /* Check if found */

      if(loc1 < VALID_PTR)
	Die(FUNCTION_NAME, "Could not find universe");

      /* Add pointer to fuep */

      WDB[fib + FINIX_PTR_FUEP] = (double)loc1;
  
      /* Add pointer to interface filename */

      WDB[fib + FINIX_PTR_IFC_FNAME] = PutText(tmpstr);

      /* Add pointer to interface*/

      WDB[fib + FINIX_PTR_IFC] = (double)loc0;

      fib = NextItem(fib);
    }
	  
}	  

#endif

/*****************************************************************************/
