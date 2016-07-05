/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcfunc.c                               */
/*                                                                           */
/* Created:       2015/02/02 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes functional multi-physics interfaces                */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCFunc:"

/*****************************************************************************/

void ProcessIFCFunc(long loc0, long update)
{
  long mat0, mat, ptr, found;

  /***********************************************************************/

  /***** Link material ***************************************************/

  mat = -1;

  /* Check material pointer (not set in all types) */

  if ((long)RDB[loc0 + IFC_PTR_MAT] < VALID_PTR)
    Die(FUNCTION_NAME, "material poiner not set");
  else if(!update)
    {
      /* Link material and set TMS limits */

      /* Reset found flag */

      found = NO;
	  
      /* Loop over materials and find match */
	  
      mat0 = (long)RDB[DATA_PTR_M0];
      while (mat0 > VALID_PTR)
	{
	  /* Reset pointer */
	      
	  mat = -1;
	      
	  /* Compare name */
	      
	  if (CompareStr(mat0 + MATERIAL_PTR_NAME, loc0 + IFC_PTR_MAT))
	    mat = mat0;
	      
	  /* Check if material was divided for burnup calculation */
	      
	  if ((ptr = (long)RDB[mat0 + MATERIAL_DIV_PTR_PARENT]) 
	      > VALID_PTR)
	    if (CompareStr(ptr + MATERIAL_PTR_NAME, loc0 + IFC_PTR_MAT))
	      mat = mat0;
	      
	  /* Check material */
	      
	  if (mat > VALID_PTR)
	    {
	      /* Link interface and put flag */
	      
	      WDB[mat + MATERIAL_PTR_IFC] = (double)loc0;
	      WDB[mat + MATERIAL_USE_IFC] = (double)YES;		  

	      WDB[loc0 + IFC_PTR_MAT] = (double)mat;
 
	      /* Check if temperature is given */
		  
	      if (RDB[loc0 + IFC_MAX_TEMP] > 0.0)
		{
		  /* Put maximum temperature */
		      
		  if (RDB[loc0 + IFC_MAX_TEMP] > 
		      RDB[mat + MATERIAL_TMS_TMAX])
		    {
		      if(!update)
			{
			  WDB[mat + MATERIAL_TMS_TMAX] = 
			    RDB[loc0 + IFC_MAX_TEMP];

			}
		      else
			Die(FUNCTION_NAME,
			    "Material temperature above TMS majorant for material %s",
			    GetText(mat + MATERIAL_PTR_NAME));
		    }
		  /* Put minimum temperature */
		      
		  if (RDB[loc0 + IFC_MIN_TEMP] < 
		      RDB[mat + MATERIAL_TMS_TMIN])
		    {
		      if(!update)
			{
			  WDB[mat + MATERIAL_TMS_TMIN] = 
			    RDB[loc0 + IFC_MIN_TEMP];

			}
		      else
			Die(FUNCTION_NAME,
			    "Material temperature below TMS minorant for material %s",
			    GetText(mat + MATERIAL_PTR_NAME));

		    }

		  /* Set on-the-fly Doppler-broadening mode */
		      
		  if ((RDB[loc0 + IFC_MIN_TEMP] !=
		       RDB[loc0 + IFC_MAX_TEMP] )  && (!update))
		    WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

		}


	      /* Set found flag */
		  
	      found = YES;

	    }
	      
	  /* Next material */
	      
	  mat0 = NextItem(mat0);
	}
      
      /* Check that material was found */
      
      if (found == NO)
	Error(loc0, "Material %s in distribution file not defined", 
	      GetText(loc0 + IFC_PTR_MAT));

    }
  else
    {
      /* Check TMS limits */

      mat = (long)RDB[loc0 + IFC_PTR_MAT];

      /* Check if temperature is given */
		  
      if (RDB[loc0 + IFC_MAX_TEMP] > 0.0)
	{
	  /* Check maximum temperature */
		      
	  if (RDB[loc0 + IFC_MAX_TEMP] > 
	      RDB[mat + MATERIAL_TMS_TMAX])
		Die(FUNCTION_NAME,
		    "Material temperature above TMS majorant for material %s",
		    GetText(mat + MATERIAL_PTR_NAME));

	  /* Checkminimum temperature */
		      
	  if (RDB[loc0 + IFC_MIN_TEMP] < 
	      RDB[mat + MATERIAL_TMS_TMIN])
	    Die(FUNCTION_NAME,
		"Material temperature below TMS minorant for material %s",
		GetText(mat + MATERIAL_PTR_NAME));	   

	}

    }
	  
  return;

}
