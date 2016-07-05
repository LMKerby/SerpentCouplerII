/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : findstlsolid.c                                 */
/*                                                                           */
/* Created:       2014/03/05 (JLe)                                           */
/* Last modified: 2014/12/26 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Finds STL solid located in position                          */
/*                                                                           */
/* Comments: - The direction vector defines the direction in which the       */
/*             nearest boundary or known cell is searched.                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FindSTLSolid:"

/*****************************************************************************/

long FindSTLSolid(long stl, double x, double y, double z, 
		  double u, double v, double w, long pre, long id)
{
  long sld, msh, ptr, pts, lst, ok, i, mode;
  
  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(stl)", DATA_ARRAY, stl);

  /* Check coordinates and direction cosines */

  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);

  CheckValue(FUNCTION_NAME, "u", "", u, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "v", "", v, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "w", "", w, -1.0, 1.0);

  /***************************************************************************/

  /***** Check for preassigned information in facet mesh *********************/

  /* Pointer to facet search mesh */

  msh = (long)RDB[stl + STL_PTR_FACET_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get pointer to search mesh cell */

  if ((ptr = MeshPtr(msh, x, y, z)) < VALID_PTR)
    return NULLPTR;

  /* Check preassigned information */      
  
  if ((lst = (long)RDB[ptr]) > VALID_PTR)
    if ((sld = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT]) < VALID_PTR)
      return -sld;

  /* Enforce Delta-tracking */

  if ((long)RDB[DATA_STL_ENFORCE_DT] == YES)
    {
      ptr = (long)RDB[DATA_DT_ENFORCE_NEXT_TRACK];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      PutPrivateData(ptr, YES, id);
    }

  /***************************************************************************/
  
  /****** Short list of solids in facet mesh *********************************/

  /* Check that short list exists and loop over list */

  if (lst > VALID_PTR)
    {
      /* Check */

      if (pre == NO)
	Die(FUNCTION_NAME, "WTF?");

      while (lst > VALID_PTR)
	{
	  /* Pointer to solid */
	  
	  if ((sld = (long)RDB[lst + SEARCH_MESH_PTR_CELL_COUNT]) < VALID_PTR)
	    break;
	  
	  /* Get mode */
	  
	  mode = (long)RDB[stl + STL_SEARCH_MODE];
	  
	  /* Pointer to fail statistics */

	  pts = (long)RDB[RES_STL_RAY_TEST];
	  CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);
	  
	  /* Resampling loop */
	  
	  for (i = 0; i < 100; i++)
	    {
	      /* Score total */

	      AddBuf1D(1.0, 1.0, pts, id, 0);

	      /* Check solid */
	      
	      if ((ok = STLRayTest(sld, msh, x, y, z, u, v, w, mode, id)) 
		  == YES)
		return sld;
	      else if (ok == NO)
		break;
	      else if (ok == STL_FACET_OVERLAP)
		{
		  /* Print error */
		  
		  if (mode == STL_SEARCH_MODE_FAST)
		    Error(stl, "Overlapping facets, try ray test mode 2");
		  else
		    Die(FUNCTION_NAME, "Overlap in safe mode");
		}
	      else if (ok <= STL_RAY_TEST_FAIL_PARA)
		{
		  /* Score failure */
		  
		  AddBuf1D(1.0, 1.0, pts, id, -ok/1000);

		  /* Resample direction */
		  
		  IsotropicDirection(&u, &v, &w, id);
		}
	      else
		Die(FUNCTION_NAME, "WTF?");
	    }
		  
	  /* Check for infinite loop */
	  
	  if (i == 100)
	    {
	      /* Record error */

	      AddBuf1D(100.0, 1.0, pts, id, -STL_RAY_TEST_FAIL_STUCK/1000);

	      /* Print warning */

	      if ((long)RDB[DATA_STL_ENFORCE_DT] == NO)	      
		Note(stl, "Particle stuck on boundary, try enforcing DT");
	      else
		Note(stl, "Particle stuck on boundary");

	      /* Put point outside */

	      return NULLPTR;
	    }

	  /* Next in content list */
	  
	  lst = NextItem(lst);
	}
    }

  /* No list or search failed --> toi failure pitäis pystyä jotenkin */
  /* tunnistamaan. */
  
  /***************************************************************************/

  /***** Loop over all cells *************************************************/

  /* Pointer to solid search mesh */

  ptr = (long)RDB[stl + STL_PTR_SOLID_MESH];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, msh);

  /* Get pointer to search mesh cell */

  if ((lst = MeshPtr(ptr, x, y, z)) < VALID_PTR)
    return NULLPTR;

  /* Loop over list */

  lst = (long)RDB[lst];
  while (lst > VALID_PTR)
    {
      /* Pointer to solid */
      
      sld = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
      CheckPointer(FUNCTION_NAME, "(sld)", DATA_ARRAY, sld);

      /* Get mode */
      
      mode = (long)RDB[stl + STL_SEARCH_MODE];
      
      /* Pointer to fail statistics */
      
      pts = (long)RDB[RES_STL_RAY_TEST];
      CheckPointer(FUNCTION_NAME, "(pts)", DATA_ARRAY, pts);

      /* Resampling loop */
      
      for (i = 0; i < 100; i++)
	{
	  /* Score total */

	  AddBuf1D(1.0, 1.0, pts, id, 0);

	  /* Check solid */
	  
	  if ((ok = STLRayTest(sld, msh, x, y, z, u, v, w, mode, id)) == YES)
	    {
	      /* Check plotter mode */

	      if ((long)RDB[DATA_PLOTTER_MODE] == NO)
		{
		  /* Pointer to counter */
	      
		  ptr = (long)RDB[lst + SEARCH_MESH_PTR_CELL_COUNT];
		  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
		  
		  /* Add counter */
		  
		  AddPrivateData(ptr, 1, id);
		}
	      
	      /* Return pointer */

	      return sld;
	    }
	  else if (ok == NO)
	    break;
	  else if (ok == STL_FACET_OVERLAP)
	    {
	      /* Print error */
	      
	      if (mode == STL_SEARCH_MODE_FAST)
		Error(stl, "Overlapping facets, try ray test mode 2");
	      else
		Die(FUNCTION_NAME, "Overlap in safe mode");
	    }
	  else if (ok <= STL_RAY_TEST_FAIL_PARA)
	    {
	      /* Score failure */
	      
	      AddBuf1D(1.0, 1.0, pts, id, -ok/1000);

	      /* Resample direction */
	      
	      IsotropicDirection(&u, &v, &w, id);
	    }
	  else
	    Die(FUNCTION_NAME, "WTF?");
	}
      
      /* Check for infinite loop */
      
      if (i == 100)
	{
	  /* Record error */
	  
	  AddBuf1D(100.0, 1.0, pts, id, -STL_RAY_TEST_FAIL_STUCK/1000);
	  
	  /* Print warning */
	  
	  if ((long)RDB[DATA_STL_ENFORCE_DT] == NO)	      
	    Note(stl, "Particle stuck on boundary, try enforcing DT");
	  else
	    Note(stl, "Particle stuck on boundary");

	  /* Put point outside */
	  
	  return NULLPTR;
	}
      
      /* Next item in list */
      
      lst = NextItem(lst);
    }

  /***************************************************************************/

  /* Point is outside */

  return NULLPTR;
}

/*****************************************************************************/
