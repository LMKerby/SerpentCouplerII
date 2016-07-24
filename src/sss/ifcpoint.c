#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ifcpoint.c                                     */
/*                                                                           */
/* Created:       2013/03/17 (JLe)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Retrieves material density factor and temperature at a point */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IFCPoint:"

/*****************************************************************************/

void IFCPoint(long mat, double *f0, double *T0, long id)
{
  long loc0, loc1, msh, type, dim, ptr, lst, np, uni, ncol;
  long tbi, ang, nr, rad, i, ptr0, ptr1;
  double f, T, wgt, r, r2, d2, d, ex, w, x, y, z, t, dx, dy, dz;
  double Temp1, Temp0, phi, phi2;
  /* Check if interfaces are defined */

  if ((long)RDB[DATA_PTR_IFC0] < VALID_PTR)
    return;

  /* Check material pointer */

  if (mat < VALID_PTR)
    return;

  /* Check pointer to material-wise interface */

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_IFC]) > VALID_PTR)
    {
      /* Get pointer to root universe */
      
      uni = (long)RDB[DATA_PTR_ROOT_UNIVERSE];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
 
      /* Get interface type */
  
      type = (long)RDB[loc0 + IFC_TYPE];
    }
  else
    {
      /* Get collision universe */
      
      ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
      uni = GetPrivateData(ptr, id);
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Check pointer to universe-wise interface */

      if ((loc0 = (long)RDB[uni + UNIVERSE_PTR_IFC_FUEP]) < VALID_PTR)
	return;

      /* Put type */

      type = (long)RDB[loc0 + IFC_FUEP_TYPE];
    }

  /* Get coordinates */
		  
  ptr = RDB[uni + UNIVERSE_PTR_PRIVA_X];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  x = GetPrivateData(ptr, id);
  
  ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Y];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  y = GetPrivateData(ptr, id);

  ptr = RDB[uni + UNIVERSE_PTR_PRIVA_Z];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  z = GetPrivateData(ptr, id);

  /* Get time */

  ptr = RDB[uni + UNIVERSE_PTR_PRIVA_T];
  CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
  t = GetPrivateData(ptr, id);

  /* Reset density factor and temperature */

  f = -1.0;
  T = -1.0;

  /* Check type */
  
  if (type == IFC_TYPE_PT_AVG)
    {
      /***********************************************************************/
      
      /***** Average of point-wise values ************************************/

      /* Get pointer to search mesh */
  
      msh = (long)RDB[loc0 + IFC_PTR_SEARCH_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

      /* Get dimensions */
  
      dim = (long)RDB[loc0 + IFC_DIM];
          
      /* Reset mean density factor, temperature and weight */
	  
      f = 0.0;
      T = 0.0;
      wgt = 0.0;
      
      /* Get exclusion radius and square of radius */
      
      r = RDB[loc0 + IFC_EXCL_RAD];
      r2 = r*r;
      
      /* Get exponent */
      
      ex = RDB[loc0 + IFC_EXP];
      
      /* Get pointer to search mesh */
      
      if ((lst = MeshPtr(msh, x, y, z)) > VALID_PTR)
	lst = (long)RDB[lst];
      
      /* Loop over content */
      
      while (lst > VALID_PTR)
	{
	  /* Pointer to point */
	  
	  loc1 = (long)RDB[lst + SEARCH_MESH_CELL_CONTENT];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  
	  /* Get parameters */
	  
	  dx = x - RDB[loc1 + IFC_PT_X];
	  dy = y - RDB[loc1 + IFC_PT_Y];
	  dz = z - RDB[loc1 + IFC_PT_Z];
	  
	  /* Avoid compiler warning */
	  
	  d2 = -1.0;
	  
	  /* Calculate square distance */
	  
	  if (dim == 3)
	    d2 = dx*dx + dy*dy + dz*dz;
	  else if (dim == 2)
	    d2 = dx*dx + dy*dy;
	  else if (dim == 1)
	    d2 = dz*dz;
	  else
	    Die(FUNCTION_NAME, "Invalid dimension");

	  /* Compare to exclusion radius */
	  
	  if (d2 < r2)
	    {
	      /* Calculate distance */
	      
	      d = sqrt(d2);
	      
	      /* Calculate weight factor */

	      w = pow(d, ex);
	      CheckValue(FUNCTION_NAME, "w", "", w, 0.0, INFTY);
	      
	      /* Invert */

	      if (w < 1E-3)
		w = 1E+3;
	      else
		w = 1.0/w;

	      /* Add to values */

	      f = f + RDB[loc1 + IFC_PT_DF]*w;
	      T = T + RDB[loc1 + IFC_PT_TMP]*w;
	      wgt = wgt + w;
	    }
	  
	  /* Next */
	  
	  lst = NextItem(lst);
	}

      /* Calculate mean */
      
      if (wgt != 0.0)
	{
	  f = f/wgt;
	  T = T/wgt;
	}

      /* Avoid round-off errors */

      if (T > RDB[loc0 + IFC_MAX_TEMP])
	T = RDB[loc0 + IFC_MAX_TEMP];
      else if (T < RDB[loc0 + IFC_MIN_TEMP])
	T = RDB[loc0 + IFC_MIN_TEMP];

      /***********************************************************************/
    }
  else if (type == IFC_TYPE_REG_MESH)
    {
      /***********************************************************************/
      
      /***** Regular mesh based distribution *********************************/

      /* Get pointer to mesh */
  
      msh = (long)RDB[loc0 + IFC_PTR_SEARCH_MESH];
      CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);
            
      /* Get pointer to mesh cell */
      
      if ((ptr = MeshPtr(msh, x, y, z)) > VALID_PTR)
	{
	  /* Get pointer to data */

	  loc1 = (long)RDB[ptr];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	  
	  /* Get values */
	  
	  f = RDB[loc1 + IFC_REG_MSH_DF];
	  T = RDB[loc1 + IFC_REG_MSH_TMP];
	}

      /***********************************************************************/
    }
  else if (type == IFC_TYPE_FUNC)
    {
      /***********************************************************************/
      
      /***** User-defined functional dependence ******************************/

      /* Get number of parameters */

      np = (long)RDB[loc0 + IFC_FUNC_NP];

      /* Get pointer to parameters */

      ptr = (long)RDB[loc0 + IFC_FUNC_PTR_PARAM];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get density factor and temperature */

      UserIFC(mat, &f, &T, x, y, z, t, np, &RDB[ptr]);

      /***********************************************************************/
    }
  else if (type == IFC_TYPE_TET_MESH)
    {
      /***********************************************************************/
      
      /***** Unstructured tetrahedral mesh ***********************************/

      /* Find region */

      if ((loc1 = FindTetCell(loc0, x, y, z, id)) > VALID_PTR)
	{
	  /* Get collision number */

	  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
	  ncol = (long)GetPrivateData(ptr, id);

	  /* Put collision flag */

	  StoreValuePair(loc0 + IFC_PTR_PREV_COL_CELL, ncol, loc1, id);

	  /* Get density factor and temperature */
		  
	  f = RDB[loc1 + IFC_TET_MSH_DF];
	  T = RDB[loc1 + IFC_TET_MSH_TMP];
	}

      /***********************************************************************/
    }
  else if ((type == IFC_TYPE_FUEP) || (type == IFC_TYPE_FPIP))
    {
      /***********************************************************************/
      
      /***** Interface to fuel performance codes *****************************/

      tbi = (long)RDB[loc0 + IFC_FUEP_PTR_T];

      /* Loop over time intervals */
      while(tbi > VALID_PTR)
	{
	  /* Compare coordinates */

	  if((t >= RDB[tbi + IFC_FUEP_T_TMIN]) &&
	     (t < RDB[tbi + IFC_FUEP_T_TMAX]))
	    break;

	  tbi = NextItem(tbi);

	}

      /* Check if found */
      if (tbi < VALID_PTR)
	{
	  if(RDB[mat + MATERIAL_USE_IFC] == YES)
	    Warn(FUNCTION_NAME,"Time not found");
	  return;
	}

      loc1 = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];

      /* Loop over axial zones */

      while (loc1 > VALID_PTR)
	{    
	  /* Compare coordinates */
	      
	  if ((z >= RDB[loc1 + IFC_FUEP_AX_ZMIN]) &&
	      (z < RDB[loc1 + IFC_FUEP_AX_ZMAX]))
	    break;

	  /* Next */
	      
	  loc1 = NextItem(loc1);
	}

      /* Check pointer */

      if (loc1 < VALID_PTR)
	{
	  if(RDB[mat + MATERIAL_USE_IFC] == YES)
	    Warn(FUNCTION_NAME,"Ax not found");
	  return;
	}

      /* Find angular zone */

      ang = (long)RDB[loc1 + IFC_FUEP_AX_PTR_ANG];

      /* Get polar angle */

      phi = PolarAngle(x,y);

      while (ang > VALID_PTR)
	{

	  /* Rotate if needed */
	  if(phi > 2.0*PI+RDB[ang + IFC_FUEP_ANG_AMIN])
	    phi2 = phi - 2.0*PI;
	  else
	    phi2 = phi;

	  /* Compare coordinates */

	  if ((phi2 >= RDB[ang + IFC_FUEP_ANG_AMIN]) &&
	      (phi2 < RDB[ang + IFC_FUEP_ANG_AMAX]))
	    break;
      
	  /* Next */
      
	  ang = NextItem(ang);
	}

      /* Check if found */

      if (ang < VALID_PTR)
	{
	  if(RDB[mat + MATERIAL_USE_IFC] == YES)
	    Warn(FUNCTION_NAME,"Ang not found");
	  return;
	}

      nr = (long)RDB[ang + IFC_FUEP_ANG_N_RAD];

      /* Calculate square radius */

      r2 = x*x + y*y;

      /* Find radial zone */

      rad = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];

      CheckPointer(FUNCTION_NAME, "(rad)", DATA_ARRAY, rad);

      i = SearchArray(&RDB[rad], r2, nr);

      if(i < 0)
	{
	  if(RDB[mat + MATERIAL_USE_IFC] == YES)
	    {

	      Warn(FUNCTION_NAME, "Rad not found %E", sqrt(r2));
	      for(i=0; i < nr; i++)
		{
		  fprintf(out, "%E ", sqrt(RDB[rad + i]));
		}
	      fprintf(out, "\n nr = %ld \n", nr);
	    }
	  return;
	}

      /* Set density factor and temperature */
 
      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      f = RDB[ptr + i + 1];

      ptr0 = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP0];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      ptr1 = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP1];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      if((RDB[tbi + IFC_FUEP_T_TMIN] < (double)-INFTY/100.0 &&
	  RDB[tbi + IFC_FUEP_T_TMAX] > (double)INFTY/100.0))
	{
	  T = RDB[ptr0 + i + 1];

	  /* Interpolate with FPIP */

	  if ((type == IFC_TYPE_FPIP) && (i >= 0))
	    {

	      T = RDB[ptr0 + i] + 
		(RDB[ptr0 + i + 1] - RDB[ptr0 + i])*
		(sqrt(r2) - sqrt(RDB[rad + i]))/
		(sqrt(RDB[rad + i + 1]) - sqrt(RDB[rad + i]));
	    }

	}
      else
	{
	  /* Interpolate between timesteps */

	  T = RDB[ptr0 + i + 1] + (RDB[ptr1 + i + 1] - RDB[ptr0 + i + 1])*
	    (t - RDB[tbi + IFC_FUEP_T_TMIN])/
	    (RDB[tbi + IFC_FUEP_T_TMAX] - RDB[tbi + IFC_FUEP_T_TMIN]);

	  /* Interpolate with FPIP */

	  if((type == IFC_TYPE_FPIP) && (i >= 0))
	    {

	      Temp1 = RDB[ptr0 + i + 1] + (RDB[ptr1 + i + 1] - RDB[ptr0 + i + 1])*
		(t - RDB[tbi + IFC_FUEP_T_TMIN])/
		(RDB[tbi + IFC_FUEP_T_TMAX] - RDB[tbi + IFC_FUEP_T_TMIN]);

	      Temp0 = RDB[ptr0 + i] + (RDB[ptr1 + i] - RDB[ptr0 + i])*
		(t - RDB[tbi + IFC_FUEP_T_TMIN])/
		(RDB[tbi + IFC_FUEP_T_TMAX] - RDB[tbi + IFC_FUEP_T_TMIN]);


	      T = Temp0 + 
		(Temp1 - Temp0)*
		(sqrt(r2) - sqrt(RDB[rad + i]))/
		(sqrt(RDB[rad + i + 1]) - sqrt(RDB[rad + i]));

	    }
	}

      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Invalid interface type");

  /* Put density factor */

  *f0 = f;

  /* Check temperature */
  
  if ((T < ZERO)  && (RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE))
    {
      /* Use TMS minimum (this is used if the point */
      /* lies outide the distribution) */

      T = RDB[mat + MATERIAL_TMS_TMIN];
      CheckValue(FUNCTION_NAME, "T", "", T, ZERO, 1E+12);      

    }

  /* Put temperature */

  *T0 = T;
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
