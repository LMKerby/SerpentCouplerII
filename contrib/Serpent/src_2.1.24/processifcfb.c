/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processifcfb.c                                 */
/*                                                                           */
/* Created:       2015/02/02 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes fuel behavior multi-physics interfaces             */
/*                                                                           */
/* Comments: - Stats are allocated in allocinterfacestat.c                   */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessIFCFB:"

/*****************************************************************************/

void ProcessIFCFB(long loc0, long update)
{
  long loc1, loc2, mat0, mat, ptr;
  long uni, nst, reg, cell, axi, surf, ang, tbi, m, nt, ntbi;
  long naxi, nang, nr, ptr1, i;
  double rad, f, r1, r0;
  double T, phi;

  /********************************************************************/

  /***** Interface for fuel performance codes *************************/

  /* First part is not executed when updating the interface */
  /* It is only for setting the TMS limits */

  if(!update)
    {
      /* Pointer to structure */

      loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

      /* Loop */

      while (loc1 > VALID_PTR)
	{
	  /* Take care of axial boundaries in 2D geometries */

	  if((long)RDB[DATA_GEOM_DIM] == 2)
	    {
	      tbi = (long)RDB[loc1 + IFC_FUEP_PTR_T];
	      while(tbi > VALID_PTR)
		{		  
		  if((long)RDB[tbi + IFC_FUEP_T_N_AX] > 1)
		    Die(FUNCTION_NAME, 
			"Multiple axial segments in 2D geometry?");
		  
		  loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
		  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);
		      
		  /* Only one axial output bin */

		  WDB[loc2 + FUEP_NZ] = (double)1;

		  /* Put axial output limits */

		  WDB[loc2 + FUEP_ZMIN] = -INFTY;
		  WDB[loc2 + FUEP_ZMAX] = INFTY;
		      
		  /* Limits for axial flux binning */
		      
		  loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];
		  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

		  /* Only one axial output bin */
		      
		  WDB[loc2 + FUEP_NZ] = (double)1;

		  /* Put axial output limits */

		  WDB[loc2 + FUEP_ZMIN] = -INFTY;
		  WDB[loc2 + FUEP_ZMAX] = INFTY;

		  /* Put limits for the axial input zone */ 

		  axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];
		  CheckPointer(FUNCTION_NAME, "axi", DATA_ARRAY, axi);

		  /* Put lower limit */

		  WDB[axi + IFC_FUEP_AX_ZMIN] = -INFTY;

		  /* Put upper limit */

		  WDB[axi + IFC_FUEP_AX_ZMAX] = INFTY;

		  /* Next */

		  tbi = NextItem(tbi);
		}
	    }
	      
	  for (m = 0; m < (long)RDB[loc1 + IFC_FUEP_N_UNI]; m++)
	    {
	      /* Get pointer to ifc nests */
		 
	      ptr = (long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST];
	      CheckPointer(FUNCTION_NAME, "ptr", DATA_ARRAY, ptr);

	      /* Find pin universe */
	      
	      uni = (long)RDB[DATA_PTR_U0];
	      while (uni > VALID_PTR)
		{
		  /* Compare */
		  
		  if (CompareStr(uni + UNIVERSE_PTR_NAME, 
				 ptr + m))
		    break;
		  
		  /* Next */
		  
		  uni = NextItem(uni);
		}
	      
	      /* Check pointer */

	      if (uni < VALID_PTR)
		Error(loc0, "Universe %s does not exist", 
		      GetText((long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST] + m));

	      /* Check universe type (tää testaa vaan että on nesti, ei */
	      /* sitä onko pinnan tyyppinä sylinteri) */

	      if ((long)RDB[uni + UNIVERSE_TYPE] != UNIVERSE_TYPE_NEST)
		Error(loc0, "Universe %s is not pin type", 
		      GetText(uni + UNIVERSE_PTR_NAME));
	      
	      /* Check that universe is not associated with another ifc */

	      if ((long)RDB[uni + UNIVERSE_PTR_IFC_FUEP] > VALID_PTR)
		Error(loc0, "Multiple interfaces for universe %s",
		      GetText(uni + UNIVERSE_PTR_NAME));

	      /* Put pointers */

	      WDB[loc1 + IFC_FUEP_PTR_UNI] = (double)uni;
	      WDB[uni + UNIVERSE_PTR_IFC_FUEP] = (double)loc1;

	      /* Get pointer to nest structure */
	      
	      nst = (long)RDB[uni + UNIVERSE_PTR_NEST];
	      CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

	      /* Pointer to time intervals */

	      if ((tbi = (long)RDB[loc1 + IFC_FUEP_PTR_T]) < VALID_PTR)
		Error(loc0, "Interface has no time intervals");

	      /* Loop over time intervals */

	      nt = 0;
	      while (tbi > VALID_PTR)
		{
		  /* Pointer to axial zones */
		      
		  if ((axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX]) < 
		      VALID_PTR)
		    Error(loc0, "Interface %ld %ld has no axial zones", 
			  tbi, nt );	      	      

		  /* Check order of radii */		 

		  /* Loop over axial zones */

		  while (axi > VALID_PTR)
		    {
		      /* Pointer to angular zones */

		      if ((ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG]) 
			  < VALID_PTR)
			Error(loc0, "Interface has no angular zones");

		      /* Loop over angular zones */
			  
		      while(ang > VALID_PTR)
			{

			  /* Pointer to radial points (hot) */
		  
			  if ((loc2 = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2]) < VALID_PTR)
			    Error(loc0, "Interface has no radial zones");

			  /* Check order (hot) */

			  rad = -1E-6;
			  for (i = 0; i < (long)RDB[ang + IFC_FUEP_ANG_N_RAD]; i++)
			    {
			      /* Compare radii */
				  
			      if (RDB[loc2 + i] >= rad)
				rad = RDB[loc2 + i];
			      else
				Error(loc0, "Radii must be in ascending order %f %f", 
				      rad, RDB[loc2 + i]);
			    }
			      
			  /* Pointer to radial points (cold) */
			      
			  if ((loc2 = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2]) < VALID_PTR)
			    Error(loc0, "Interface has no radial zones");

			  /* Check order (cold) */

			  rad = -1E-6;
			  for(i = 0; i < (long)RDB[ang + IFC_FUEP_ANG_N_RAD]; i++)
			    {
			      /* Compare radii */

			      if (RDB[loc2 + i] >= rad)
				rad = RDB[loc2 + i];
			      else
				Error(loc0, "Radii must be in ascending order %f %f", 
				      rad,RDB[loc2 + i]);
			    }

			  ang = NextItem(ang);
			}

		      axi = NextItem(axi);
		    }

		  tbi = NextItem(tbi);
		  nt++;
		}

	      /* Pointer to time intervals */

	      if ((tbi = (long)RDB[loc1 + IFC_FUEP_PTR_T]) < VALID_PTR)
		Error(loc0, "Interface has no axial zones");	      	      
	      /* Loop over time intervals */

	      while(tbi > VALID_PTR)
		{
		  axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];

		  /* Loop over axial zones to set TMS Min/Max temperatures for nest materials*/
		  /* Pilkoin näitä tarkistuksia ja prosessointia vähän useampaan eri looppiin */
		  /* Selkeyden vuoksi 9-Jul-13 (VVa) */

		  while(axi > VALID_PTR)
		    {
		      /* Pointer to angular zones */

		      ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG];
		      CheckPointer(FUNCTION_NAME, "(ang)", DATA_ARRAY, ang);		      

		      /* Loop over angular zones */

		      while(ang > VALID_PTR)
			{
			  /* Get pointer to nest regions */

			  reg = (long)RDB[nst + NEST_PTR_REGIONS];
			  while (reg > VALID_PTR)
			    {
			      r0 = 0.0;
			      r1 = 0.0;

			      /* Pointer to outer surface*/

			      surf = (long)RDB[reg + NEST_REG_PTR_SURF_IN];

			      /* Do not process the outermost (non-bounded) region */
			      /* which is usually coolant or such */

			      if(surf < VALID_PTR)
				{
				  reg = NextItem(reg);
				  continue;
				}

			      ptr=(long)RDB[surf + SURFACE_PTR_PARAMS];

			      /* Outer radius */

			      r1 = RDB[ptr + 2];

			      /* Pointer to inner surface */

			      surf = (long)RDB[reg + NEST_REG_PTR_SURF_OUT];

			      /* Get inner radius if available*/

			      if(surf > VALID_PTR)
				{
				  ptr = (long)RDB[surf + SURFACE_PTR_PARAMS];
				  r0 = RDB[ptr + 2];
				}

			      /* Pointer to cell */

			      cell = (long)RDB[reg + NEST_REG_PTR_CELL];
			      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
			  
			      /* Pointer to material */

			      if ((mat = (long)RDB[cell + CELL_PTR_MAT]) > VALID_PTR)
				{

				  loc2 = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];

				  for(i = 0; i < (long)RDB[ang + IFC_FUEP_ANG_N_RAD]; i++)
				    {
				      /* if this node is not in this nest region, */
				      /* skip NOTE: There might not be nodes in   */
				      /* every region */
			      
				      if((RDB[loc2 + i] < r0) || (RDB[loc2 + i] > r1))
					continue;

				      /* Check temperature */
				      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP0];
				      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);      

				      /* Get temperature at this node */
					  
				      T = RDB[ptr + i];
			      
				      if (T > 0.0)
					{

					  /* Adjust material temperatures */
				      
					  if (T > RDB[mat + MATERIAL_TMS_TMAX])
					    WDB[mat + MATERIAL_TMS_TMAX] = T;
				      
					  if (T < RDB[mat + MATERIAL_TMS_TMIN])
					    WDB[mat + MATERIAL_TMS_TMIN] = T;

					  /* Put limit also to parent material*/
					  if((mat0=(long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
					    {

					      /* Adjust parent material temperatures         */
					      /* Tää tarvitaan vaan, jos ei käytä div-korttia*/
					      /* vaan vanhaa burn <nr> tyyliä                */
				      
					      if (T > RDB[mat0 + MATERIAL_TMS_TMAX])
						WDB[mat0 + MATERIAL_TMS_TMAX] = T;
				      
					      if (T < RDB[mat0 + MATERIAL_TMS_TMIN])
						WDB[mat0 + MATERIAL_TMS_TMIN] = T;

					    }
					}
				      /* End of radial for loop */
				    }
			  
				  /* If there is no node in this region, the maximum */
				  /* and minimum temperatures are at the inner/outer */
				  /* surfaces */

				  /*************************************/
				  /* Check temperature at outer radius */
				  /*************************************/

				  /* Get index */

				  loc2 = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];

				  i = SearchArray(&RDB[loc2], r1, 
						  (long)RDB[ang + IFC_FUEP_ANG_N_RAD]);

				  if(i > -1)
				    {
				      /* Get pointer to temperatures  */
				      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP0];
				  
				      /* Get temperature */

				      T = RDB[ptr + i + 1];

				      /* Interpolate temperature if IFC_TYPE_FPIP */

				      if((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FPIP) 
					{
					  if(RDB[loc2 + i+1] - RDB[loc2 + i] == 0)
					    T = RDB[ptr + i];
					  else
					    T = RDB[ptr + i] +
					      (RDB[ptr + i+1] - RDB[ptr + i])*
					      (r1 - RDB[loc2 + i])/
					      (RDB[loc2 + i+1] - RDB[loc2 + i]);
					}
				  
				      if(T < 0)
					Die(FUNCTION_NAME,
					    "Negative temperature when interpolating between: %E %E", 
					    RDB[ptr + i+1], RDB[ptr + i]);

				      /* Adjust TMS temperatures */
				  
				      if (T > RDB[mat + MATERIAL_TMS_TMAX])
					WDB[mat + MATERIAL_TMS_TMAX] = T;
			      
				      if (T < RDB[mat + MATERIAL_TMS_TMIN])
					WDB[mat + MATERIAL_TMS_TMIN] = T;

				      if((mat0=(long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
					{

					  /* Adjust parent material temperatures         */
					  /* Tää tarvitaan vaan, jos ei käytä div-korttia*/
					  /* vaan vanhaa burn <nr> tyyliä                */
				      
					  if (T > RDB[mat + MATERIAL_TMS_TMAX])
					    WDB[mat0 + MATERIAL_TMS_TMAX] = T;
			      
					  if (T < RDB[mat + MATERIAL_TMS_TMIN])
					    WDB[mat0 + MATERIAL_TMS_TMIN] = T;

					}
				    
				    }
			      
				  /*************************************/
				  /* Check temperature at inner radius */
				  /* redundant if type != FPIP         */
				  /*************************************/

				  i = SearchArray(&RDB[loc2], r0, 
						  (long)RDB[ang + IFC_FUEP_ANG_N_RAD]);

				  if(i > -1)
				    {
				      /* Get pointer to temperatures  */
				      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP0];
				  
				      /* Get temperature */

				      T = RDB[ptr + i];

				      /* Interpolate temperature if IFC_TYPE_FPIP */

				      if((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FPIP) 
					{
					  if(RDB[loc2 + i+1] - RDB[loc2 + i] == 0)
					    T = RDB[ptr + i];
					  else
					    T = RDB[ptr + i] +
					      (RDB[ptr + i + 1] - RDB[ptr + i])*
					      (r0 - RDB[loc2 + i])/
					      (RDB[loc2 + i+1] - RDB[loc2 + i]);
					}
				  
				      /* Adjust TMS temperatures */
				  
				      if (T > RDB[mat + MATERIAL_TMS_TMAX])
					WDB[mat + MATERIAL_TMS_TMAX] = T;
			      
				      if (T < RDB[mat + MATERIAL_TMS_TMIN])
					WDB[mat + MATERIAL_TMS_TMIN] = T;

				      if((mat0=(long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
					{

					  /* Adjust parent material temperatures         */
					  /* Tää tarvitaan vaan, jos ei käytä div-korttia*/
					  /* vaan vanhaa burn <nr> tyyliä                */
				      
					  if (T > RDB[mat + MATERIAL_TMS_TMAX])
					    WDB[mat0 + MATERIAL_TMS_TMAX] = T;
			      
					  if (T < RDB[mat + MATERIAL_TMS_TMIN])
					    WDB[mat0 + MATERIAL_TMS_TMIN] = T;

					}
				    
				    }

				  /* Put flag */

				  WDB[mat + MATERIAL_USE_IFC] = (double)YES;

				  /* Set on-the-fly Doppler-broadening mode */
		  
				  WDB[mat + MATERIAL_TMS_MODE] = (double)YES;

				  if((mat0=(long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
				    {
				      /* Put flag */

				      WDB[mat0 + MATERIAL_USE_IFC] = (double)YES;

				      /* Set on-the-fly Doppler-broadening mode */
			  
				      WDB[mat0 + MATERIAL_TMS_MODE] = (double)YES;
			  
				    }

				}
			  
			      /* Next region */

			      reg = NextItem(reg);
			    }


			  /* Next angular zone */

			  ang = NextItem(ang);
			}
		      		      
		      /* Next axial zone */

		      axi = NextItem(axi);
		    }

		  tbi = NextItem(tbi);
		}

	    }
	  /* Next pin */

	  loc1 = NextItem(loc1);
	}

    }

  /***************************************************************************/

  /***** Additional processing for fuep type *********************************/

  /* Loop over pins */

  loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
  while (loc1 > VALID_PTR)
    {

      /* Loop over time intervals */

      tbi = (long)RDB[loc1 + IFC_FUEP_PTR_T];
      while(tbi > VALID_PTR)
	{

	  /* Loop over axial zones */

	  axi = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];

	  if((ntbi = NextItem(tbi)) > VALID_PTR)
	    naxi = (long)RDB[ntbi + IFC_FUEP_T_PTR_AX];
	  else
	    naxi = -1;

	  while(axi > VALID_PTR)
	    {
	      /* Loop over angular zones */

	      ang = (long)RDB[axi + IFC_FUEP_AX_PTR_ANG];

	      if(naxi > VALID_PTR)
		nang = (long)RDB[naxi + IFC_FUEP_AX_PTR_ANG];
	      else
		nang = -1;

	      while(ang > VALID_PTR)
		{

		  /* Calculate the parameters for the limiting planes of */
		  /* the angular segment */

		  phi = RDB[ang + IFC_FUEP_ANG_AMIN];
		  if(cos(phi)==0)
		    WDB[ang + IFC_FUEP_ANG_CMIN] = -sin(phi)*INFTY;
		  else
		    WDB[ang + IFC_FUEP_ANG_CMIN] = -sin(phi)/cos(phi);

		  phi = RDB[ang + IFC_FUEP_ANG_AMAX];
		  if(cos(phi)==0)
		    WDB[ang + IFC_FUEP_ANG_CMAX] = -sin(phi)*INFTY;
		  else
		    WDB[ang + IFC_FUEP_ANG_CMAX] = -sin(phi)/cos(phi);
		    

		  /* Loop over radial zones */

		  nr = (long)RDB[ang + IFC_FUEP_ANG_N_RAD];

		  for(i = 0; i < nr; i++)
		    {

		      /* Put EOS temperatures */

		      if(nang > VALID_PTR)
			{
			  ptr1 = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP1];
			  ptr = (long)RDB[nang + IFC_FUEP_ANG_PTR_TEMP0];
			      
			  /* Next steps TEMP 0 is this steps TEMP 1*/
			  WDB[ptr1 + i] = RDB[ptr + i];
			}
		      else
			{
			  /* Last step Temp distribution constant in time */
			  ptr1 = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP1];
			  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_TEMP0];
			      
			  /* This steps TEMP 0 is this steps TEMP 1*/
			  WDB[ptr1 + i] = RDB[ptr + i];

			}

		      /* Calculate square radius */

		      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
		      WDB[ptr + i] = RDB[ptr + i]*RDB[ptr + i];

		      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];
		      WDB[ptr + i] = RDB[ptr + i]*RDB[ptr + i];

		      /* Calculate density factor */

		      if(i > 0)
			{

			  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];
			  f = RDB[ptr + i - 1] - RDB[ptr + i];

			  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];
			  f = f/(RDB[ptr + i - 1] - RDB[ptr + i]);

			}
		      else /* First zone */
			{

			  f = 1.0;

			  ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_COLD_R2];

			  if(RDB[ptr + i] != 0.0)
			    {
			      f = RDB[ptr + i];

			      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_HOT_R2];			      
			      f = f/RDB[ptr + i];

			    }			 			  
			}

		      /* Check value */

		      if (f > 1.0)
			{
			  /*Warn(FUNCTION_NAME, "Density factor larger than 1.0, setting 1.0, outradius of zone %E", 
			    sqrt(RDB[ptr + i]));*/
			  f=1.0;
			}

		      /* Put density factor */
		      ptr = (long)RDB[ang + IFC_FUEP_ANG_PTR_DF];
		      WDB[ptr + i] = f;

		      /* Next radial zone */
	      
		    }

		  /* Next angular zone */
		  ang = NextItem(ang);
		  if(nang > VALID_PTR)
		    nang = NextItem(nang);
		}

	      /* Next axial zone */

	      axi = NextItem(axi);
	      if(naxi > VALID_PTR)
		naxi = NextItem(naxi);
	    }

	  /* Next time interval */

	  tbi = NextItem(tbi);
	      
	}
      /* Next pin */
	  
      loc1 = NextItem(loc1);
    }   
	  
  fprintf(out, "OK.\n\n");

  /***************************************************************************/
}

/*****************************************************************************/
