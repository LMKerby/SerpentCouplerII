/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : iteratetfb.c                                   */
/*                                                                           */
/* Created:       2012/01/11 (JLe)                                           */
/* Last modified: 2012/11/04 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Calculates material temperatures for temperature feedback    */
/*                                                                           */
/* Comments: - Tää ei välttämättä vielä toimi reproducible MPI -moodissa     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "IterateTFB:"

#define MAX_ITER 1000

/*****************************************************************************/

void IterateTFB()
{
  long tfb, reg, nst, mat, ptr, i, print, tb;
  double Tlim, Tout, Tin, r0, r1, Lhr, k, PowIn, c0, c1, c2, T, gtot, Tp, Ti;
  double ef, ec, hrad, hgas, Pres, dens;

  /* Check flag */

  if ((long)RDB[DATA_USE_TFB] == NO)
    return;

  /* Start timer */
  
  StartTimer(TIMER_TFB);

  /* Printing for debug */

  print = 0;

  if (print > 0) 
    fprintf(out,"\n"); 

  /* Get time bin index */

  tb = (long)RDB[DATA_DYN_TB];

  /* Loop over feedbacks */

  tfb = (long)RDB[DATA_PTR_TFB0];
  while (tfb > VALID_PTR)
    {
      /* Get boundary condition */

      Tlim = RDB[tfb + TFB_TLIM];
      CheckValue(FUNCTION_NAME, "Tlim", "", Tlim, 273.15, 1E+4);

      /* Get pointer to nest */

      nst = RDB[tfb + TFB_PTR_NST];
      CheckPointer(FUNCTION_NAME, "(nst)", DATA_ARRAY, nst);

      /* Update geometry */

      UpdateTFBgeom(tfb);

      /* Update heat conductivities */

      UpdateTFBhc(tfb);

      /* Print basic data concerning the current pin-type */
      
      if(print > 0)
        fprintf(out, "\nPin %s - Tlim = %1.2f: - P =%4.2f \n",
		GetText(nst + NEST_PTR_NAME), Tlim, RDB[tfb + TFB_ITER_POW]);

      /* The limit temperature is the outer temperature of the outermost */
      /* region */

      Tout = Tlim;

      /* Loop over regions starting from outside */

      reg = (long)RDB[tfb + TFB_PTR_REG_LIST];
      reg = LastItem(reg);
      
      while (reg > VALID_PTR)
        {
	  /*******************************************************************/

	  /***** Solve heat conduction in region *****************************/

          /* Pointer to material */

          mat = (long)RDB[reg + TFB_REG_PTR_MAT];
          CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

          /* Get inner and outer radii */
          
	  r0 = RDB[reg + TFB_REG_ITER_R0];
          r1 = RDB[reg + TFB_REG_ITER_R1];

          /* Get the volume average temperature (K) from previous step */
          
	  T = RDB[reg + TFB_REG_ITER_TEMP];

          /* Get the temperature distribution constants from previous step   */
          
	  c0 = RDB[reg + TFB_REG_ITER_C0];
          c1 = RDB[reg + TFB_REG_ITER_C1];
          c2 = RDB[reg + TFB_REG_ITER_C2];

          /* Get power at current region (linear power W/cm) */
          
	  Lhr = RDB[reg + TFB_REG_ITER_POW];
          PowIn = RDB[reg + TFB_REG_ITER_POWIN];

          /* Get heat conductivity (W/cmK)*/
          
	  k = RDB[reg + TFB_REG_HC];

          /* The solution                                               */
          /* Different relation has to be used for the innermost region */
          /* This calculates also the constants for                     */
          /* T(r) = c0*r^2 + c1*ln(r) + c2 */

          if (PrevItem(reg) > VALID_PTR)
            {
              /* non-central region */

              Tin = Tout + Lhr/(4*PI*k) + (PowIn/(2*PI*k) - 
		 (Lhr*pow(r0,2))/((2*PI*k)*(pow(r1,2)-pow(r0,2))))
                *log(r1/r0);

              c0 = -Lhr/(4*PI*k*(pow(r1,2)-pow(r0,2)));
              c1 = (Lhr/(4*PI*k) + Tout -Tin)/log(r1/r0);
              c2 = ((Tin - c0*pow(r0,2))*log(r1) - (Tout - c0*pow(r1,2))
                    *log(r0))/log(r1/r0);

              /* Radiative transfer to be taken in account for gas */

              if (RDB[mat + MATERIAL_MDENS] < 0.01)
		{
		  /* First guess */

		  Ti = Tout + 100;
		  Tp = Tout;

		  /* Update the emissivities */

		  ef = 0.78558 + 1.5263E-5*Ti;
		  ec = 0.809;

		  /* Iterate until convergence */
		  
		  for (i = 0; i <MAX_ITER; i++)
		    {
		      /* Update constants C1 & C2 */
		      
		      c1 = -(Ti-Tout)/(log(r1/r0));
		      c2 = (Ti*log(r1) - Tout*log(r0))/log(r1/r0);
		      
		      /* calculate the pressure in the gas gap in Pa */
		      
		      Pres = 10E6*
			(pow(RDB[reg + TFB_REG_R1],2) - 
			 pow(RDB[reg + TFB_REG_R0],2))/
			(pow(r1,2)-pow(r0,2))*(T/300.0);
		      
		      /* calculate the sum of jump distances (cm) */
		      
		      gtot = 0.7816*(k*100.0)*pow((Ti - Tout)/2.0, 0.5)/
			Pres*pow(2.0*4.002602, 0.5)/0.06*0.01; 
		      
		      /* calculate the heat transfer coefficients*/
		      
		      hgas = k/((r1 - r0)+(1E-4 + 1E-4) - 1.397E-4 + 1.8*gtot);
		      hrad = 1.0/(1.0/ef + 1.0/ec - 1.0)*5.6704E-12*
			(pow(Ti,2.0) + pow(Tout,2.0))*(Ti + Tout);
		      
		      /* Save previous inner temperature and calculate new */
		      
		      Tp = Ti;
		      Ti = Tout + PowIn/((hrad + hgas)*2.0*PI*r0);
		      
		      /* Break if converged */
		      
		      if (fabs(Tp - Ti) < 0.01)  
			break; 
		    }

		  /* Check convergence */

                  if (i == MAX_ITER)
                    Warn(FUNCTION_NAME,
			 "Gap heat transfer may not have converged.");

		  /* Set temperature */
		  
		  Tin = Ti;

		  /* Finally calculate the constants*/

		  c0 = 0.0;
		  c1 = -(Tin - Tout)/(log(r1/r0));
		  c2 = (Tin*log(r1) - Tout*log(r0))/log(r1/r0);
		}
	    }
          else
	    {
	      /* Central region */

	      Tin = Tout + Lhr/(4.0*PI*k);
	      c0 = -Lhr/(4.0*PI*k*pow(r1, 2.0));
	      c1 = 0;
	      c2 = Tout - c0*pow(r1, 2.0);
	    }

          /* Calculate the volume averaged temperature */
	  
          if (r0 == 0.0)
            T = (c0*r1*r1 + c1*(2.0*log(r1) - 1.0) + 2.0*c2)/2.0;
          else
            T = (r1*r1/2.0*(c0*r1*r1 + c1*(2.0*log(r1) - 1.0) + 2.0*c2) 
		 - r0*r0/2.0*(c0*r0*r0 + c1*(2.0*log(r0) - 1.0) + 2.0*c2))
	      /(r1*r1 - r0*r0);

          /* Print some values for control */

          if (print == 2)
	    {
	      fprintf(out,"%-8s r1 = %8.5f, k = %6.4f,",
		      GetText(mat + MATERIAL_PTR_NAME), r1, k); 
	      fprintf(out," Tin = %9.3f, Pow = %9.3f\n", 
		      Tin, Lhr); 
	    }
          else if (print == 1)
            fprintf(out, "%9.3f %9.3f\n", Tout, Tin);

          /* Write the results into memory */

          WDB[reg + TFB_REG_ITER_TEMP] = T;
          WDB[reg + TFB_REG_ITER_C0] = c0;
          WDB[reg + TFB_REG_ITER_C1] = c1;
          WDB[reg + TFB_REG_ITER_C2] = c2;

	  /*******************************************************************/

	  /***** Parameters to neutronics model ******************************/

	  /* Check radii */

	  if (r1 < r0)
	    Die(FUNCTION_NAME, "Error in radii");
	  else if (r0 < RDB[reg + TFB_REG_R0])
	    Die(FUNCTION_NAME, "Error in expansion of inner surface");
	  else if (r1 < RDB[reg + TFB_REG_R1])
	    Die(FUNCTION_NAME, "Error in expansion of outer surface");

	  /* Inside radius */
	  
	  if ((ptr = (long)RDB[reg + TFB_REG_PTR_RAD_IN]) > VALID_PTR)
	    WDB[ptr] = r1;
	    
	  /* Outside radius */
	  
	  if ((ptr = (long)RDB[reg + TFB_REG_PTR_RAD_OUT]) > VALID_PTR)
	    WDB[ptr] = r0;

	  /* Density factor */

	  WDB[mat + MATERIAL_TFB_DF] = 
	    (RDB[reg + TFB_REG_R1]*RDB[reg + TFB_REG_R1]-
	     RDB[reg + TFB_REG_R0]*RDB[reg + TFB_REG_R0])/(r1*r1 - r0*r0);

	  /* Set flag */

	  WDB[DATA_USE_DENSITY_FACTOR] = (double)YES;

	  /* Check factor */

	  if (RDB[mat + MATERIAL_TFB_DF] > 1.0)
	    {
	      Warn(FUNCTION_NAME, "df = %1.2f in %s", 
		   RDB[mat + MATERIAL_TFB_DF], 
		   GetText(mat + MATERIAL_PTR_NAME));

	      WDB[mat + MATERIAL_TFB_DF] = 1.0;
	    }

	  /*******************************************************************/

	  /***** Store values to statistics **********************************/
	  
          /* Get index */
          
	  i = (long)RDB[reg + TFB_REG_IDX];

          /* Volume-averated temperature */
          
	  ptr = (long)RDB[tfb + TFB_PTR_MEAN_VTEMP];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(T, 1.0, ptr, 0, -1, i, tb);

          /* Maximum temperature */
          
	  ptr = (long)RDB[tfb + TFB_PTR_MAX_TEMP];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(Tin, 1.0, ptr, 0, -1, i, tb);

          /* Minimum temperature */
          
	  ptr = (long)RDB[tfb + TFB_PTR_MIN_TEMP];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(Tout, 1.0, ptr, 0, -1, i, tb);

	  /* Update density */
        
	  dens = RDB[mat + MATERIAL_MDENS]*
	    (RDB[reg + TFB_REG_R1]*RDB[reg + TFB_REG_R1]-
	     RDB[reg + TFB_REG_R0]*RDB[reg + TFB_REG_R0])/(r1*r1 - r0*r0);

	  /* Density */
	  
          ptr = (long)RDB[tfb + TFB_PTR_MEAN_MDENS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(dens, 1.0, ptr, 0, -1, i, tb);

          /* Radius */

	  ptr = (long)RDB[tfb + TFB_PTR_MEAN_RAD];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

          if (i == 0)
            AddBuf(RDB[reg + TFB_REG_ITER_R0], 1.0, ptr, 0, -1, i, tb);

          AddBuf(RDB[reg + TFB_REG_ITER_R1], 1.0, ptr, 0, -1, i + 1, tb);

          /* Heat conductivity */

	  ptr = (long)RDB[tfb + TFB_PTR_MEAN_HC];
          CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
          AddBuf(RDB[reg + TFB_REG_HC], 1.0, ptr, 0, -1, i, tb);

	  /*******************************************************************/

	  /****** Store values for next cycle initial guess ******************/

	  /* JLE 12/11/04 - 2.1.10 */

	  if ((long)RDB[DATA_DYN_TB] == 0)
	    {
	      WDB[reg + TFB_REG_INI_R0] = RDB[reg + TFB_REG_ITER_R0];
	      WDB[reg + TFB_REG_INI_R1] = RDB[reg + TFB_REG_ITER_R1];
	      WDB[reg + TFB_REG_INI_C0] = RDB[reg + TFB_REG_ITER_C0];
	      WDB[reg + TFB_REG_INI_C1] = RDB[reg + TFB_REG_ITER_C1];
	      WDB[reg + TFB_REG_INI_C2] = RDB[reg + TFB_REG_ITER_C2];
	    }

	  /*******************************************************************/

          /* The outer temperature for the next ring is Tin of this ring */

          Tout = Tin;
	  
          /* Next region */

          reg = PrevItem(reg);

	  /*******************************************************************/
	}

      /* Next feedback */
      
      tfb = NextItem(tfb);
    }

  /* Stop timer */
  
  StopTimer(TIMER_TFB);
}

/*****************************************************************************/
