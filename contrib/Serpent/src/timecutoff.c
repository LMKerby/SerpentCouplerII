/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : timecutoff.c                                   */
/*                                                                           */
/* Created:       2012/09/22 (JLe)                                           */
/* Last modified: 2016/02/17 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Terminates sampled track by time                             */
/*                                                                           */
/* Comments: - Used for time cut-off and dynamic criticality source          */
/*             simulation                                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TimeCutoff:"

/*****************************************************************************/

long TimeCutoff(long trk, long part, long *cell, double *dt, double *x, 
		double *y, double *z, double u, double v, double w, double E, 
		double t, double wgt, double spd, long mode, long id)
{
  long mat, ptr, uni;
  double dl;

  /* Check particle pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Check velocity and time step */
  
  CheckValue(FUNCTION_NAME, "spd", "", spd, ZERO, SPD_C);
  CheckValue(FUNCTION_NAME, "dt", "", *dt, ZERO, INFTY);

  /* Check upper boundary */
  
  if (t + *dt > RDB[DATA_TIME_CUT_TMAX])
    {
      /* Check that tracks are stopped at outer boundary */

      if ((long)RDB[DATA_STOP_AT_BOUNDARY] == NO)
	Die(FUNCTION_NAME, "Tracks not stopped at outer boundary");

      /* Check simulation mode */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	Die(FUNCTION_NAME, "Time cut-off in criticality source mode");

      /* Calculate distance exceeding cut-off (suhtis) */
      
      dl = (t + *dt - RDB[DATA_TIME_CUT_TMAX])*spd;

      /* Adjust time step */
      
      *dt = RDB[DATA_TIME_CUT_TMAX] - t;

      /* Move particle back */
      
      *x = *x - dl*u;
      *y = *y - dl*v;
      *z = *z - dl*w;
 
      /* Get new location */
	      
      *cell = WhereAmI(*x, *y, *z, u, v, w, id);
      CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, *cell);

      /* Apply boundary conditions */
      /* NOTE: Tämän voi ehkä tarvita jos stopatboundary ei toimi niin */
      /* kuin pitäisi. */

      /*
      bc = BoundaryConditions(cell, x, y, z, &u, &v, &w, &wgt, id);
      */

      /* Get material pointer */
      
      mat = (long)RDB[*cell + CELL_PTR_MAT];
      mat = MatPtr(mat, id);

      /* Put coordinates */
	  
      WDB[part + PARTICLE_X] = *x;
      WDB[part + PARTICLE_Y] = *y;
      WDB[part + PARTICLE_Z] = *z;
      
      WDB[part + PARTICLE_U] = u;
      WDB[part + PARTICLE_V] = v;
      WDB[part + PARTICLE_W] = w;
      
      WDB[part + PARTICLE_E] = E;
      WDB[part + PARTICLE_T] = RDB[DATA_TIME_CUT_TMAX];
      WDB[part + PARTICLE_WGT] = wgt;
      WDB[part + PARTICLE_PTR_MAT] = mat;
      
      /* Check number of time bins */

      if ((long)RDB[DATA_DYN_NB] == 1)
	{
	  /* Single bin, score mean generation number */
	  
	  ptr = (long)RDB[RES_MEAN_NGEN];
	  AddBuf1D(RDB[part + PARTICLE_GEN_IDX], wgt, ptr, id, 0);

	  /* Put particle to bank or stack */

	  /* The banks might be written to file at savedynsrc.c */
	  /* so we can't move it directly to stack, when using  */
	  /* delayed neutrons (JLe: bankataanko kaikki vai vain */
	  /* viivästyneet neutronit?) */

	  if (RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DELDYN)
	    ToBank(part, id);
	  else
	    ToStack(part, id);
	}
      else
	{
	  /* Check if last bin and score generation number */
	  
	  if (RDB[DATA_DYN_TMAX] == RDB[DATA_TIME_CUT_TMAX])
	    {
	      ptr = (long)RDB[RES_MEAN_NGEN];
	      AddBuf1D(RDB[part + PARTICLE_GEN_IDX], wgt, ptr, id, 0);
	    }

	  /* Multiple bins, bank particle */

	  ToBank(part, id);
	}

      /* Score time cut-off */

      if ((long)RDB[part + PARTICLE_TYPE] == PARTICLE_TYPE_NEUTRON)
	{
	  ptr = (long)RDB[RES_TOT_NEUTRON_CUTRATE];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(1.0, wgt, ptr, id, 0);	      
	}
      else
	{
	  ptr = (long)RDB[RES_TOT_PHOTON_CUTRATE];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(1.0, wgt, ptr, id, 0);	      
	}

      /* Return time cutoff */
	  
      return TRACK_END_TCUT;  
    }
  else
    {
      /* No time cutoff but collision, store collision time */

      /* Get collision universe */
	      
      ptr = (long)RDB[DATA_PTR_COLLISION_UNI];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      uni = GetPrivateData(ptr, id);

      /* Check pointer */
  
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      /* Put time to private data*/

      ptr = RDB[uni + UNIVERSE_PTR_PRIVA_T];
      CheckPointer(FUNCTION_NAME, "(ptr)", PRIVA_ARRAY, ptr);
      PutPrivateData(ptr, t + *dt, id);
    }

  /* Exit subroutine */

  return trk;
}

/*****************************************************************************/
