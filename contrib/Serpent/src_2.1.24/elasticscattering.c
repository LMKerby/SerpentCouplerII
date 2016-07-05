/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : elasticscattering.c                            */
/*                                                                           */
/* Created:       2011/02/28 (JLe)                                           */
/* Last modified: 2014/02/25 (JLe)                                           */
/* Version:       2.1.18                                                     */
/*                                                                           */
/* Description: Handles elastic scattering for neutrons                      */
/*                                                                           */
/* Comments: - From Serpent 1.1.13 (10.8.2010)                               */
/*                                                                           */
/*           - Noi nuclide temperaturet vois olla datassa jo kT:nä           */
/*                                                                           */
/*           - Vai pitäisikö kT hakea vasta targetvelocity.c:ssä?            */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ElasticScattering:"

/*****************************************************************************/

void ElasticScattering(long mat, long rea, double *E, double *u, double *v, 
		       double *w, long id)
{
  double awr, muc, Etot;
  double V, Vx, Vy, Vz, Vt, Vtx, Vty, Vtz, Px, Py, Pz;
  double Vcx, Vcy, Vcz, kT;
  long nuc;

  /***************************************************************************/

  /***** Get initial values and check ****************************************/

  /* Check reaction and material pointers */

  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Pointer to nuclide */

  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Check initial energy */

  CheckValue(FUNCTION_NAME, "E", "", *E, ZERO, INFTY);

  /* Get temperature */

  if(RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE)
    kT = GetTemp(mat, id)*KELVIN;
  else
    kT = RDB[nuc + NUCLIDE_TEMP]*KELVIN;

  /* Get awr */

  awr = RDB[nuc + NUCLIDE_AWR];

  /* Initial neutron velocity */

  V = sqrt(*E);
  
  Vx = *u*V;
  Vy = *v*V;
  Vz = *w*V;

  /* Sample initial target velocity */

  TargetVelocity(rea, *E, &Vtx, &Vty, &Vtz, *u, *v, *w, kT, id);
  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);  

  /***************************************************************************/

  /***** Remember some values before the collision ***************************/

  /* Total energy */

  if ((Etot = *E + awr*Vt*Vt) < ZERO)
    Die(FUNCTION_NAME, "Error in energy");

  /* Total momentum */

  Px = Vx + awr*Vtx;
  Py = Vy + awr*Vty;
  Pz = Vz + awr*Vtz;

  /* Check */

  if (Px*Px + Py*Py + Pz*Pz < ZERO)
    Die(FUNCTION_NAME, "Error in momentum");

  /***************************************************************************/

  /***** Transformation to C-frame *******************************************/

  /* Calculate velocity of centre-of-mass */
  
  Vcx = (Vx + awr*Vtx)/(awr + 1.0);
  Vcy = (Vy + awr*Vty)/(awr + 1.0);
  Vcz = (Vz + awr*Vtz)/(awr + 1.0);

  /* Neutron velocities in C-frame */
  
  Vx = Vx - Vcx;
  Vy = Vy - Vcy;
  Vz = Vz - Vcz;

  V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

  /* Target velocities in C-frame */

  Vtx = Vtx - Vcx;
  Vty = Vty - Vcy;
  Vtz = Vtz - Vcz;

  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);
  
  /***************************************************************************/

  /***** Scattering in C-frame ***********************************************/

  if (*E < 1E-6)
    {
      /* Assume isotropic scattering if neutron energy is low. */

      IsotropicDirection(u, v, w, id);
    }
  else
    {
      /* Sample scattering cosine in C-frame. NOTE: Varmista että energia */
      /* on L-framessa. Muutos C-frameen kertoimella awr/(1.0 + awr).     */
      
      muc = SampleMu(rea, *E, id);
      
      /* Calculate direction cosines and rotate */
      
      *u = Vx/V;
      *v = Vy/V;
      *w = Vz/V;

      AziRot(muc, u, v, w, id);
    }

  /* Velocities after collision. */
  
  Vx = *u*V;
  Vy = *v*V;
  Vz = *w*V;

  Vtx = -*u*Vt;
  Vty = -*v*Vt;
  Vtz = -*w*Vt;

  /***************************************************************************/

  /***** Transformation back to L-frame **************************************/

  /* Neutron velocities in L-frame */

  Vx = Vx + Vcx;
  Vy = Vy + Vcy;
  Vz = Vz + Vcz;

  V = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);

  /* Target velocities in L-frame */

  Vtx = Vtx + Vcx;
  Vty = Vty + Vcy;
  Vtz = Vtz + Vcz;

  Vt = sqrt(Vtx*Vtx + Vty*Vty + Vtz*Vtz);  

  /* Set neutron energy */
  
  *E = (Vx*Vx + Vy*Vy + Vz*Vz);

  /* Set direction cosines */
  
  *u = Vx/V;
  *v = Vy/V;
  *w = Vz/V;

  /***************************************************************************/

  /***** Check conservation of energy and momentum ***************************/

#ifdef DEBUG

  /* Check relative change in total energy */

  Etot = 1.0 - (*E + awr*Vt*Vt)/Etot;
  CheckValue(FUNCTION_NAME, "E", "", Etot, -1E-5, 1E-4);

  /* Check relative change in momentum components */

  if (Px > 0.0)
    {
      Px = 1.0 - (Vx + awr*Vtx)/Px;
      CheckValue(FUNCTION_NAME, "Px", "", Px, -1E-3, 1E-3);
    }

  if (Py > 0.0)
    {
      Py = 1.0 - (Vy + awr*Vty)/Py;
      CheckValue(FUNCTION_NAME, "Py", "", Py, -1E-3, 1E-3);
    }

  if (Pz > 0.0)
    {
      Pz = 1.0 - (Vz + awr*Vtz)/Pz;
      CheckValue(FUNCTION_NAME, "Pz", "", Pz, -1E-3, 1E-3);
    }

#endif
  
  /***************************************************************************/
}

/*****************************************************************************/



