/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : pairproduction.c                               */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2016/02/16 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Handles pair production of photons                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PairProduction:"

/*****************************************************************************/

void PairProduction(long mat, long rea, long part, double E, double x, double y,
		    double z, double u0, double v0, double w0, double wgt,
		    double t, long id)
{
  long ptd, new1, new2;
  double E2, eps, eps2, epsmin, epsmax, tepsmin, stepsmin, xseps, xsepsmax,
      F0, G0, G, G2, F0fc, phi1, phi2, r1, epsE, Ee, Ep, betae, betap, mue,
      mup, ue, ve, we, up, vp, wp, u, v, w, sine2, a, b, Ede, Edp, Ed;
  const double *F0c, *mdxsfc;

  /***************************************************************************/

  /* Check reaction pointer */
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Pointer to coefficient data F0 */
  F0c = &RDB[(long)RDB[ptd + PHOTON_DIST_PTR_PP_F0]];
  CheckPointer(FUNCTION_NAME, "(F0c)", DATA_ARRAY,
               (long)RDB[ptd + PHOTON_DIST_PTR_PP_F0]);

  /* Pointer to the maximum differential xs */
  mdxsfc = &RDB[(long)RDB[ptd + PHOTON_DIST_PTR_PP_MDXSFC]];
  CheckPointer(FUNCTION_NAME, "(mdxsfc)", DATA_ARRAY,
               (long)RDB[ptd + PHOTON_DIST_PTR_PP_MDXSFC]);

  /***************************************************************************/


  /***** Sample the electron and positron energy *****************************/

  eps = 0.0;
  epsmin = E_RESTMASS/E;
  epsmax = 1.0 - epsmin;

  /* Sample the electron reduced energy eps */

  if (E < 2.0) {
    /* Sample eps from uniform distribution */
    eps = epsmin + RandF(id)*(epsmax - epsmin);
  }

  else if (E <= 100.0){
    /* Sample eps from the differential xs */

    tepsmin = 2.0*epsmin;
    stepsmin = sqrt(tepsmin);

    F0 = F0c[0]*stepsmin + tepsmin*(F0c[1] + F0c[2]*stepsmin + F0c[3]*tepsmin);
    F0fc = F0 + RDB[ptd + PHOTON_DIST_PP_FC];
    G0 = RDB[ptd + PHOTON_DIST_PP_G0]*epsmin;

    /* Rational function fit for the maximum differential xs */
    E2 = E*E;
    xsepsmax = (mdxsfc[0]*E2 + mdxsfc[1]*E + mdxsfc[2])
               /(E2 + mdxsfc[3]*E + mdxsfc[4]);

    /* Rejection sampling loop */
    do {

      /* Sample from a uniform distribution */
      eps = epsmin + RandF(id)*(epsmax - epsmin);
      eps2 = eps*eps;
      G = G0/(eps-eps2);

      if (G < 1.0) {
        G2 = G*G;
        phi1 = 20.867 - 3.242*G + 0.625*G2;
        phi2 = 20.209 - 1.930*G - 0.086*G2;
        xseps = (2.0*(eps2 - eps) + 1.0)*(phi1 + F0fc) +
            2.0/3.0*(eps-eps2)*(phi2 + F0fc);
      }
      else {
        phi1 = 21.12 - 4.184*log(G + 0.952);
        xseps = (4.0/3.0*(eps2-eps) + 1.0)*(phi1 + F0fc);
      }

    } while (xsepsmax*RandF(id) > xseps);

  }
  else
    Die(FUNCTION_NAME, "Photon energy exceeds the pair production limit");


  CheckValue(FUNCTION_NAME, "eps", "", eps, 0, 1.0);

  /* Calculate the electron and positron kinetic energies */
  epsE = eps*E;
  Ee = epsE - E_RESTMASS;
  Ep = E - epsE - E_RESTMASS;

  /***************************************************************************/


  /***** Calculate direction cosines for electron and positron ***************/

  /* Sample electron scattering cosine */
  do {
    betae = sqrt(Ee*(Ee + 2.0*E_RESTMASS))/(Ee + E_RESTMASS);
    r1 = RandF(id);
    mue = (2.0*r1 + betae - 1.0)/(2.0*betae*r1 - betae + 1.0);
    sine2 = 1.0 - mue*mue;

  } while (sine2 == 0.0);

  /* Sample positron scattering cosine */
  betap = sqrt(Ep*(Ep + 2.0*E_RESTMASS))/(Ep + E_RESTMASS);
  r1 = RandF(id);
  mup = (2.0*r1 + betap - 1.0)/(2.0*betap*r1 - betap + 1.0);

  /* Check */
  CheckValue(FUNCTION_NAME, "mue", "", mue, -1.0, 1.0);
  CheckValue(FUNCTION_NAME, "mup", "", mup, -1.0, 1.0);

  /* Initialize electron direction cosines */

  ue = u0;
  ve = v0;
  we = w0;

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mue", "", mue, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "ue", "", ue, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "ve", "", ve, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "we", "", we, -1.01, 1.01);

  /* Rotate */

  AziRot(mue, &ue, &ve, &we, id);

  /* Calculate positron direction cosines */
  a = sqrt((1.0 - mup*mup)/sine2);
  b = mup + a*mue;
  up = b*u0 - a*ue;
  vp = b*v0 - a*ve;
  wp = b*w0 - a*we;

  CheckValue(FUNCTION_NAME, "positron direction cosines", "",
             up*up + vp*vp + wp*wp - 1.0, -1E-4, 1E-4);

  /***************************************************************************/

  /* TTB-approximation for the electron */
  Ede = TTB(mat, part, Ee, x, y, z, ue, ve, we, wgt, t, 0, id);

  CheckValue(FUNCTION_NAME, "Ede", "", Ede, 0.0, Ee);

  /* TTB-approximation for the positron */
  Edp = TTB(mat, part, Ep, x, y, z, up, vp, wp, wgt, t, 1, id);

  CheckValue(FUNCTION_NAME, "Edp", "", Edp, 0.0, Ep);


  /***** Positron annihilation ***********************************************/

  /* Make two new photons */
  new1 = DuplicateParticle(part, id);
  new2 = DuplicateParticle(part, id);

  /* Put variables */
  WDB[new1 + PARTICLE_X] = x;
  WDB[new1 + PARTICLE_Y] = y;
  WDB[new1 + PARTICLE_Z] = z;
  WDB[new1 + PARTICLE_WGT] = wgt; /*TODO: Onko tää oikein? */
  WDB[new1 + PARTICLE_T] = t;

  WDB[new2 + PARTICLE_X] = x;
  WDB[new2 + PARTICLE_Y] = y;
  WDB[new2 + PARTICLE_Z] = z;
  WDB[new2 + PARTICLE_WGT] = wgt; /*TODO: Onko tää oikein? Tätä ei ole vanhassa versiossa... */
  WDB[new2 + PARTICLE_T] = t;

  /* Put energies */
  WDB[new1 + PARTICLE_E] = E_RESTMASS;
  WDB[new2 + PARTICLE_E] = E_RESTMASS;

  /* Sample direction isotropically (mup can be ignored) */
  IsotropicDirection(&u, &v, &w, id);

  /* Put direction cosines (opposite directions) */
  WDB[new1 + PARTICLE_U] = u;
  WDB[new1 + PARTICLE_V] = v;
  WDB[new1 + PARTICLE_W] = w;

  WDB[new2 + PARTICLE_U] = -u;
  WDB[new2 + PARTICLE_V] = -v;
  WDB[new2 + PARTICLE_W] = -w;

  /* Put photons in que */
  ToQue(new1, id);
  ToQue(new2, id);

  /* Put incident photon back in stack */
  ToStack(part, id);

  /***************************************************************************/

  /* Total deposited energy*/
  Ed = Ede + Edp;

  /* Score pulse-height detector */
  
  PulseDet(part, mat, Ed, x, y, z, wgt, id); 
}

/*****************************************************************************/
