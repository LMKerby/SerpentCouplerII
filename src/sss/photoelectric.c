/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : photoelectric.c                                */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2016/02/16 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Handles photoelectric effect for photons                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Photoelectric:"

/*****************************************************************************/

void Photoelectric(long mat, long rea, long part, double E, double x, double y,
		   double z, double u, double v, double w, double  wgt,
		   double  t, long id) {

  long ptd, ptr, Nssd, i, ss, ssmin, lo, Nxstotd, Nxsmin, lomin;
  double lE, Ee, beta, r1, ue, ve, we, totpexs, mue, if0, rtot, xsi, EdAR,
         EdTTB, Ed;
  const double *Ud, *Nxsd, *Exsd, *xsd, *xsssd, *Exsssd, *Exstotd, *xstotd;

  /* Avoid compiler warnings */
  Nxstotd = 0;
  Ud = NULL;
  Nxsd = NULL;
  Exsd = NULL;
  xsd = NULL;
  Exstotd = NULL;
  xstotd = NULL;

  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Size of the subshell probability array */
  Nssd = (long)RDB[ptd + PHOTON_DIST_N_PE_SS];

  if (Nssd > 0) {

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_PE_EB];
    CheckPointer(FUNCTION_NAME, "(Ud)", DATA_ARRAY, ptr);
    Ud = &RDB[ptr];

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_PE_NXS];
    CheckPointer(FUNCTION_NAME, "(Nxsd)", DATA_ARRAY, ptr);
    Nxsd = &RDB[ptr];

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_PE_SSE];
    CheckPointer(FUNCTION_NAME, "(Exsd)", DATA_ARRAY, ptr);
    Exsd = &RDB[ptr];

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_PE_SSXS];
    CheckPointer(FUNCTION_NAME, "(xsd)", DATA_ARRAY, ptr);
    xsd = &RDB[ptr];

    Nxstotd = (long)RDB[ptd + PHOTON_DIST_N_PE_TOT];

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_PE_ETOT];
    CheckPointer(FUNCTION_NAME, "(Exstotd)", DATA_ARRAY, ptr);
    Exstotd = &RDB[ptr];

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_PE_XSTOT];
    CheckPointer(FUNCTION_NAME, "(xstotd)", DATA_ARRAY, ptr);
    xstotd = &RDB[ptr];

  }

  /* Initialize */
  Ee = EdAR = EdTTB = 0.0;

  if (Nssd == 0 || E <= Ud[Nssd-1]) {
    /* Photon energy is smaller than any ionization energy above
     * DATA_PHOTON_EMIN.
     * Photoelectron is assumed to obtain all the photon energy. */
    Ee = E;
  }
  else {

    lE = log(E);

    /* Find the minimum subhell (K=0, L1=1, L2=2, etc) */
    for (ssmin = 0; ssmin < Nssd; ssmin++)
      if (E > Ud[ssmin])
        break;

    /* Find energy interval from the energy grid corresponding to minimum
     * subshell index */
    Nxsmin = (long)Nxsd[ssmin];
    Exsssd = &RDB[(long)Exsd[ssmin]];
    lomin = SearchArray(Exsssd, lE, Nxsmin);

    if (lomin == -1)
      Die(FUNCTION_NAME, "energy not found in the subshell xs data");

    /* Interpolation factor */
    if0 = (lE - Exsssd[lomin])/ (Exsssd[lomin+1] - Exsssd[lomin]);

    /* Find log energy index, the maximum index is lomin - Nxsmin + Nxstotd.
     * The idex can be smaller than the maximum because the edge energies are stored
     * twice. */
    lo = lomin - Nxsmin + Nxstotd + 1;
    while ((lo > 0) && (Exstotd[--lo] > lE));

    /* Log-log linear interpolation of the total xs */
    totpexs = exp(xstotd[lo] + (xstotd[lo+1] - xstotd[lo])*(lE - Exstotd[lo])
      / (Exstotd[lo+1] - Exstotd[lo]));

    CheckValue(FUNCTION_NAME, "lE", ":log energy not found in the interval of the total xs data", lE, Exstotd[lo], Exstotd[lo+1]);

    rtot = RandF(id)*totpexs;

    /* Calculate the shell xs */
    /* NOTE: assuming uniform energy grid! */
    xsi = 0.0;
    ss = -1;

    for (i = ssmin; i < Nssd; i++) {

      xsssd = &RDB[(long)xsd[i]];

      /* Energy index */
      lo = (long)Nxsd[i] - Nxsmin + lomin;

      /* Log-log linear interpolation of the subshell xs */
      xsi += exp(xsssd[lo] + (xsssd[lo+1] - xsssd[lo])*if0);

      if (rtot <= xsi) {
        ss = i;
        break;
      }
    }

    if (ss == -1) {
      /* Photoelectric effect occurred with an electron for which
       * ionization energy is smaller than DATA_PHOTON_EMIN.
       * Photoelectron is assumed to obtain all the photon energy. */
      Ee = E;
    }
    else {

      /* Set the electron energy */
      Ee = E - Ud[ss];

      /* Atomic relaxation, store locally deposited energy */
      EdAR = AtomicRelaxation(mat, rea, part, ss, x, y, z, wgt, t, id);

      CheckValue(FUNCTION_NAME, "EdAR", "", EdAR, 0.0, Ud[ss]);
    }
  }

  CheckValue(FUNCTION_NAME, "Ee", "", Ee, 0.0, E);


  if (Ee < RDB[DATA_PHOTON_EMIN]) {
    /* Electron energy is deposited locally */
    Ed = Ee + EdAR;
  }
  else {

    /***** Sample the direction of the photoelectron *************************/

    /* Rejection loop */
    do {
      r1 = RandF(id);
    } while (4.0*(1.0 - r1)*r1 < RandF(id));

    beta = sqrt(Ee*(Ee + 2.0*E_RESTMASS))/(Ee + E_RESTMASS);

    /* Cosine of the direction angle */
    mue = (2.0*r1 + beta - 1.0)/(2.0*beta*r1 - beta + 1.0);

    CheckValue(FUNCTION_NAME, "mue", "", mue, -1.0, 1.0);

    /* Electron direction cosines */
    ue = u;
    ve = v;
    we = w;

    /* Sanity check for mu and direction vectors (for NAN's etc.) */

    CheckValue(FUNCTION_NAME, "mue", "", mue, -1.01, 1.01);
    CheckValue(FUNCTION_NAME, "ue", "", ue, -1.01, 1.01);
    CheckValue(FUNCTION_NAME, "ve", "", ve, -1.01, 1.01);
    CheckValue(FUNCTION_NAME, "we", "", we, -1.01, 1.01);

    /* Rotate direction cosines */
    AziRot(mue, &ue, &ve, &we, id);

    /*************************************************************************/


    /* Use TTB-approximation for the photoelectron, store locally deposited energy */
    EdTTB = TTB(mat, part, Ee, x, y, z, ue, ve, we, wgt, t, 0, id);

    CheckValue(FUNCTION_NAME, "EdTTB", "", EdTTB, 0.0, Ee);

    /* Total deposited energy */
    Ed = EdTTB + EdAR;
  }

  /* Put particle back in stack */
  ToStack(part, id);

  /* Score pulse-height detector */

  PulseDet(part, mat, Ed, x, y, z, wgt, id); 
}

/*****************************************************************************/
