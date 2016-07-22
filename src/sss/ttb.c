/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ttb.c                                          */
/*                                                                           */
/* Created:       2014/06/15 (TKa)                                           */
/* Last modified: 2015/07/03 (TKa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Thick-target bremsstrahlung approximation for electrons and  */
/*              positrons                                                    */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/


#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TTB:"

/*****************************************************************************/

double TTB(long mat, long part, double Te, double x, double y, double z, double u,
         double v, double w, double wgt, double t, long pflag, long id) {

  long ptd, ptr, nTe0, new1, brcdfptr, brpdfptr, Ncdf, Nk, idx, i;
  double  sumEk, lTe, Yk, rcdf, Ek, cdfmax, a, Ed;
  const double *Ted, *lTed, *lYkd, *brcdf, *brpdf;

  if ((Te < RDB[DATA_PHOTON_EMIN]) || (RDB[DATA_PHOTON_USE_TTB] == NO)) {
    /* Electron/positron energy is deposited locally */
    /* TODO: (!RDB[DATA_PHOTON_USE_TTB]) vois siirtää rutiineihin, joissa elektroneja
      luodaan */
    return Te;
  }

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);
  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  ptd = (long)RDB[mat + MATERIAL_PTR_TTB];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  nTe0 = (long)RDB[ptd + TTB_NE];

  ptr = (long)RDB[ptd + TTB_E];
  CheckPointer(FUNCTION_NAME, "(Te0)", DATA_ARRAY, ptr);
  Ted = &RDB[ptr];

  ptr = (long)RDB[ptd + TTB_LE];
  CheckPointer(FUNCTION_NAME, "(lTe0)", DATA_ARRAY, ptr);
  lTed = &RDB[ptr];

  /* Select positron or electron data */
  if (pflag && ((long)RDB[DATA_PHOTON_TTBPM] == YES)) {
    /* Positron */
    ptr = (long)RDB[ptd + TTB_LYP];
    CheckPointer(FUNCTION_NAME, "(lYph0)", DATA_ARRAY, ptr);
    lYkd = &RDB[ptr];

    ptr = (long)RDB[ptd + TTB_BRPCDF];
    CheckPointer(FUNCTION_NAME, "(brpcdf)", DATA_ARRAY, ptr);
    brcdfptr = ptr;

    ptr = (long)RDB[ptd + TTB_BRPPDF];
    CheckPointer(FUNCTION_NAME, "(brppdf", DATA_ARRAY, ptr);
    brpdfptr = ptr;
  }
  else {
    /* Electron */
    ptr = (long)RDB[ptd + TTB_LYE];
    CheckPointer(FUNCTION_NAME, "(lYph0)", DATA_ARRAY, ptr);
    lYkd = &RDB[ptr];

    ptr = (long)RDB[ptd + TTB_BRECDF];
    CheckPointer(FUNCTION_NAME, "(brecdf)", DATA_ARRAY, ptr);
    brcdfptr = ptr;

    ptr = (long)RDB[ptd + TTB_BREPDF];
    CheckPointer(FUNCTION_NAME, "(brepdf)", DATA_ARRAY, ptr);
    brpdfptr = ptr;
  }

  /* Find log energy NOTE: assuming log interpolated energy array */
  lTe = log(Te);
  idx = (long)((lTe - lTed[0])/(lTed[1] - lTed[0]));

  /* Check index and energy TODO: Nämä vois laittaa checkeiksi */
  if ((idx < 0) || (idx > nTe0))
    Die(FUNCTION_NAME, "Electron/positron energy not found for index %ld", idx);

  if ((Te < Ted[idx]) || (Te > Ted[idx+1]))
    Die(FUNCTION_NAME, "Energy not found in the interval: Ted[%ld] = %.5E, Te = %.5E, Ted[%ld] = %.5E", idx, Ted[idx], Te, idx+1, Ted[idx+1]);

  /* Interpolate photon yield */
  Yk = exp(lYkd[idx] + (lYkd[idx+1] - lYkd[idx])*(lTe - lTed[idx])
      /(lTed[idx+1] - lTed[idx]));

  /* Sample number of photons */
  Nk = (long)(Yk + RandF(id));

  CheckValue(FUNCTION_NAME, "Nk", "", Nk, 0, INFTY);

  if (Nk == 0) {
    /* No bremsstrahlung photons created, deposit energy locally */
    return Te;
  }

  /* Set bremsstrahlung energy cdf and pdf, the grid is selected using
   * interpolation weights */
  if ((idx == 0) || (lTed[idx] + RandF(id)*(lTed[idx+1] - lTed[idx]) <= lTe)) {
    brcdf = &RDB[(long)RDB[brcdfptr + idx + 1]];
    brpdf = &RDB[(long)RDB[brpdfptr + idx + 1]];
    Ncdf = idx + 2; /* +2 due to the index change */

    /* Interpolate the maximum cdf assuming the cdf is integrated from linearly
     * interpolated pdf (pdf(x) = ax + b) */
    a = (brpdf[idx+1] - brpdf[idx]) / (Ted[idx+1] - Ted[idx]);
    cdfmax = brcdf[idx] + 0.5*(2.0*brpdf[idx] + a*(Te - Ted[idx]))*(Te - Ted[idx]);
  }
  else {
    /* NOTE: The lower index is selected, meaning that the maximum photon energy
     * will be below the electron energy */
    brcdf = &RDB[(long)RDB[brcdfptr + idx]];
    brpdf = &RDB[(long)RDB[brpdfptr + idx]];
    Ncdf = idx + 1;
    cdfmax = brcdf[idx];
  }


  /***** Sample photon energies **********************************************/

  sumEk = 0;

  for (i = 0; i < Nk; i++) {

    /* Sample from  */
    rcdf = RandF(id)*cdfmax;
    idx = SearchArray(brcdf, rcdf, Ncdf);

    if (idx == -1)
      Die(FUNCTION_NAME, "rcdf not found");

    /* Solve the photon energy assuming linearly interpolated pdf */
    a = (brpdf[idx+1] - brpdf[idx]) / (Ted[idx+1] - Ted[idx]);
    Ek = Ted[idx] + (sqrt(brpdf[idx]*brpdf[idx] + 2.0*a*(rcdf - brcdf[idx]))
          - brpdf[idx])/a;

    /* Check Ek limits */
    CheckValue(FUNCTION_NAME, "bremsstrahlung photon energy", "",
               Ek, Ted[idx], Ted[idx+1]);

    if (Ek < RDB[DATA_PHOTON_EMIN])
      Die(FUNCTION_NAME, "Photon energy %.5E smaller than DATA_PHOTON_EMIN", Ek);
    if (Ek > Te)
      Die(FUNCTION_NAME, "Photon energy %.5E larger than electron/positron energy", Ek);

    sumEk += Ek;

    if ((sumEk > Te) && ((long)RDB[DATA_PHOTON_TTBEC] == YES)) {
      /* Sum of photon energies is above the electron energy. Sampling is stopped */

      /* Residual energy */
      sumEk -= Ek;
      Ek = Te - sumEk;

      if (Ek < RDB[DATA_PHOTON_EMIN]) {
        /* The last photon is not created */
        Nk = i;
      }
      else {
        /* Create a new photon having the residual energy */
        sumEk = Te;

        /* Correct the number of photons */
        Nk = i+1;

        new1 = DuplicateParticle(part, id);
        WDB[new1 + PARTICLE_X] = x;
        WDB[new1 + PARTICLE_Y] = y;
        WDB[new1 + PARTICLE_Z] = z;
        WDB[new1 + PARTICLE_WGT] = wgt; /*TODO: Onko oikein? */
        WDB[new1 + PARTICLE_T] = t;
        WDB[new1 + PARTICLE_E] = Ek;
        WDB[new1 + PARTICLE_U] = u;
        WDB[new1 + PARTICLE_V] = v;
        WDB[new1 + PARTICLE_W] = w;

        ToQue(new1, id);
      }

      break;
    }

    /* Create new photon */
    new1 = DuplicateParticle(part, id);

    WDB[new1 + PARTICLE_X] = x;
    WDB[new1 + PARTICLE_Y] = y;
    WDB[new1 + PARTICLE_Z] = z;
    WDB[new1 + PARTICLE_WGT] = wgt; /*TODO: Onko oikein? */
    WDB[new1 + PARTICLE_T] = t;
    WDB[new1 + PARTICLE_E] = Ek;

    /* NOTE: Photon direction is not sampled */
    WDB[new1 + PARTICLE_U] = u;
    WDB[new1 + PARTICLE_V] = v;
    WDB[new1 + PARTICLE_W] = w;

    ToQue(new1, id);
  }

  /***************************************************************************/


  /* Locally deposited energy */
  Ed = Te - sumEk;

  /* Ed can be negative if RDB[DATA_PHOTON_TTBEC]==NO.
   * TODO: Should energy deposition be used only when RDB[DATA_PHOTON_TTBEC]==YES?
   * */
  if (Ed < 0.0)
    Ed = 0.0;

  /* TODO: Statistics: Total number of bremsstrahlung photons created */
  /* WDB[DATA_PHOTON_BREM_TOT] += Nk; */

  return Ed;
}

/*****************************************************************************/
