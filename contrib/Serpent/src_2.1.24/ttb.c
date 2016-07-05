/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ttb.c                                          */
/*                                                                           */
/* Created:       2014/06/15 (TKa)                                           */
/* Last modified: 2015/05/03 (TKa)                                           */
/* Version:       2.1.24                                                     */
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

void TTB(long rea, long part, double Te, double x, double y, double z, double u,
         double v, double w, double wgt, double t, long pflag, long id) {

  long ptd, ptr, nTe0, new1, brcdfptr, brpdfptr, Ncdf, Nk, idx, i;
  double  sumEk, lTe, Yk, rcdf, Ek, cdfmax, a;
  /*double Edeploc;*/
  const double *Ted, *lTed, *lYkd, *brcdf, *brpdf;

  if ((Te < RDB[DATA_PHOTON_EMIN]) || (!RDB[DATA_PHOTON_USE_TTB])) {
    /* Electron/positron energy is deposited locally */
    /* TODO: (!RDB[DATA_PHOTON_USE_TTB]) vois siirtää rutiineihin, joissa elektroneja
      luodaan */
    /*Edeploc = Te;*/
    return;
  }

  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  nTe0 = (long)RDB[ptd + PHOTON_DIST_N_TTB_NE];

  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_E];
  CheckPointer(FUNCTION_NAME, "(Te0)", DATA_ARRAY, ptr);
  Ted = &RDB[ptr];

  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_LE];
  CheckPointer(FUNCTION_NAME, "(lTe0)", DATA_ARRAY, ptr);
  lTed = &RDB[ptr];

  if (pflag) {
    /* Positron */
    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_LYP];
    CheckPointer(FUNCTION_NAME, "(lYph0)", DATA_ARRAY, ptr);
    lYkd = &RDB[ptr];

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_BRPCDF];
    CheckPointer(FUNCTION_NAME, "(brpcdf)", DATA_ARRAY, ptr);
    brcdfptr = ptr;

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_BRPPDF];
    CheckPointer(FUNCTION_NAME, "(brppdf", DATA_ARRAY, ptr);
    brpdfptr = ptr;
  }
  else {
    /* Electron */
    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_LYE];
    CheckPointer(FUNCTION_NAME, "(lYph0)", DATA_ARRAY, ptr);
    lYkd = &RDB[ptr];

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_BRECDF];
    CheckPointer(FUNCTION_NAME, "(brecdf)", DATA_ARRAY, ptr);
    brcdfptr = ptr;

    ptr = (long)RDB[ptd + PHOTON_DIST_PTR_TTB_BREPDF];
    CheckPointer(FUNCTION_NAME, "(brepdf)", DATA_ARRAY, ptr);
    brpdfptr = ptr;
  }

  /* Find log energy */
  lTe = log(Te);
  idx = SearchArray(lTed, lTe, nTe0);

  if (idx == -1)
    Die(FUNCTION_NAME, "Electron/positron energy not found");

  /* Interpolate photon yield */
  Yk = exp(lYkd[idx] + (lYkd[idx+1] - lYkd[idx])*(lTe - lTed[idx])
      /(lTed[idx+1] - lTed[idx]));

  /* Sample number of photons */
  Nk = (long)(Yk + RandF(id));

  CheckValue(FUNCTION_NAME, "bremsstrahlung photon number yield", "",
             Nk, 0, INFTY);

  if (Nk == 0) {
    /* No bremsstrahlung photons, deposit energy locally */
    /*Edeploc = Te;*/
    return;
  }

  /* Set bremsstrahlung energy cdf and pdf
   * NOTE: interpolation is omitted, the cdf and pdf of idx+1 is used */
  brcdf = &RDB[(long)RDB[brcdfptr + idx + 1]];
  brpdf = &RDB[(long)RDB[brpdfptr + idx + 1]];

  /* Interpolate the maximum cdf assuming the cdf is integrated from linearly
   * interpolated pdf (pdf(x) = ax + b) */
  a = (brpdf[idx+1] - brpdf[idx]) / (Ted[idx+1] - Ted[idx]);
  cdfmax = brcdf[idx] + 0.5*(2.0*brpdf[idx] + a*(Te - Ted[idx]))*(Te - Ted[idx]);

  Ncdf = idx + 2; /* +2 due to the index change above */
  sumEk = 0;


  /***** Sample photon energies **********************************************/

  for (i = 0; i < Nk; i++) {

    /* Sample from  */
    rcdf = RandF(id)*cdfmax;
    idx = SearchArray(brcdf, rcdf, Ncdf);

    if (idx == -1)
      Die(FUNCTION_NAME, "rcdf not found");

    /* Solve the photon energy assuming linearly interolated pdf */
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

    if (sumEk > Te) {
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
  /*Edeploc = Te - sumEk; */

  /* TODO: Statistics: Total number of bremsstrahlung photons created */
  /* WDB[DATA_PHOTON_BREM_TOT] += Nk; */


}

/*****************************************************************************/
