#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : comptonscattering.c                            */
/*                                                                           */
/* Created:       2011/04/15 (JLe)                                           */
/* Last modified: 2016/02/16 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Handles incoherent Compton scattering of photons             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ComptonScattering:"

/* TODO:
 * - High-energy Compton? Is Klein-Nishina (with or without scattering
 *   function) adequate then?
*/

/* Local function definitions */
void PartialPivotGauss(double[3][3], double[3], double[3]);


/*****************************************************************************/

void ComptonScattering(long mat, long rea, long part, double *E0, double x, double y,
                       double z, double *u, double *v, double *w, double wgt,
                       double t, long id) {

  long ptd, ptr, Nxtd, Nssd, lo, Npzd, i, zeroflag, ss, eli,
      Nssar, ptrelncdfar, pznegative, pzminidx, ssmin;
  double mu, E, Ek0, Ee0, Ee, Sxt, xtmax, Sxtmax, xt, pzmaxd, pzmax,
      tmp, cpintmax, cpint, pzmaxabs, Uss, reln, a,
      a1, a2, a3, pz, pz2, rcpint, Ediscr, mue, mue0, ue0, ve0, we0, ue, ve,
      we, uq, vq, wq, u0, v0, w0, Ed, EdTTB, EdAR, q, pe0c, pec, cbeta, cdelta,
      calpha;
  const double *xtd, *Sxtd, *uoccupd, *Ud, *cpd, *cpintd, *pzd, *cpssd,
      *cpintssd, *extad, *cdfeln;
  static const double zeta = 57.0320106155672; /* 1e6/(sqrt(2)*h*c)*angstrom */
  static const double fscerm = FS_CONST*E_RESTMASS;
  static const double warntol = 1.e-10;      /* Warning tolerance for Doppler broadened energy */
  double Auvw[3][3];
  double buvw[3];
  double uvw[3];

  /***************************************************************************/

  /* TODO: tähän joku if jos dataa ei löydy? */

  /* Check reaction pointer */
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
  
  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Number of scattering function data */
  Nxtd = (long)RDB[ptd + PHOTON_DIST_MINC];

  /* Momentum transfers (x = sin(theta/2)/lambda) */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_VIC];  /* xt data */
  CheckPointer(FUNCTION_NAME, "(VIC)", DATA_ARRAY, ptr);
  xtd = &RDB[ptr];

  /* Incoherent scattering functions */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_INC_FF];
  CheckPointer(FUNCTION_NAME, "(INC_FF)", DATA_ARRAY, ptr);
  Sxtd = &RDB[ptr];

  /* Number of Compton profile data */
  Npzd = (long)RDB[ptd + PHOTON_DIST_N_INC_CP];

  /* Index of the minimum pz in the data */
  pzminidx = (long)RDB[ptd + PHOTON_DIST_INC_PZMINIDX];

  /* Projection of the intial electron momentum on xt */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_INC_CPPZ];
  CheckPointer(FUNCTION_NAME, "(CPPZ)", DATA_ARRAY, ptr);
  pzd = &RDB[ptr];

  /* Compton profiles */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_INC_CP];
  CheckPointer(FUNCTION_NAME, "(CP)", DATA_ARRAY, ptr);
  cpd = &RDB[ptr];

  /* Integrated Compton profiles */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_INC_CPINT];
  CheckPointer(FUNCTION_NAME, "(CPINT)", DATA_ARRAY, ptr);
  cpintd = &RDB[ptr];

  /* Number of subshells */
  Nssd = (long)RDB[ptd + PHOTON_DIST_N_UOCCUP];
  CheckValue(FUNCTION_NAME, "Nssd", "", Nssd, 1, 40);

  /* Number of electrons per subshell */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_UOCCUP];
  CheckPointer(FUNCTION_NAME, "(UOCCUP)", DATA_ARRAY, ptr);
  uoccupd = &RDB[ptr];

  /* Electron binding energies */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_UI];
  CheckPointer(FUNCTION_NAME, "(UIONIZ)", DATA_ARRAY, ptr);
  Ud = &RDB[ptr];

  /* Last Compton profile as log */
  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_INC_EXTA];
  CheckPointer(FUNCTION_NAME, "(EXTA)", DATA_ARRAY, ptr);
  extad = &RDB[ptr];

  ptr = (long)RDB[ptd + PHOTON_DIST_PTR_INC_ELNCDF];
  CheckPointer(FUNCTION_NAME, "(EXTA)", DATA_ARRAY, ptr);
  cdfeln = &RDB[ptr];

  /***************************************************************************/

  /* Store initial photon energy */
  Ek0 = *E0;

  /* Total energy of the electron before the scattering (approximately) */
  Ee0 = E_RESTMASS;

  /* The binding energy of the electron (zero when Doppler broadening is not
   * used) */
  Uss = 0.0;

  /* Check the minimum pz in the data */
  CheckValue(FUNCTION_NAME, "minimum pz", "", pzd[pzminidx], Ee0/fscerm, Ee0/fscerm);

  /* Avoid compiler warnings */
  ss = pz = cpint = 0;
  
  
  if (Ek0 > RDB[DATA_PHOTON_EKN]) {
    /* The photon energy and scattering cosine are sampled from Klein-Nishina
     * formula. */
    KleinNishina(Ek0, &E, &mu, id);
  }

  else {
    /* Use incoherent scattering function with Klein-Nishina formula */

    /*************************************************************************/

    xtmax = log(zeta*Ek0*SQRT2);

    zeroflag = (xtd[0] == 0.0 || Sxtd[0] == 0.0);

    /* Calculate the maximum incoherent scattering function */

    if (xtmax < xtd[1] && zeroflag) {
      /* Lin-lin linear interpolation, when the first element is zero */
      /* NOTE: The first elements are stored normally, the second ones as
       * log-log */
      Sxtmax = Sxtd[0] + (exp(Sxtd[1]) - Sxtd[0])*(exp(xtmax) - xtd[0])
               /(exp(xtd[1]) - xtd[0]);
    }
    else if (xtmax >= xtd[Nxtd-1]) {
      Sxtmax = exp(xtd[Nxtd-1]);
    }
    else {
      /* Find the lower boundary, the first non-log element is excluded */
      lo = SearchArray(&xtd[1], xtmax, Nxtd - 1) + 1;

      if (lo == -1)
        Die(FUNCTION_NAME, "Maximum incoherent scattering function %.5E below data minimum", exp(xtmax));

      /* Log-log linear interpolation */
      Sxtmax = exp(Sxtd[lo] + (Sxtd[lo+1] - Sxtd[lo])*(xtmax - xtd[lo])
               /(xtd[lo+1] - xtd[lo]));
    }

    /* Rejection loop */
    do {
      /* Sample scattering cosine from Klein-Nishina formula */
      KleinNishina(Ek0, &E, &mu, id);

      /* Check scattering cosine */
      CheckValue(FUNCTION_NAME, "mu", "", mu, -1.0, 1.0);

      /***** Use the incoherent form factor as a rejection function **/
      xt = log(zeta*Ek0*sqrt(1.0 - mu));

      /*TODO: check the boundaries of xt? */

      if ((xt < xtd[1]) && zeroflag) {
        Sxt = Sxtd[0] + (exp(Sxtd[1]) - Sxtd[0])*(exp(xt) - xtd[0])
                 /(exp(xtd[1]) - xtd[0]);
      }
      else {

        /* Find the lower boundary, exclude the first non-log element */
        lo = SearchArray(&xtd[1], xt, Nxtd - 1) + 1;

        if (lo == -1)
          Die(FUNCTION_NAME, "Incoherent scattering function %.5E below data minimum", exp(xtmax));

        /* Interpolate the scattering function */
        Sxt = exp(Sxtd[lo] + (Sxtd[lo+1] - Sxtd[lo])*(xt - xtd[lo])/(xtd[lo+1] - xtd[lo]));

      }

    } while (RandF(id)*Sxtmax > Sxt);


    /*************************************************************************/


    /*************************************************************************/

    if ((long)RDB[DATA_PHOTON_USE_DOPPLER]) {

      /* Doppler broadening of the photon energy, i.e., the initial momentum
       * and the binding energy of the electron are taken into account */

      /* Maximum pz in the data TODO: Tää on turha */
      pzmaxd = RDB[ptd + PHOTON_DIST_INC_PZMAX];

      /* Find the minimum shell number */
      for (ssmin = 0; ssmin < Nssd; ssmin++)
        if (Ud[ssmin] < Ek0)
          break;

      /***** Sample the subshell *********************************************/

      while (1) {

        /* Sample the shell according to the number of electrons per shell */
        reln = RandF(id)*(cdfeln[Nssd] - cdfeln[ssmin]) + cdfeln[ssmin];
        ss = SearchArray(&cdfeln[ssmin], reln, Nssd + 1 - ssmin) + ssmin;

        CheckValue(FUNCTION_NAME, "ss", "", ss, ssmin, Npzd - 1);

        /* Set subshell data */
        cpssd = &RDB[(long)cpd[ss]];
        cpintssd = &RDB[(long)cpintd[ss]];

        /* Calculate the maximum pz of the subshell (E = E0 - Ud[i]) */
        tmp = Ek0*(Ek0 - Ud[ss])*(1.0 - mu);
        pzmax = (tmp - Ud[ss]*Ee0)/(sqrt(2.0*tmp + Ud[ss]*Ud[ss]) * fscerm);

        /* Calculate the integral between [0, pzmax] */
        if (pzmax >= pzmaxd) {
          /* Integral of the exponential extrapolation */
          cpintmax = cpintssd[Npzd-1] +
              cpssd[Npzd-1]/extad[ss]*(exp(extad[ss]*(pzmax - pzmaxd)) - 1.0);
        }
        else {
          /* Linear interpolation */
          pzmaxabs = fabs(pzmax);

          /* Get the lover boundary of pzmax */
          if (pzmaxabs == pzmaxd)
            lo = Npzd - 2;
          else
            lo = SearchArray(pzd, pzmaxabs, Npzd);

          CheckValue(FUNCTION_NAME, "lo", "", lo, 0, Npzd - 2);

          a = (cpssd[lo+1] - cpssd[lo]) / (pzd[lo+1] - pzd[lo]);

          /* Calculate the integral */
          cpintmax = cpintssd[lo] + (0.5*a*(pzmaxabs - pzd[lo]) + cpssd[lo])
              *(pzmaxabs - pzd[lo]);
        }

        /* Calculate the integral between pzmax and pzmin */
        if (pzmax < 0.0)
          cpint = cpintssd[pzminidx] - cpintmax;
        else
          cpint = cpintssd[pzminidx] + cpintmax;

        /* TODO: nuo lasketut integraalit vois säilöä, nopeuttais ehkä vähäsen */

        CheckValue(FUNCTION_NAME, "cpint", "", cpint, 0.0, 1.0);

        /* Sampling not needed for single shell atoms */
        if (Nssd == 1)
          break;

        /* Accept or reject the shell */
        if (RandF(id) <= cpint)
          break;

      }

      /* Binding energy of the sampled shell */
      Uss = Ud[ss];

     /***********************************************************************/


      /***** Sample pz and solve the corresponding photon energy *************/
      while (1) {
	
         /***** Sample pz from the Compton profile of the sampled subshell ***/
	 
        pznegative = 0;

        if (pzmax > 0.0) {
          rcpint = RandF(id)*cpint;
          if (rcpint < cpintssd[pzminidx]) {
            /* pz between [pzmin, 0] */
            pznegative = 1;
          }
          else {
            /* pz between [0, pzmax] */
            rcpint = rcpint - cpintssd[pzminidx];
            pznegative = 0;
          }
        }
        else {
          /* pz between [pzmin, pzmax], both are negative */
          rcpint = cpintmax + RandF(id)*cpint;
          pznegative = 1;
        }

        if (rcpint > cpintssd[Npzd-1]) {
          /* Tail */
          lo = Npzd - 1;

          if (pznegative) /* TODO: Tää checkiksi */
            Die(FUNCTION_NAME, "negative pz");

          pz = pzmaxd + log(extad[ss]*(rcpint - cpintssd[lo])/cpssd[lo] + 1.0)/extad[ss];
        }
        else {
          /* Get the lower boundary */
          lo = SearchArray(cpintssd, rcpint, Npzd);

          CheckValue(FUNCTION_NAME, "lo", ": rcpint not found", lo, 0, Npzd - 2);

          /* First order polynomial coefficient (p(x) = a*x + b) */
          a = (cpssd[lo+1] - cpssd[lo]) / (pzd[lo+1] - pzd[lo]);

          /* Calculate pz */
          if (a == 0.0) {
            /* The slope is zero (note: cpss[lo] is always greater than zero)*/
            pz = pzd[lo] + (rcpint - cpintssd[lo]) / cpssd[lo];
          }
          else {
            /* Here we use the integral of the linearly interpolated pdf */
            pz = pzd[lo] + (sqrt(cpssd[lo]*cpssd[lo] + 2.0*a*(rcpint - cpintssd[lo]))
                  - cpssd[lo])/a;
          }
        }

        pz *= fscerm;

        /*********************************************************************/


        /***** Solve energy **************************************************/

        pz2 = pz*pz;
        a1 = Ek0*(1.0 - mu) + Ee0;
        a2 = pz2 - a1*a1;
        a3 = mu*pz2 - a1*Ee0;
        Ediscr = a3*a3 - a2*(pz2 - Ee0*Ee0);

        if (Ediscr < 0.0 || a2 == 0.0) {
          /* Sample new pz to avoid negative discriminant or division by zero.
           * Ediscr can be negative due to floating point accuracy. These are rare
           * events. */
          continue;
        }

        /* Calculate energy */
        if (pznegative)
          E = Ek0*(a3 + sqrt(Ediscr))/a2;
        else
          E = Ek0*(a3 - sqrt(Ediscr))/a2;

        /* Check energy */
        if (E > Ek0 - Uss) {
          if ((E - Ek0 + Uss)/(Ek0 - Uss) < warntol) {
            /* Photon energy can be above the maximum due to floating point accuracy. */
            continue;
          }
          else {
            Warn(FUNCTION_NAME, "Sampled photon energy %.15E above the maximum %.15E.\n"
                                "Sampling new energy... ", E, Ek0 - Uss);
            continue;
          }
        }
        else if (E < 0.0) {
          Warn(FUNCTION_NAME, "Sampled photon energy %.15E below zero.\n"
                              "Sampling new energy... ", E);
          continue;
        }
        
#ifdef DEBUG
        if (isnan(E))
          Die(FUNCTION_NAME, "NaN photon energy in Doppler broadening.");
#endif

        /* Rejection sampling (approximation to Ribberfors DDCS) */
        if (RandF(id)*(Ek0 - Uss) <= E)
          break;

        /*********************************************************************/

      }

      if (pznegative)
        pz = -pz;

      /***********************************************************************/

    }
    /*************************************************************************/

  }
  /***************************************************************************/



  /***** Set photon energy and direction *************************************/

  /* Check energy */
  CheckValue(FUNCTION_NAME, "E", "", E, 0.0, Ek0 - Uss);

  /* Set the photon energy */
  *E0 = E;

  /* Store direction cosines of incident photon */
  u0 = *u;
  v0 = *v;
  w0 = *w;

  /* Sanity check for mu and direction vectors (for NAN's etc.) */

  CheckValue(FUNCTION_NAME, "mu", "", mu, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "u", "", u0, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "v", "", v0, -1.01, 1.01);
  CheckValue(FUNCTION_NAME, "w", "", w0, -1.01, 1.01);

  /* Rotate direction cosines of the photon */
  AziRot(mu, u, v, w, id);

  /***************************************************************************/


  /***** Calculate electron energy and direction *****************************/

  /* Set the electron kinetic energy */
  Ee = Ek0 + Ee0 - E - Uss - E_RESTMASS;

  if (Ee < RDB[DATA_PHOTON_EMIN]) {
    /* Direction not needed, TTB deposited energy is equal to the electron
     * energy */
    EdTTB = Ee;
  }
  else {

    /* Magnitude and direction cosines of momentum transfer vector */
    q = sqrt(Ek0*Ek0 + E*E - 2.0*Ek0*E*mu);

    if (q < ZERO)
      Die(FUNCTION_NAME, "Zero momentum transfer vector magnitude");

    uq = (Ek0*u0 - E*(*u))/q;
    vq = (Ek0*v0 - E*(*v))/q;
    wq = (Ek0*w0 - E*(*w))/q;


    /* Electron direction */

    if ((long)RDB[DATA_PHOTON_USE_DOPPLER] && (long)RDB[DATA_PHOTON_COMP_EANG]) {

      /* Use the "pe sampling method" described in Toni's thesis */

      /* TODO:
       * - faster sampling of the pre-collision electron direction
       * - rotation matrix would probably be faster than PartialPivotGauss */

      CheckValue(FUNCTION_NAME, "momentum transfer vector direction cosines", "",
                 uq*uq + vq*vq + wq*wq - 1.0, -1E-4, 1E-4);

      /* Sample the diection of the pre-collision electron */
      do {
        IsotropicDirection(&ue0, &ve0, &we0, id);
        cdelta = ue0*uq + ve0*vq + we0*wq;
        pe0c = -pz/cdelta;
      } while (pe0c < 0);


      mue0 = ue0*u0 + ve0*v0 + we0*w0;
      pec = sqrt(q*q + pe0c*pe0c - 2*pz*q);
      mue = (Ek0 - E*mu + pe0c*mue0)/pec;
      calpha = (q*cdelta + pe0c)/pec;
      cbeta = (Ek0*mue + pe0c*calpha - pec)/E;

      if (mue < -1 || mue > 1)
        Die(FUNCTION_NAME, "mue out of bounds");

      if (calpha < -1 || calpha > 1)
        Die(FUNCTION_NAME, "calpha out of bounds");

      if (cbeta < -1 || cbeta > 1)
        Die(FUNCTION_NAME, "cbeta out of bounds");

      /* Set the matrix */
      Auvw[0][0] = ue0;
      Auvw[0][1] = ve0;
      Auvw[0][2] = we0;
      Auvw[1][0] = u0;
      Auvw[1][1] = v0;
      Auvw[1][2] = w0;
      Auvw[2][0] = *u;
      Auvw[2][1] = *v;
      Auvw[2][2] = *w;

      buvw[0] = calpha;
      buvw[1] = mue;
      buvw[2] = cbeta;

      /* Solve electron direction cosines */
      PartialPivotGauss(Auvw, buvw, uvw);
      ue = uvw[0];
      ve = uvw[1];
      we = uvw[2];

      /* TODO: tarkista kulmat ym*/
    }

    else {
      /* Approximation: electron travels in the direction of the momentum
       * transfer vector. This is equal to the free-electron scattering angle when
       * Doppler broadening is not used. */
      ue = uq;
      ve = vq;
      we = wq;
    }

    CheckValue(FUNCTION_NAME, "electron direction cosines", "",
               ue*ue + ve*ve + we*we - 1.0, -1E-4, 1E-4);

    /* TTB-approximation for the Compton electron, store deposited energy */
    EdTTB = TTB(mat, part, Ee, x, y, z, ue, ve, we, wgt, t, 0, id);
  }

  CheckValue(FUNCTION_NAME, "EdTTB", "", EdTTB, 0.0, Ee);

  /***************************************************************************/


  /***** Atomic relaxation ***************************************************/

  if (Uss < RDB[DATA_PHOTON_EMIN]) {
    /* Binding energy below cutoff, energy is deposited locally */
    EdAR = Uss;
  }
  else {

    /* Pointer to atomic relaxation data */
    ptr = (long)RDB[(long)RDB[rea + REACTION_PTR_NUCLIDE] + NUCLIDE_PTR_RELAX];
    CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, ptr);

    /* Number of subshells in the atomic relaxation data */
    Nssar = (long)RDB[ptr + RELAX_NSS];

    /* Pointer to shell eletron number CDF of the relaxation model */
    ptrelncdfar = (long)RDB[ptr + RELAX_ELNCDF];
    CheckPointer(FUNCTION_NAME, "(ptrareln)", DATA_ARRAY, ptrelncdfar);

    /* Check that the relaxation and Compton electron numbers agree */
    if (Nssar <= ss)
      ss = -1;
    else if ((long)cdfeln[ss + 1] != (long)RDB[ptrelncdfar + ss + 1]) {

      /* Sample the electron number (note: cdfeln[0] = 0) */
      eli = (long)cdfeln[ss] + ceil(RandF(id)*uoccupd[ss]);

      /* Find the subshell in the relaxation data which corresponds to the
       * electron number */
      ss = -1;
      for (i = 0; i < Nssar; i++) {
        if ((long)RDB[ptrelncdfar + i + 1] >= eli) {
          ss = i;
          break;
        }
      }
    }

    if (ss == -1) {
      /* Relaxtion data not available for the subshell, energy is deposited
       * locally */
      /* NOTE: This can also happen when different binding energies are used
       * here and in AtomicRelaxation */
      EdAR = Uss;
    }
    else {
      /* Atomic relaxation */
      EdAR = AtomicRelaxation(mat, rea, part, ss, x, y, z, wgt, t, id);

      if (EdAR == -1) {
        /* Relaxtion data not available for the subshell, energy is deposited
         * locally */
        EdAR = Uss;
      }
    }
  }

  /* NOTE: AtomicRelaxation uses ENDF binding energies, whereas
   * ComptonScattering uses the ones given by Biggs. Therefore, EdAR may exceed
   * Ebss and its value is not checked.
   */
  /* CheckValue(FUNCTION_NAME, "EdAR", "", EdAR, 0.0, Ebss); */

  /***************************************************************************/

  /* Total deposited energy */
  Ed = EdTTB + EdAR;

  /* Score pulse-height detector */
  
  PulseDet(part, mat, Ed, x, y, z, wgt, id); 
}

/*****************************************************************************/


/*****************************************************************************/
void PartialPivotGauss(double A[3][3], double b[3], double x[3]) {
  /* Solves linear equation Ax=b with Gaussian elimination using scaled
   * partial pivoting.
   * NOTE: The elements of A and b are changed! */
  long i, j, k;
  long idx[3];
  double smax, aij, r, rmax, tmp, sum;
  double s[3];
  long n = 3;

  /* Calculate a scale factor array */
  for (i = 0; i < n; i++) {
    smax = 0.0;
    for (j = 0; j < n; j++) {
      aij = fabs(A[i][j]);
      if (aij > smax)
        smax = aij;
    }
    s[i] = smax;

    /* Initialize the index array */
    idx[i] = i;
  }

  /* Forward elimination */
  for (i = 0; i < n - 1; i++) {

    rmax = 0;
    k = i;

    /* Find pivot equation */
    for (j = i; j < n; j++) {
      r = fabs(A[idx[j]][i])/s[idx[j]];
      if (r > rmax) {
        rmax = r;
        k = j;
      }
    }

    /* Update the index array */
    tmp = idx[k];
    idx[k] = idx[i];
    idx[i] = tmp;

    /* Elimination */
    for (j = i+1; j < n; j++) {
      r = A[idx[j]][i]/A[idx[i]][i];
      A[idx[j]][i] = r;
      for (k = i+1; k < n; k++)
        A[idx[j]][k] = A[idx[j]][k] - r*A[idx[i]][k];

      /* Calculate b */
      b[idx[j]] = b[idx[j]] - A[idx[j]][i]*b[idx[i]];
    }
  }

  /* Back substitution */
  x[n-1] = b[idx[n-1]]/A[idx[n-1]][n-1];
  for (i = n-2; i >= 0; i--) {
    sum = b[idx[i]];
    for (j = i + 1; j < n; j++)
      sum = sum - A[idx[i]][j]*x[j];
    x[i] = sum/A[idx[i]][i];
  }

}
/*****************************************************************************/



#ifdef __cplusplus 
} 
#endif 
