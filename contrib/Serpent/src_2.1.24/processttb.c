/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processttb.c                                   */
/*                                                                           */
/* Created:       2014/08/11 (TKa)                                           */
/* Last modified: 2015/05/21 (TKa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes thick target bremsstrahlung data                   */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessTTB:"


/* Local function definitions */
static double LogIntegral(const double *, const double *, long, long);
static void CumLogIntegral(const double *, const double *, double *, long, long);
static double TTBInterp(double *, double *, long, double x1, long);

/*****************************************************************************/


/*****************************************************************************/

void ProcessTTB(long loc0, long nuc)
{
  long i, i0, j, Nrow, Ncolumn, nS, nkappa, nTekhi, Z, ptr, ptr1, nTe0;
  const long Nbuf = 500;  /* NOTE: same as linebuf length */
  double Nmat, imfp, kappacr, kpkcdfcr, beta, tp, Fp;
  double *kappa, *Tekhi, *TeS, *Sbr, *Sbri, *Stot, *Te0, *kpk, *imfppStote,
      *imfppStotp, *Stoti;
  double **stoppow, **khi, **kpkcdf, **khiT, **khiiT, **coeff, **brecdf,
      **brpcdf, **brepdf, **brppdf, **Yke, **Ykp;
  char fname[MAX_STR], linebuf[500], tmpstr[100], strZ[10], *buf0;
  FILE *fp;

  
  /* Check pointers */
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Initialize stuff */
  Z = (long)RDB[nuc + NUCLIDE_Z];
  sprintf(strZ, "%ld", Z);
  Nmat = N_AVOGADRO/RDB[nuc + NUCLIDE_AW];

  /* Avoid compiler warnings */
  kappa = NULL;
  Tekhi = NULL;
  stoppow = NULL;
  khi = NULL;


  /***** Read stopping power data ********************************************/

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_ELSP_FNAME));

  nS = 0;

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for reading");

  while (fgets(linebuf, Nbuf, fp)) {

    if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
       (strcmp(tmpstr, strZ) == 0)) {

      if (!fgets(linebuf, Nbuf, fp))
        Die(FUNCTION_NAME, "Can't read stopping power data");

      /* Read the data array size */
      sscanf(linebuf, "Ndata %ld %ld", &Nrow, &Ncolumn);

      nS = Nrow;

      /* Allocate memory for data array */
      stoppow = (double **)Mem(MEM_ALLOC, Nrow, sizeof(double*));

      for (i = 0; i < Nrow; i++) {
        if (!fgets(linebuf, Nbuf, fp))
          Die(FUNCTION_NAME, "Can't read stopping power data");

        stoppow[i] = (double *)Mem(MEM_ALLOC, Ncolumn, sizeof(double));

        buf0 = linebuf;
        for (j = 0; j < Ncolumn; j++)
          stoppow[i][j] = strtod(buf0, &buf0);
      }

      /* Exit loop */
      break;
    }
  }

  fclose(fp);

  if (!stoppow)
    Die(FUNCTION_NAME, "Stopping power data couldn't be found or read for Z=%s", strZ);

  /***************************************************************************/


  /***** Read scaled bremsstrahlung xs data **********************************/

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_ELBR_FNAME));

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for reading");

  /* Read kappa first */
  if (!fgets(linebuf, Nbuf, fp))
    Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");
  if (sscanf(linebuf, "Kappa %ld", &nkappa) == 0)
    Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data: Kappa");

  /* Allocate memory for data array */
  kappa = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));

  for (i = 0; i < nkappa; i++) {
    if (!fgets(linebuf, Nbuf, fp))
      Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");
    sscanf(linebuf, "%lf", &kappa[i]);
  }

  /* Find and read the data for the element */
  while (fgets(linebuf, Nbuf, fp)) {

    if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
       (strcmp(tmpstr, strZ) == 0)) {

      if (!fgets(linebuf, Nbuf, fp))
        Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");

      /* Read the data array size */
      sscanf(linebuf, "NTe %ld", &nTekhi);

      /* Allocate memory for data arrays */
      Tekhi = (double *)Mem(MEM_ALLOC, nTekhi, sizeof(double));
      khi = (double **)Mem(MEM_ALLOC, nTekhi, sizeof(double*));

      for (i = 0; i < nTekhi; i++) {
        if (!fgets(linebuf, Nbuf, fp))
          Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");

        sscanf(linebuf, "Te %lf", &Tekhi[i]);

        /* Allocate memory */
        khi[i] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));

        for (j = 0; j < nkappa; j++) {
          if (!fgets(linebuf, Nbuf, fp))
            Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");
          sscanf(linebuf, "%lf", &khi[i][j]);
        }
      }

      /* Exit loop */
      break;
    }
  }

  fclose(fp);

  if (!khi || !kappa || !Tekhi)
    Die(FUNCTION_NAME, "Bremsstrahlung xs data couldn't be found or read for Z=%s", strZ);

  /***************************************************************************/


  /***** Data processing *****************************************************/

  /* Stopping power data S:
   * - Electron kinetic energy (MeV)
   * - Stopping power (MeV cm2/g)
   *
  * Columns: Energy, collision S, radiative S, total S, density effect
   * parameter.
   * */
  TeS = (double *)Mem(MEM_ALLOC, nS, sizeof(double));
  Sbr = (double *)Mem(MEM_ALLOC, nS, sizeof(double));
  Stot = (double *)Mem(MEM_ALLOC, nS, sizeof(double));

  for (i = 0; i < nS; i++) {
    TeS[i] = stoppow[i][0];
    Sbr[i] = stoppow[i][2];
    Stot[i] = stoppow[i][3];
  }

  /* Scaled bremsstrahlung data
   * - kappa
   * - Te (eV)
   * - khi (millibarn)
   * */
  for (i = 0; i < nTekhi; i++) {
    Tekhi[i] = 1.0e-6*Tekhi[i];  /* eV to MeV */
    for (j = 0; j < nkappa; j++) {
      khi[i][j] = 1e-3*khi[i][j]; /* millibarn to barn */
    }
  }

  /***************************************************************************/


  /***** Check data **********************************************************/

  if (Tekhi[0] > RDB[DATA_PHOTON_EMIN])
    Die(FUNCTION_NAME, "Minimum electron energy in TTB xs data larger than DATA_PHOTON_EMIN");

  if (Tekhi[nTekhi - 1] < RDB[DATA_PHOTON_EMAX])
    Die(FUNCTION_NAME, "Maximum electron energy in TTB xs data smaller than DATA_PHOTON_EMAX");


  if (TeS[0] > RDB[DATA_PHOTON_EMIN])
    Die(FUNCTION_NAME, "Minimum electron energy in TTB stopping power data larger than DATA_PHOTON_EMIN");

  if (TeS[nS - 1] < RDB[DATA_PHOTON_EMAX])
    Die(FUNCTION_NAME, "Maximum electron energy in TTB stopping power data smaller than DATA_PHOTON_EMAX");

  /***************************************************************************/


  /***** Create energy grid **************************************************/

  /* Length of the energy grid */
  nTe0 = 200;

  /* Create electron energy array */
  Te0 = MakeArray(RDB[DATA_PHOTON_EMIN], RDB[DATA_PHOTON_EMAX], nTe0, 2);

  /***************************************************************************/


  /***** Allocate memory *****************************************************/

  imfppStote = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  imfppStotp = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  Stoti = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  Sbri = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  kpkcdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  coeff = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brecdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brpcdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brepdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brppdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  Yke = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  Ykp = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));

  for (i = 0; i < nTe0; i++) {
    kpkcdf[i] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
    coeff[i] = (double *)Mem(MEM_ALLOC, 4, sizeof(double));
    brecdf[i] = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
    brpcdf[i] = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
    brepdf[i] = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
    brppdf[i] = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
    Yke[i] = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
    Ykp[i] = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));

    for (j = 0; j < nTe0; j++) {
      Yke[i][j] = 0;
      Ykp[i][j] = 0;
    }
  }

  kpk = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
  khiT = (double **)Mem(MEM_ALLOC, nkappa, sizeof(double*));
  khiiT = (double **)Mem(MEM_ALLOC, nkappa, sizeof(double*));

  for (i = 0; i < nkappa; i++) {
    khiT[i] = (double *)Mem(MEM_ALLOC, nTekhi, sizeof(double));
    khiiT[i] = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  }

  /***************************************************************************/


  /* Set transpose of khi */
  for (i = 0; i < nkappa; i++)
    for (j = 0; j < nTekhi; j++)
      khiT[i][j] = khi[j][i];

  /* Interpolate the transpose of khi (log-log interpolation using splines) */
  for (i = 0; i < nkappa; i++)
    CSplineInterpolate0(Tekhi, khiT[i], nTekhi, 0.0, 0.0, Te0, khiiT[i], nTe0, 2);


  for (i = 0; i < nTe0; i++) {

    /* Calculate khii per kappa */
    for (j = 0; j < nkappa; j++)
      kpk[j] = khiiT[j][i]/kappa[j];

    /* Calculate cdf of khi per kappa (log-log) */
    CumLogIntegral(kappa, kpk, kpkcdf[i], nkappa, 0);
  }

  /* Log-log interpolate total stopping power */
  CSplineInterpolate0(TeS, Stot, nS, 0, 0, Te0, Stoti, nTe0, 2);

  /* Log-log interpolate radiative stopping power */
  CSplineInterpolate0(TeS, Sbr, nS, 0, 0, Te0, Sbri, nTe0, 2);


  /***** Integrate photon yields *********************************************/

  for (i0 = 0; i0 < nTe0-2; i0++) {

    /* Calculate inverse mean free path divided by total stopping power for
     * bremsstrahlung above Te0[i0] */

    for (i = i0; i < nTe0; i++) {

      kappacr = Te0[i0]/Te0[i];

      if (kappacr == 1.0) {
        imfppStote[i] = 0.0;
        imfppStotp[i] = 0.0;
      }
      else {
        /* Interpolate the cdf of khi per kappa at kappacr */
        kpkcdfcr = TTBInterp(kappa, kpkcdf[i], nkappa, kappacr, 3);

        beta = sqrt(Te0[i]*(Te0[i] + 2*E_RESTMASS))/(Te0[i] + E_RESTMASS);

        /* Inverse mean free path*/
        imfp = Nmat*Z*Z/(beta*beta)*(kpkcdf[i][nkappa-1] - kpkcdfcr);

        imfppStote[i] = imfp/Stoti[i];

        /* Positrons */
        tp = log(1.0 + 1.0e6/(Z*Z)*Te0[i]/E_RESTMASS);
        Fp = 1.0 - exp(tp*(-1.2359e-1 + tp*(6.1274e-2 + tp*(-3.1516e-2 + tp*(7.7446e-3
                    + tp*(-1.0595e-3 + tp*(7.0568e-5 - 1.8080e-6*tp)))))));
        imfppStotp[i] = Fp*imfp/(Stoti[i] + Sbri[i]*(Fp - 1.0));
      }

    }

    /* Calculate the photon number yield for electrons*/
    CSplineConstruct(&Te0[i0], &imfppStote[i0], nTe0-i0, 0.0, 0.0, coeff, 1,
                     &Yke[i0][i0 + 1]);

    /* Calculate the photon number yield for positrons*/
    CSplineConstruct(&Te0[i0], &imfppStotp[i0], nTe0-i0, 0.0, 0.0, coeff, 1,
                     &Ykp[i0][i0 + 1]);
  }

  /* Calculate bremsstrahlung cdf for electrons and positrons */
  for (i = 0; i < nTe0; i++) {
    for (j = 0; j < i+1; j++) {
      brecdf[i][j] = fabs(Yke[j][i] - Yke[0][i]);
      brpcdf[i][j] = fabs(Ykp[j][i] - Ykp[0][i]);
    }
  }

  /* The bremsstrahlung pdf is constructed 'linearly', i.e. assuming linear
   * interpolation between the data points. The last value of the pdf is set to
   * zero, meaning that the probability P(Te = Ek) = 0. The linear assumption
   * violates the cspline integration above, but the introduced error should be
   * small. */
  for (i = 0; i < nTe0; i++) {
    brepdf[i][i] = 0.0;
    brppdf[i][i] = 0.0;
    for (j = i-1; j >= 0; j--) {
      brepdf[i][j] = 2.0*(brecdf[i][j+1] - brecdf[i][j])/(Te0[j+1] - Te0[j])
          - brepdf[i][j+1];
      brppdf[i][j] = 2.0*(brpcdf[i][j+1] - brpcdf[i][j])/(Te0[j+1] - Te0[j])
          - brppdf[i][j+1];
    }
  }

  /* Set the yield of the first energy to be non-zero in order to avoid
   * log(0) = Infty */
  Yke[0][0] = 1.e-2*Yke[0][1];
  Ykp[0][0] = 1.e-2*Ykp[0][1];

  /***************************************************************************/


  /***** Set data ************************************************************/

  WDB[loc0 + PHOTON_DIST_N_TTB_NE] = nTe0;

  /* energy array */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_E] = (double)ptr;
  for (i = 0; i < nTe0; i++)
    WDB[ptr++] = Te0[i];

  /* log energy array */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_LE] = (double)ptr;
  for (i = 0; i < nTe0; i++)
    WDB[ptr++] = log(Te0[i]);

  /* bremsstrahlung cdf for electrons */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_BRECDF] = (double)ptr;
  for (i = 0; i < nTe0; i++) {
    ptr1 = ReallocMem(DATA_ARRAY, i+1);
    WDB[ptr++] = (double)ptr1;
    for (j = 0; j < i+1; j++)
      WDB[ptr1++] = brecdf[i][j];
  }

  /* bremsstrahlung cdf for positrons */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_BRPCDF] = (double)ptr;
  for (i = 0; i < nTe0; i++) {
    ptr1 = ReallocMem(DATA_ARRAY, i+1);
    WDB[ptr++] = (double)ptr1;
    for (j = 0; j < i+1; j++)
      WDB[ptr1++] = brpcdf[i][j];
  }

  /* bremsstrahlung pdf for electrons */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_BREPDF] = (double)ptr;
  for (i = 0; i < nTe0; i++) {
    ptr1 = ReallocMem(DATA_ARRAY, nTe0);
    WDB[ptr++] = (double)ptr1;
    for (j = 0; j < i+1; j++)
      WDB[ptr1++] = brepdf[i][j];
  }

  /* bremsstrahlung pdf for positrons */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_BRPPDF] = (double)ptr;
  for (i = 0; i < nTe0; i++) {
    ptr1 = ReallocMem(DATA_ARRAY, nTe0);
    WDB[ptr++] = (double)ptr1;
    for (j = 0; j < i+1; j++)
      WDB[ptr1++] = brppdf[i][j];
  }

  /* log photon number yield for electrons above the first energy */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_LYE] = (double)ptr;
  for (i = 0; i < nTe0; i++)
    WDB[ptr++] = log(Yke[0][i]);

  /* log photon number yield for positrons above the first energy */
  ptr = ReallocMem(DATA_ARRAY, nTe0);
  WDB[loc0 + PHOTON_DIST_PTR_TTB_LYP] = (double)ptr;
  for (i = 0; i < nTe0; i++)
    WDB[ptr++] = log(Ykp[0][i]);

  /***************************************************************************/


  /***** Free allocated memory ***********************************************/

  for (i = 0; i < Nrow; i++)
    Mem(MEM_FREE, stoppow[i]);
  Mem(MEM_FREE, stoppow);

  Mem(MEM_FREE, Tekhi);
  Mem(MEM_FREE, kappa);
  Mem(MEM_FREE, Te0);
  Mem(MEM_FREE, TeS);
  Mem(MEM_FREE, Sbr);
  Mem(MEM_FREE, Sbri);
  Mem(MEM_FREE, Stot);
  Mem(MEM_FREE, Stoti);
  Mem(MEM_FREE, imfppStote);
  Mem(MEM_FREE, imfppStotp);
  Mem(MEM_FREE, kpk);

  for (i = 0; i < nTekhi; i++) {
    Mem(MEM_FREE, khi[i]);
  }
  Mem(MEM_FREE, khi);

  for (i = 0; i < nkappa; i++) {
    Mem(MEM_FREE, khiT[i]);
    Mem(MEM_FREE, khiiT[i]);
  }
  Mem(MEM_FREE, khiT);
  Mem(MEM_FREE, khiiT);

  for (i = 0; i < nTe0; i++) {
    Mem(MEM_FREE, kpkcdf[i]);
    Mem(MEM_FREE, coeff[i]);
    Mem(MEM_FREE, Yke[i]);
    Mem(MEM_FREE, Ykp[i]);
    Mem(MEM_FREE, brecdf[i]);
    Mem(MEM_FREE, brpcdf[i]);
    Mem(MEM_FREE, brepdf[i]);
    Mem(MEM_FREE, brppdf[i]);
  }

  Mem(MEM_FREE, kpkcdf);
  Mem(MEM_FREE, coeff);
  Mem(MEM_FREE, Yke);
  Mem(MEM_FREE, Ykp);
  Mem(MEM_FREE, brecdf);
  Mem(MEM_FREE, brpcdf);
  Mem(MEM_FREE, brepdf);
  Mem(MEM_FREE, brppdf);

  /***************************************************************************/

}

/*****************************************************************************/


/*****************************************************************************/

static double LogIntegral(const double *x, const double *y, long N, long mode) {
  /* Helper function */

  long ninteg, i;
  double x1, x2, y1, y2, a, F;

  if (N < 2)
    Die(FUNCTION_NAME, "LogIntegral: Number of elements less than two");

  ninteg = N - 1;
  F = 0;

  if (mode == 0) {
    /* Integral of log-log interpolation on linear scale */
    for (i = 0; i < ninteg; i++) {
      x1 = x[i];
      x2 = x[i+1];
      y1 = y[i];
      y2 = y[i+1];


      if (x1 == x2)
         Die(FUNCTION_NAME, "LogIntegral: x1 == x2");
      if (x1 == 0)
        Die(FUNCTION_NAME, "LogIntegral: x1 == 0");
      if (x2 == 0)
        Die(FUNCTION_NAME, "LogIntegral: x2 == 0");
      if (y1 == 0)
        Die(FUNCTION_NAME, "LogIntegral: y1 == 0");
      if (y2 == 0)
        Die(FUNCTION_NAME, "LogIntegral: y2 == 0");

      a = log(y2/y1)/log(x2/x1);

      if (a == -1)
        F = y1*x1*log(x2/x1);
      else
        F = F + y1/(a + 1.0)*(pow(x2/x1, a)*x2 - x1);
    }
  }
  else if (mode == 1) {
    /* Integral of log-lin interpolation on linear scale */
    for (i = 0; i < ninteg; i++) {
      x1 = x[i];
      x2 = x[i+1];
      y1 = y[i];
      y2 = y[i+1];

      if (x1 == x2)
          Die(FUNCTION_NAME, "LogIntegral: x1 == x2");
      if (x1 == 0)
        Die(FUNCTION_NAME, "LogIntegral: x1 == 0");
      if (x2 == 0)
        Die(FUNCTION_NAME, "LogIntegral: x2 == 0");

      F = F + y1*(x2 - x1) + (y2 - y1)*(x2 + (x1 - x2)/(log(x2/x1)));
    }
  }
  else
    Die(FUNCTION_NAME, "LogIntegral: unknown rule");

  return F;

}

/*****************************************************************************/


/*****************************************************************************/

static void CumLogIntegral(const double *x, const double *y, double *F,
                             long N, long mode) {
  /* Helper function */

  long i;
  double xtmp[2];
  double ytmp[2];

  if (N < 2)
    Die(FUNCTION_NAME, "CumLogIntegral: Number of elements less than two");

  F[0] = 0;

  for (i = 1; i < N; i++) {
    xtmp[0] = x[i-1];
    xtmp[1] = x[i];
    ytmp[0] = y[i-1];
    ytmp[1] = y[i];
    F[i] = F[i-1] + LogIntegral(xtmp, ytmp, 2, mode);
  }

}

/*****************************************************************************/


/*****************************************************************************/

static double TTBInterp(double *x, double *y, long n, double x1, long type) {
  /* Helper function */

  long lo;

  lo = SearchArray(x, x1, n);

  if (lo == -1) {
    if (x[n-1] == x1)
      lo = n-2;
    else {
      printf("x[0]: %.5E\n", x[0]);
      printf("y[0]: %.5E\n", y[0]);
      printf("x[n-1]: %.5E\n", x[n-1]);
      printf("y[n-1]: %.5E\n", y[n-1]);
      printf("a: %.5E\n", x1);
      Die(FUNCTION_NAME, "TTBInterp: \"a\" out of boundaries of x");
    }
  }

  /* TODO: pitäiskö tää ottaa pois? */
  if (x[lo] == x[lo+1])
    Die(FUNCTION_NAME, "TTBInterp: x[lo] == x[lo+1]");

  return ENDFInterp(type, x1, x[lo],x[lo+1], y[lo], y[lo+1]);

}

/*****************************************************************************/
