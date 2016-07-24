#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processttb.c                                   */
/*                                                                           */
/* Created:       2014/08/11 (TKa)                                           */
/* Last modified: 2015/10/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
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
static void KhiPerKappaCDF(double **, double *, double *, long, long, double **,
                           double *, long);
static void PhotonNumberYield(double **, double *, double *, long, long,
                              double *, double **);
static double PositronFp(double, long);
static double LogIntegral(const double *, const double *, long, long);
static void CumLogIntegral(const double *, const double *, double *, long, long);
static double TTBInterp(double *, double *, long, double, long);

/*****************************************************************************/


/*****************************************************************************/

void ProcessTTB()
{
  long i, j, loc0, Nrow, Ncolumn, nS, nkappa, nTekhi, Z, ptr, ptr1, nTe0, photonmat,
       mat, iso, nuc, Nnuc, Zidx;
  const long Nbuf = 5000;  /* NOTE: same as linebuf length */
  const long Zmax = 200;
  long ZidxMap[200]; /* Size given by Zmax */
  long *nTekhiData, *nSData;
  double Nmat, Fp, mfrac, afrac, Temin, Temax;
  double *kappa, *Tekhi, *TeS, *Sbre, *Sbrei, *Stote, *Te0, *Stotei, *Sbrp,
         *Stotp, *Stotpi, *Sbrpi;
  double **khie, **kpkecdf, **brecdf, **brpcdf, **brepdf, **brppdf, **Yke,
         **Ykp, **kpkpcdf, **TekhiData, **khip;
  double ***khiData, ***stoppowData;
  char fname[MAX_STR], linebuf[5000], tmpstr[100], strZ[10], *buf0;
  FILE *fp;


  /*******************************************************/
  /* Read nuclide data */

  /* Number of nuclides in photon transport mode */
  Nnuc = (long)RDB[DATA_N_PHOTON_NUCLIDES];

  /* Initialize map */
  for (i = 0; i < Zmax; i++)
    ZidxMap[i] = -1;


  /***** Read scaled bremsstrahlung xs data **********************************/

  /* Scaled bremsstrahlung data
   * - kappa
   * - Te (eV)
   * - khi (millibarn)
   * */

  fprintf(out, "Reading bremsstrahlung data...\n");

  /* Create temporary data arrays */
  khiData = (double ***)Mem(MEM_ALLOC, Nnuc, sizeof(double**));
  TekhiData = (double **)Mem(MEM_ALLOC, Nnuc, sizeof(double*));
  nTekhiData = (long *)Mem(MEM_ALLOC, Nnuc, sizeof(long));

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_ELBR_FNAME));

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for reading");

  /* Read number of kappa */
  if (!fgets(linebuf, Nbuf, fp))
    Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");
  if (sscanf(linebuf, "Kappa %ld", &nkappa) == 0)
    Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data: Kappa");

  /* Allocate memory */
  kappa = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));

  /* Read kappa */
  for (i = 0; i < nkappa; i++) {
    if (!fgets(linebuf, Nbuf, fp))
      Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");
    sscanf(linebuf, "%lf", &kappa[i]);
  }

  if (!kappa)
    Die(FUNCTION_NAME, "Bremsstrahlung kappa data couldn't be found or read");


  /* Loop over nuclides */
  Zidx = 0;
  nuc = RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR) {

    if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON) {
      Z = (long)RDB[nuc + NUCLIDE_Z];
      sprintf(strZ, "%ld", Z);

      /* Check if Z has already been read */
      if (ZidxMap[Z] != -1) {
        Die(FUNCTION_NAME, " Z=%ld already read from the bremsstrahlung data", Z);
        continue;
      }

      /* Map Z and index*/
      ZidxMap[Z] = Zidx;

      /* Rewind the file because the data might not be ordered by Z */
      rewind(fp);


      /* Find and read the data for the element */
      while (fgets(linebuf, Nbuf, fp)) {

        if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
           (strcmp(tmpstr, strZ) == 0)) {

          if (!fgets(linebuf, Nbuf, fp))
            Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");

          /* Read the data array size */
          sscanf(linebuf, "NTe %ld", &nTekhiData[Zidx]);

          /* Allocate memory for data arrays */
          TekhiData[Zidx] = (double *)Mem(MEM_ALLOC, nTekhiData[Zidx], sizeof(double));
          khiData[Zidx] = (double **)Mem(MEM_ALLOC, nTekhiData[Zidx], sizeof(double*));

          for (i = 0; i < nTekhiData[Zidx]; i++) {
            if (!fgets(linebuf, Nbuf, fp))
              Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");

            sscanf(linebuf, "Te %lf", &TekhiData[Zidx][i]);
            TekhiData[Zidx][i] = 1.0e-6*TekhiData[Zidx][i];  /* eV to MeV */

            /* Allocate memory */
            khiData[Zidx][i] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));

            for (j = 0; j < nkappa; j++) {
              if (!fgets(linebuf, Nbuf, fp))
                Die(FUNCTION_NAME, "Can't read bremsstrahlung xs data");
              sscanf(linebuf, "%lf", &khiData[Zidx][i][j]);
              khiData[Zidx][i][j] = 1.0e-3*khiData[Zidx][i][j];   /* millibarn to barn */
            }

          }
          if (!khiData[Zidx][0])
            Die(FUNCTION_NAME, "Bremsstrahlung xs data couldn't be found or read for Z=%s", strZ);

          /* Exit loop */
          break;
        }
      }

      if (!khiData[Zidx] || !TekhiData[Zidx])
        Die(FUNCTION_NAME, "Bremsstrahlung xs data couldn't be found or read for Z=%s", strZ);
    }

    Zidx++;
    nuc = NextItem(nuc);
  }

  /* Close the file */
  fclose(fp);

  /* Check the number of read nuclides */
  if (Zidx != Nnuc)
    Die(FUNCTION_NAME, "Difference in number of nuclides, photon nuclides "
                       "= %ld, bremsstrahlung data nuclides = %ld", Nnuc, Zidx);

  /* Check that the energy arrays are identical */
  for (Zidx = 1; Zidx < Nnuc; Zidx++) {
    if (nTekhiData[0] != nTekhiData[Zidx])
      Die(FUNCTION_NAME, "Difference in bremsstrahlung energy array sizes");

    for (i = 0; i < nTekhiData[0]; i++) {
      if (TekhiData[0][i] != TekhiData[Zidx][i])
        Die(FUNCTION_NAME, "Difference in bremsstrahlung energy array");
    }
  }

  /* Set single energy array */
  nTekhi = nTekhiData[0];
  Tekhi = (double *)Mem(MEM_ALLOC, nTekhi, sizeof(double));
  for (i = 0; i < nTekhi; i++)
    Tekhi[i] = TekhiData[0][i];

  fprintf(out, "OK.\n\n");

  /***************************************************************************/


  /***** Read stopping power data ********************************************/

  /* Stopping power data S:
   * - Electron kinetic energy (MeV)
   * - Stopping power S (MeV cm2/g)
   *
   * Columns: Energy, collision S, radiative S, total S, density effect parameter.
   * */

  fprintf(out, "Reading electron stopping power data...\n");

  /* Create temporary data arrays */
  stoppowData = (double ***)Mem(MEM_ALLOC, Nnuc, sizeof(double**));
  nSData = (long *)Mem(MEM_ALLOC, Nnuc, sizeof(long));

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_ELSP_FNAME));

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for reading");

  /* Loop over nuclides */
  Zidx = 0;
  nuc = RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR) {

    if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON) {
      Z = (long)RDB[nuc + NUCLIDE_Z];
      sprintf(strZ, "%ld", Z);

      /* Check that the map is correct */
      if (ZidxMap[Z] != Zidx)
        Die(FUNCTION_NAME, "Error in Zmap for Z=%ld, index=%ld", Z, Zidx);

      /* Rewind the file because the data might not be ordered by Z */
      rewind(fp);

      while (fgets(linebuf, Nbuf, fp)) {

        if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
           (strcmp(tmpstr, strZ) == 0)) {

          if (!fgets(linebuf, Nbuf, fp))
            Die(FUNCTION_NAME, "Can't read stopping power data");

          /* Read the data array size */
          sscanf(linebuf, "Ndata %ld %ld", &Nrow, &Ncolumn);
          nSData[Zidx] = Nrow;

          /* Allocate memory for data array */
          stoppowData[Zidx] = (double **)Mem(MEM_ALLOC, Nrow, sizeof(double*));

          for (i = 0; i < Nrow; i++) {
            if (!fgets(linebuf, Nbuf, fp))
              Die(FUNCTION_NAME, "Can't read stopping power data");

            stoppowData[Zidx][i] = (double *)Mem(MEM_ALLOC, Ncolumn, sizeof(double));

            buf0 = linebuf;
            for (j = 0; j < Ncolumn; j++)
              stoppowData[Zidx][i][j] = strtod(buf0, &buf0);
          }

          if (!stoppowData[Zidx][0])
            Die(FUNCTION_NAME, "Stopping power data couldn't be found or read for Z=%s", strZ);

          /* Exit loop */
          break;
        }
      }

      if (!stoppowData[Zidx])
        Die(FUNCTION_NAME, "Stopping power data couldn't be found or read for Z=%s", strZ);
    }

    Zidx++;
    nuc = NextItem(nuc);
  }

  /* Close the file */
  fclose(fp);

  /* Check the number of read nuclides */
  if (Zidx != Nnuc)
    Die(FUNCTION_NAME, "Difference in number of nuclides, photon nuclides "
                       "= %ld, stopping power data nuclides = %ld", Nnuc, Zidx);

  /* Check that the energy arrays of stopping powers are identical */
  for (Zidx = 1; Zidx < Nnuc; Zidx++) {
    if (nSData[0] != nSData[Zidx])
      Die(FUNCTION_NAME, "Difference in stopping power data array sizes");

    for (i = 0; i < nSData[0]; i++) {
      if (stoppowData[0][i][0] != stoppowData[Zidx][i][0])
        Die(FUNCTION_NAME, "Difference in topping power energy array");
    }
  }

  /* Set single energy array */
  nS = nSData[0];
  TeS = (double *)Mem(MEM_ALLOC, nS, sizeof(double));
  for (i = 0; i < nS; i++)
    TeS[i] = stoppowData[0][i][0];

  /* Free memory */
  Mem(MEM_FREE, nSData);

  fprintf(out, "OK.\n\n");

  /***************************************************************************/


  /***** Check energy boundaries**********************************************/

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

  /* Create electron energy array
   * NOTE: TTB() assumes that the energy array is uniformly distributed in log.
   * TTB() electron/positron energy search must be changed if uniform log is
   * not used.
   */
  Temin = RDB[DATA_PHOTON_EMIN];
  Temax = RDB[DATA_PHOTON_EMAX];
  Te0 = MakeArray(Temin, Temax, nTe0, 2);

  /***************************************************************************/


  /***** Allocate memory *****************************************************/

  Sbre = (double *)Mem(MEM_ALLOC, nS, sizeof(double));
  Stote = (double *)Mem(MEM_ALLOC, nS, sizeof(double));
  Stotei = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  Sbrei = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  Sbrp = (double *)Mem(MEM_ALLOC, nS, sizeof(double));
  Stotp = (double *)Mem(MEM_ALLOC, nS, sizeof(double));
  Stotpi = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));
  Sbrpi = (double *)Mem(MEM_ALLOC, nTe0, sizeof(double));

  kpkecdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  kpkpcdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brecdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brpcdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brepdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  brppdf = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  Yke = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));
  Ykp = (double **)Mem(MEM_ALLOC, nTe0, sizeof(double*));

  for (i = 0; i < nTe0; i++) {
    kpkecdf[i] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
    kpkpcdf[i] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
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

  khie = (double **)Mem(MEM_ALLOC, nTekhi, sizeof(double*));
  khip = (double **)Mem(MEM_ALLOC, nTekhi, sizeof(double*));
  for (i = 0; i < nTekhi; i++) {
    khie[i] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
    khip[i] = (double *)Mem(MEM_ALLOC, nkappa, sizeof(double));
  }

  /***************************************************************************/

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR) {

    /* Skip parent materials (JLe 1.10.2015 / 2.1.25) */

    if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
      {
	/* Pointer to next */

	mat = NextItem(mat);

	/* Cycle loop */

	continue;
      }
    
    /* Initialize stopping power arrays */
    for (i = 0; i < nS; i++) {
      Sbre[i] = 0.0;
      Stote[i] = 0.0;
      Sbrp[i] = 0.0;
      Stotp[i] = 0.0;
    }

    /* Initialize scaled bremsstrahlung matrix */
    for (i = 0; i < nTekhi; i++) {
      for (j = 0; j < nkappa; j++) {
        khie[i][j] = 0.0;
        khip[i][j] = 0.0;
      }
    }

    /* Set atomic density */
    Nmat = RDB[mat + MATERIAL_ADENS];

    /* Loop over composition */
    photonmat = 0;
    iso = (long)RDB[mat + MATERIAL_PTR_COMP];
    while (iso > VALID_PTR) {

      nuc = RDB[iso + COMPOSITION_PTR_NUCLIDE];

      /* Only photon nuclides are included! */
      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON) {

        /* Set flag for NUCLIDE_TYPE_PHOTON */
        photonmat = 1;

        /* Mass fraction */
        mfrac = RDB[nuc + NUCLIDE_AW]*RDB[iso + COMPOSITION_ADENS]
                       /(N_AVOGADRO*RDB[mat + MATERIAL_MDENS]);

        /* Atomic fraction */
        afrac = RDB[iso + COMPOSITION_ADENS]/RDB[mat + MATERIAL_ADENS];

        /* Z and map index */
        Z = (long)RDB[nuc + NUCLIDE_Z];
        Zidx = ZidxMap[Z];

        /* Calculate radiative and total stopping powers using Bragg's rule.
         * For compounds, this seems to hold for radiative stopping powers, at
         * least when compared to NIST data. NOTE However, it doesn't hold for
         * collisional stopping powers (tot=rad+col). */
        /* TODO: Density effect correction here? */
        for (i = 0; i < nS; i++) {
          Sbre[i] += mfrac*stoppowData[Zidx][i][2];
          Stote[i] += mfrac*stoppowData[Zidx][i][3];
          Fp = PositronFp(TeS[i], Z);
          Sbrp[i] += mfrac*Fp*stoppowData[Zidx][i][2];
          Stotp[i] += mfrac*(stoppowData[Zidx][i][3] + (Fp - 1.0)*stoppowData[Zidx][i][2]);
        }

        /* Scaled bremsstrahlung data */
        for (i = 0; i < nTekhi; i++) {
          Fp = PositronFp(Tekhi[i], Z);
          for (j = 0; j < nkappa; j++) {
            khie[i][j] += Nmat*afrac*Z*Z*khiData[Zidx][i][j];
            khip[i][j] += Fp*Nmat*afrac*Z*Z*khiData[Zidx][i][j];
          }
        }
      }

      iso = NextItem(iso);
    }

    /* Check that material has photon nuclides */
    if (!photonmat) {
      mat = NextItem(mat);
      continue;
    }

    /* Multiply stopping powers with material density (g/cm^3)*/
    for (i = 0; i < nS; i++) {
      Sbre[i] *= RDB[mat + MATERIAL_MDENS];
      Stote[i] *= RDB[mat + MATERIAL_MDENS];
      Sbrp[i] *= RDB[mat + MATERIAL_MDENS];
      Stotp[i] *= RDB[mat + MATERIAL_MDENS];
    }

    /***** TESTING: Numerically integrate the radiative stopping powers
     * using the bremsstrahlung data  */
    /*
    if (0) {
    double integral;
    for (i = 0; i < nTekhi; i++) {
      double integral = CumLogIntegral(kappa, khi[i], kpkcdf[i], nkappa, 0);
      integral = 0.0;

      integral = LogIntegral(kappa, khi[i], nkappa, 0); // log-log
      integral = LogIntegral(kappa, khi[i], nkappa, 1); // log-lin

      for (j = 1; j < nkappa; j++) {
        integral += 0.5*(khi[i][j] + khi[i][j-1])*(kappa[j] - kappa[j-1]);
      }

      integral = CSplineIntegrate0(kappa, khie[i], nkappa, 0.0, 0.0, kappa[0], kappa[nkappa-1]);

      beta = sqrt(Tekhi[i]*(Tekhi[i] + 2*E_RESTMASS))/(Tekhi[i] + E_RESTMASS);
      Sbrc[i] = Nmat/(beta*beta)*Tekhi[i]*integral;
    }
    }
    */
    /*************************************************************************/

    /*************************************************************************/

    /* Log-log interpolate total stopping power */
    CSplineInterpolate0(TeS, Stote, nS, 0, 0, Te0, Stotei, nTe0, 2);
    CSplineInterpolate0(TeS, Stotp, nS, 0, 0, Te0, Stotpi, nTe0, 2);

    /* Log-log interpolate radiative stopping power */
    CSplineInterpolate0(TeS, Sbre, nS, 0, 0, Te0, Sbrei, nTe0, 2);
    CSplineInterpolate0(TeS, Sbrp, nS, 0, 0, Te0, Sbrpi, nTe0, 2);

    /* Calculate the CDF of khi per kappa (log-log) (also interpolates the
     * transpose of khi).*/
    KhiPerKappaCDF(khie, Tekhi, kappa, nTekhi, nkappa, kpkecdf, Te0, nTe0);
    KhiPerKappaCDF(khip, Tekhi, kappa, nTekhi, nkappa, kpkpcdf, Te0, nTe0);

    /* Calculate photon number yields */
    PhotonNumberYield(kpkecdf, Te0, kappa, nTe0, nkappa, Stotei, Yke);
    PhotonNumberYield(kpkpcdf, Te0, kappa, nTe0, nkappa, Stotpi, Ykp);

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


    /***** Set data **********************************************************/

    /* Allocate memory for TTB block */
    loc0 = NewItem(mat + MATERIAL_PTR_TTB, TTB_BLOCK_SIZE);

    WDB[loc0 + TTB_NE] = nTe0;

    /* energy array */
    ptr = ReallocMem(DATA_ARRAY, nTe0);
    WDB[loc0 + TTB_E] = (double)ptr;
    for (i = 0; i < nTe0; i++)
      WDB[ptr++] = Te0[i];

    /* log energy array */
    ptr = ReallocMem(DATA_ARRAY, nTe0);
    WDB[loc0 + TTB_LE] = (double)ptr;
    for (i = 0; i < nTe0; i++)
      WDB[ptr++] = log(Te0[i]);

    /* bremsstrahlung cdf for electrons */
    ptr = ReallocMem(DATA_ARRAY, nTe0);
    WDB[loc0 + TTB_BRECDF] = (double)ptr;
    for (i = 0; i < nTe0; i++) {
      ptr1 = ReallocMem(DATA_ARRAY, i+1);
      WDB[ptr++] = (double)ptr1;
      for (j = 0; j < i+1; j++)
        WDB[ptr1++] = brecdf[i][j];
    }

    /* bremsstrahlung pdf for electrons */
    ptr = ReallocMem(DATA_ARRAY, nTe0);
    WDB[loc0 + TTB_BREPDF] = (double)ptr;
    for (i = 0; i < nTe0; i++) {
      ptr1 = ReallocMem(DATA_ARRAY, nTe0);
      WDB[ptr++] = (double)ptr1;
      for (j = 0; j < i+1; j++)
        WDB[ptr1++] = brepdf[i][j];
    }

    /* log photon number yield for electrons above the first energy */
    ptr = ReallocMem(DATA_ARRAY, nTe0);
    WDB[loc0 + TTB_LYE] = (double)ptr;
    for (i = 0; i < nTe0; i++)
      WDB[ptr++] = log(Yke[0][i]);


    /* Set separate bremsstrahlung data for positrons */
    if ((long)RDB[DATA_PHOTON_TTBPM] == YES) {

      /* bremsstrahlung cdf for positrons */
      ptr = ReallocMem(DATA_ARRAY, nTe0);
      WDB[loc0 + TTB_BRPCDF] = (double)ptr;
      for (i = 0; i < nTe0; i++) {
        ptr1 = ReallocMem(DATA_ARRAY, i+1);
        WDB[ptr++] = (double)ptr1;
        for (j = 0; j < i+1; j++)
          WDB[ptr1++] = brpcdf[i][j];
      }

      /* bremsstrahlung pdf for positrons */
      ptr = ReallocMem(DATA_ARRAY, nTe0);
      WDB[loc0 + TTB_BRPPDF] = (double)ptr;
      for (i = 0; i < nTe0; i++) {
        ptr1 = ReallocMem(DATA_ARRAY, nTe0);
        WDB[ptr++] = (double)ptr1;
        for (j = 0; j < i+1; j++)
          WDB[ptr1++] = brppdf[i][j];
      }

      /* log photon number yield for positrons above the first energy */
      ptr = ReallocMem(DATA_ARRAY, nTe0);
      WDB[loc0 + TTB_LYP] = (double)ptr;
      for (i = 0; i < nTe0; i++)
        WDB[ptr++] = log(Ykp[0][i]);
    }

    /*************************************************************************/

    mat = NextItem(mat);
  }

  /***************************************************************************/


  /***** Free allocated memory ***********************************************/

  for (j = 0; j < Nnuc; j++) {
    Mem(MEM_FREE, TekhiData[j]);

    for (i = 0; i < nTekhi; i++)
      Mem(MEM_FREE, khiData[j][i]);
    Mem(MEM_FREE, khiData[j]);

    for (i = 0; i < nS; i++)
      Mem(MEM_FREE, stoppowData[j][i]);
    Mem(MEM_FREE, stoppowData[j]);
  }
  Mem(MEM_FREE, TekhiData);
  Mem(MEM_FREE, khiData);
  Mem(MEM_FREE, stoppowData);
  Mem(MEM_FREE, nTekhiData);

  Mem(MEM_FREE, Tekhi);
  Mem(MEM_FREE, kappa);
  Mem(MEM_FREE, Te0);
  Mem(MEM_FREE, TeS);
  Mem(MEM_FREE, Sbre);
  Mem(MEM_FREE, Sbrei);
  Mem(MEM_FREE, Stote);
  Mem(MEM_FREE, Stotei);
  Mem(MEM_FREE, Sbrp);
  Mem(MEM_FREE, Sbrpi);
  Mem(MEM_FREE, Stotp);
  Mem(MEM_FREE, Stotpi);

  for (i = 0; i < nTekhi; i++) {
    Mem(MEM_FREE, khie[i]);
    Mem(MEM_FREE, khip[i]);
  }
  Mem(MEM_FREE, khie);
  Mem(MEM_FREE, khip);

  for (i = 0; i < nTe0; i++) {
    Mem(MEM_FREE, kpkecdf[i]);
    Mem(MEM_FREE, kpkpcdf[i]);
    Mem(MEM_FREE, Yke[i]);
    Mem(MEM_FREE, Ykp[i]);
    Mem(MEM_FREE, brecdf[i]);
    Mem(MEM_FREE, brpcdf[i]);
    Mem(MEM_FREE, brepdf[i]);
    Mem(MEM_FREE, brppdf[i]);
  }

  Mem(MEM_FREE, kpkecdf);
  Mem(MEM_FREE, kpkpcdf);
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

static void PhotonNumberYield(double **kpkcdf, double *Te, double *kappa,
                              long NTe, long Nkappa, double *Stot, double **Yk) {
  /* Integrates photon number yields */

  long i0, i;
  double kappacr, kpkcdfcr, beta, imfp;
  double *imfppStot, **coeff;

  /* Allocate memory */
  imfppStot = (double *)Mem(MEM_ALLOC, NTe, sizeof(double));
  coeff = (double **)Mem(MEM_ALLOC, NTe, sizeof(double*));

  for (i = 0; i < NTe; i++)
    coeff[i] = (double *)Mem(MEM_ALLOC, 4, sizeof(double));


  /* Loop over energies */
  for (i0 = 0; i0 < NTe-2; i0++) {

    /* Calculate inverse mean free path divided by total stopping power for
     * bremsstrahlung above Te0[i0] */
    for (i = i0; i < NTe; i++) {

      kappacr = Te[i0]/Te[i];

      if (kappacr == 1.0) {
        imfppStot[i] = 0.0;
      }
      else {

        /* Interpolate the cdf of khi per kappa at kappacr */
        kpkcdfcr = TTBInterp(kappa, kpkcdf[i], Nkappa, kappacr, 3);

        /* Inverse mean free path above kappacr for bremsstrahlung */
        beta = sqrt(Te[i]*(Te[i] + 2*E_RESTMASS))/(Te[i] + E_RESTMASS);
        imfp = (kpkcdf[i][Nkappa-1] - kpkcdfcr)/(beta*beta);

        /* Inverse mean free path per total stopping power */
        imfppStot[i] = imfp/Stot[i];
      }
    }

    /* Calculate the photon number yield for electrons*/
    CSplineConstruct(&Te[i0], &imfppStot[i0], NTe-i0, 0.0, 0.0, coeff, 1,
                     &Yk[i0][i0 + 1]);
  }

  /* Free allocated memory */
  for (i = 0; i < NTe; i++)
    Mem(MEM_FREE, coeff[i]);

  Mem(MEM_FREE, coeff);
  Mem(MEM_FREE, imfppStot);
}

/*****************************************************************************/


/*****************************************************************************/

static void KhiPerKappaCDF(double **khi, double *Tekhi, double *kappa,
                             long NTekhi, long Nkappa, double **kpkcdf, double *Tei,
                             long NTei) {
  /* Calculates the cumulative distribution function (CDF) of khi per kappa.
   * Also interpolates the transpose of khi.
   *
   * Array sizes:
   * khi = khi[NTekhi][Nkappa]
   * kpkcdf = kpkcdf[NTei][Nkappa]
   * */
  long i, j;
  double *kpk, **khiT, **khiTi;

  /* Allocate memory */
  kpk = (double *)Mem(MEM_ALLOC, Nkappa, sizeof(double));
  khiT = (double **)Mem(MEM_ALLOC, Nkappa, sizeof(double*));
  khiTi = (double **)Mem(MEM_ALLOC, Nkappa, sizeof(double*));

  for (i = 0; i < Nkappa; i++) {
    khiT[i] = (double *)Mem(MEM_ALLOC, NTekhi, sizeof(double));
    khiTi[i] = (double *)Mem(MEM_ALLOC, NTei, sizeof(double));
  }

  /* Set transpose of khi */
  for (i = 0; i < Nkappa; i++)
    for (j = 0; j < NTekhi; j++)
      khiT[i][j] = khi[j][i];

  /* Interpolate the transpose of khi (log-log interpolation using splines) */
  for (i = 0; i < Nkappa; i++)
    CSplineInterpolate0(Tekhi, khiT[i], NTekhi, 0.0, 0.0, Tei, khiTi[i], NTei, 2);

  for (i = 0; i < NTei; i++) {
    /* Calculate khii per kappa */
    for (j = 0; j < Nkappa; j++)
      kpk[j] = khiTi[j][i]/kappa[j];

    /* Calculate cdf of khi per kappa (log-log) */
    CumLogIntegral(kappa, kpk, kpkcdf[i], Nkappa, 0);
  }

  /* Free allocated memory */
  for (i = 0; i < Nkappa; i++) {
    Mem(MEM_FREE, khiT[i]);
    Mem(MEM_FREE, khiTi[i]);
  }
  Mem(MEM_FREE, khiT);
  Mem(MEM_FREE, khiTi);
  Mem(MEM_FREE, kpk);
}

/*****************************************************************************/


/*****************************************************************************/

static double PositronFp(double Te, long Z) {
  /* Calculates the approximative ratio of the radiative stopping
   * powers for positrons and electrons for a given element Z and positron
   * kinetic energy Te.
   *
   * Source: F. Salvat et al., "PENELOPE-2011: A Code System for Monte Carlo
   * Simulation of Electron and Photon Transport" (2011) */
  double tp, Fp;
  tp = log(1.0 + 1.0e6/(Z*Z)*Te/E_RESTMASS);
  Fp = 1.0 - exp(tp*(-1.2359e-1 + tp*(6.1274e-2 + tp*(-3.1516e-2 + tp*(7.7446e-3
              + tp*(-1.0595e-3 + tp*(7.0568e-5 - 1.8080e-6*tp)))))));
  return Fp;
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
#ifdef __cplusplus 
} 
#endif 
