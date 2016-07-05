/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processphotoelectric.c                         */
/*                                                                           */
/* Created:       2014/11/21 (TKa)                                           */
/* Last modified: 2015/07/03 (TKa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Processes photoelectric data                                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessPhotoelectric:"

/*****************************************************************************/

/* NOTE:
 * - The minimum energy limit (DATA_PHOTON_EMIN) has not been implemented
 *
 * */

void ProcessPhotoelectric(long loc0, long nuc) {

  long ptrNxs, ptrU, ptrE, ptrE2, ptrxs, ptrxs2;
  double *Etmparr, *xstmparr;
  long i, j, Ndata, Nss, Nread, continueloop, MTnumber;
  const long Nbuf = 500;  /* NOTE: same as linebuf length */
  double Etmp, xstmp;
  char fname[MAX_STR], linebuf[500], strZ[5], tmpstr[100];
  FILE *fp;

  /* Check pointers */
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

  /* Initialize variables */
  Nss = 0;

  sprintf(strZ, "%ld", (long)RDB[nuc + NUCLIDE_Z]);

  
  /***** Read photoelectric shell xs data **********************************/

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_PESS_FNAME));

  if (!(fp = fopen(fname, "r")))
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  while (fgets(linebuf, Nbuf, fp)){

    if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
       (strcmp(tmpstr, strZ) == 0)) {

      if (!fgets(linebuf, Nbuf, fp))
        Die(FUNCTION_NAME, "Can't read photoelectric data, file: %s , line:\n %s", fname, linebuf);

      if (sscanf(linebuf, "NSS %ld\n", &Nss) == 0)
        Die(FUNCTION_NAME, "Can't read photoelectric data: NSS");

      CheckValue(FUNCTION_NAME, "NSS", "", Nss, 0, INFTY);

      /* Store number of subshells (only above 1 keV!) */
      WDB[loc0 + PHOTON_DIST_N_PE_SS] = Nss;

      if (Nss > 0) {

        /* Subshell data array sizes */
        ptrNxs = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_PTR_PE_NXS] = (double)ptrNxs;

        /* Binding energy array */
        ptrU = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_PTR_PE_EB] = (double)ptrU;

        /* Photon energy array */
        ptrE = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_PTR_PE_SSE] = (double)ptrE;

        /* Subshell xs array */
        ptrxs = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_PTR_PE_SSXS] = (double)ptrxs;

        /* Loop shells */
	for (i = 0; i < Nss; i++) {

          /* Read MT number */
          if (!fgets(linebuf, Nbuf, fp))
            Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n line: %s", fname, linebuf);
          if (sscanf(linebuf, "MT %ld\n", &MTnumber) == 0)
            Die(FUNCTION_NAME, "Can't read photoelectric data: MT");

          CheckValue(FUNCTION_NAME, "MTnumber", "", MTnumber, 534, 599);

          /* Read number of data */
          if (!fgets(linebuf, Nbuf, fp))
            Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n line: %s", fname, linebuf);
          if (sscanf(linebuf, "Ndata %ld\n", &Ndata) == 0)
            Die(FUNCTION_NAME, "Can't read photoelectric data: Ndata");

          CheckValue(FUNCTION_NAME, "MTnumber", "", Ndata, 3, INFTY);

          /* Temporary arrays */
          Etmparr = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
          xstmparr = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

          Nread = 0;
          continueloop = 0;

          /* Read data */
          for (j = 0; j < Ndata; j++) {

            if (!fgets(linebuf, Nbuf, fp))
              Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n line: %s", fname, linebuf);

            if (continueloop)
              continue;

            if (sscanf(linebuf, "%lf %lf\n", &Etmp, &xstmp) != 2)
              Die(FUNCTION_NAME, "Can't read photoelectric data: Etmp, xstmp");

            Etmp = Etmp/1.0e6;  /* From ev to MeV */
            Etmparr[Nread] = Etmp;
            xstmparr[Nread] = xstmp;
            Nread++;

            /* Maximum limit, NOTE: assuming the data is in ascending
             * order */
            if (Etmp > RDB[DATA_PHOTON_EMAX])
              continueloop = 1;
          }

          if (Nread == 0)
            Die(FUNCTION_NAME, "No photoelectric shell data below DATA_PHOTON_EMAX");


          /* Set the ionisation energy, which is assumed to be the first
           * element of the data */
          WDB[ptrU++] = Etmparr[0];

          /* Store array size */
          WDB[ptrNxs++] = Nread;

          /* Photon energy array */
          ptrE2 = ReallocMem(DATA_ARRAY, Nread);
          WDB[ptrE++] = (double)ptrE2;

          /* Subshell xs */
          ptrxs2 = ReallocMem(DATA_ARRAY, Nread);
          WDB[ptrxs++] = (double)ptrxs2;

          /* Store data as log-log */
          for (j = 0; j < Nread; j++) {
            /* Check zeros TODO: tähän jotain fiksumpaa? */
            if (Etmparr[j] == 0.0)
              Die(FUNCTION_NAME, "Zero photon energy in photoelectric subshell cross section data. Log-log interpolation not possible.");
            if (Etmparr[j] == 0.0)
              Die(FUNCTION_NAME, "Zero cross section in photoelectric subshell cross section data. Log-log interpolation not possible.");

            WDB[ptrE2++] = log(Etmparr[j]);
            WDB[ptrxs2++] = log(xstmparr[j]);
          }

          Mem(MEM_FREE, Etmparr);
          Mem(MEM_FREE, xstmparr);

        }
      }

      /* Exit loop */
      break;

    }
  }

  fclose(fp);

  /*************************************************************************/


  /***** Read total photoelectric xs ***************************************/

  /* Here the the total xs is the sum of the interpolated xs read above */

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_PETOT_FNAME));

  if (!(fp = fopen(fname, "r")))
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  while (fgets(linebuf, Nbuf, fp)){

    if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
       (strcmp(tmpstr, strZ) == 0)) {

      if (!fgets(linebuf, Nbuf, fp))
        Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n line: %s", fname, linebuf);

      if (sscanf(linebuf, "Ndata %ld\n", &Ndata) == 0)
        Die(FUNCTION_NAME, "Can't read photoelectric data: Ndata");

      CheckValue(FUNCTION_NAME, "Ndata", "", Ndata, 0, INFTY);

      /* Temporary arrays */
      Etmparr = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));
      xstmparr = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

      Nread = 0;
      continueloop = 0;

      /* Read data */
      for (j = 0; j < Ndata; j++) {

        if (!fgets(linebuf, Nbuf, fp))
          Die(FUNCTION_NAME, "Can't read photoelectric data \n file: %s \n line: %s", fname, linebuf);

        if (continueloop)
          continue;

        if (sscanf(linebuf, "%lf %lf\n", &Etmp, &xstmp) != 2)
          Die(FUNCTION_NAME, "Can't read photoelectric data: Etmp, xstmp");

        Etmp = Etmp/1.0e6;    /* From ev to MeV */
        Etmparr[Nread] = Etmp;
        xstmparr[Nread] = xstmp;
        Nread++;

        /*TODO: Check that E > RDB[DATA_PHOTON_EMIN]*/

        /* Maximum limit, NOTE: assuming the data is in ascending
         * order */
        if (Etmp > RDB[DATA_PHOTON_EMAX])
          continueloop = 1;
      }

      if (Nread == 0)
        Die(FUNCTION_NAME, "No photoelectric xs data below DATA_PHOTON_EMAX");

      /* Store number of data */
      WDB[loc0 + PHOTON_DIST_N_PE_TOT] = Nread;

      /* Energy pointer */
      ptrE = ReallocMem(DATA_ARRAY, Nread);
      WDB[loc0 + PHOTON_DIST_PTR_PE_ETOT] = (double)ptrE;

      /* xs pointer */
      ptrxs = ReallocMem(DATA_ARRAY, Nread);
      WDB[loc0 + PHOTON_DIST_PTR_PE_XSTOT] = (double)ptrxs;

      /* Store data as log-log */
      for (j = 0; j < Nread; j++) {
        /* Check zeros TODO: tähän jotain fiksumpaa varoitusta? */
        if (Etmparr[j] == 0.0)
          Die(FUNCTION_NAME, "Zero photon energy in photoelectric subshell cross section data. Log-log interpolation not possible.");
        if (xstmparr[j] == 0.0)
          Die(FUNCTION_NAME, "Zero cross section in photoelectric subshell cross section data. Log-log interpolation not possible.");

        WDB[ptrE++] = log(Etmparr[j]);
        WDB[ptrxs++] = log(xstmparr[j]);
      }

      Mem(MEM_FREE, Etmparr);
      Mem(MEM_FREE, xstmparr);

      /* Exit loop */
      break;

    }
  }

  fclose(fp);

  /*************************************************************************/

}
/*****************************************************************************/


















