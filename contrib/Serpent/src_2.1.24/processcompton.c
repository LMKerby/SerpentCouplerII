/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processcompton.c                               */
/*                                                                           */
/* Created:       2014/11/21 (TKa)                                           */
/* Last modified: 2015/05/03 (TKa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes Compton data                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessCompton:"

/* TODO: 
 * - Read Compton-profile data only if RDB[DATA_PHOTON_USE_DOPPLER] == 1
 * 
 * */

/*****************************************************************************/

void ProcessCompton(long loc0, long nuc) {

  long i, j, Ndata, ptr0, ptr1, Nss, Npz, idx, pzminidx;
  const long Npzd = 31;
  const long Nbuf = 500;  /* NOTE: same as linebuf length */
  const double pzminabs = 1.0/FS_CONST;
  double cpintinf, dummy, Elim;
  double *pzdarr, *pzarr, *xt, *Sxt, *exta;
  double **cpdmatrix, **cpmatrix, **cpintmatrix;
  char fname[MAX_STR], linebuf[500], strZ[5], tmpstr[100], *buf0;
  FILE *fp;
  const long Npzincl = 15;
  double pzincl[15] = {12., 17., 25., 35., 45., 50., 55., 65., 70., 80.,
                             90., 110., 120., 130., -1.};
  
  /* Set the minimum pz */
  pzincl[Npzincl-1] = pzminabs;
  
  /* Check pointers */
  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);


  sprintf(strZ, "%ld", (long)RDB[nuc + NUCLIDE_Z]);


  /***** Read incohrent scattering functions *****************************/

  xt = NULL;
  Sxt = NULL;

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_INCOH_FNAME));

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  while (fgets(linebuf, Nbuf, fp)) {

    if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
       (strcmp(tmpstr, strZ) == 0)) {

      if (fgets(linebuf, Nbuf, fp)) {
        if (sscanf(linebuf, "Ndata %ld\n", &Ndata) == 0)
          Die(FUNCTION_NAME, "Can't read incoherent scattering functions: Ndata");

        CheckValue(FUNCTION_NAME, "Ndata", "", Ndata, 3, INFTY);

        /* Store array size */
        WDB[loc0 + PHOTON_DIST_MINC] = Ndata;

        /* Allocate memory for momentum transfers */
        xt = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

        /* Allocate memory for scattering functions */
        Sxt = (double *)Mem(MEM_ALLOC, Ndata, sizeof(double));

        /* Read and store data */
        for (i = 0; i < Ndata; i++) {
          if (fgets(linebuf, Nbuf, fp) != NULL) {
            if (sscanf(linebuf, "%lf %lf\n", &xt[i], &Sxt[i]) != 2)
              Die(FUNCTION_NAME, "Can't read incoherent scattering functions data: xt, Sxt");
          }
          else
            Die(FUNCTION_NAME, "Can't read incoherent data");
        }

      }
      else
        Die(FUNCTION_NAME, "Can't read incoherent data");

      /* Exit loop */
      break;
    }
  }
  fclose(fp);

  if (!xt || !Sxt)
    Die(FUNCTION_NAME, "Incoherent scattering function data not found");

  /***********************************************************************/


  /***** Process incoherent scattering functions *************************/

  /* Momentum transfers */
  ptr0 = ReallocMem(DATA_ARRAY, Ndata);
  WDB[loc0 + PHOTON_DIST_PTR_VIC] = (double)ptr0;

  /* Scattering functions */
  ptr1 = ReallocMem(DATA_ARRAY, Ndata);
  WDB[loc0 + PHOTON_DIST_PTR_INC_FF] = (double)ptr1;

  /* Check the first elements */
  if (xt[0] == 0.0 || Sxt[0] == 0.0) {
    /* Don't use log */
    WDB[ptr0++] = xt[0];
    WDB[ptr1++] = Sxt[0];
  }
  else {
    WDB[ptr0++] = log(xt[0]);
    WDB[ptr1++] = log(Sxt[0]);
  }

  /* Store the rest as log-log */
  for (i = 1; i < Ndata; i++) {
    WDB[ptr0++] = log(xt[i]);
    WDB[ptr1++] = log(Sxt[i]);
  }

  /* Free memory */
  Mem(MEM_FREE, xt);
  Mem(MEM_FREE, Sxt);

  /***********************************************************************/


  /***** Read Compton profile data  **************************************/
  
  /* Check the energy limit - below the limit all the binding energies 
   * are not in descending order (at least above Z=56), and the 
   * sampling method will fail. */
  Elim = 0.0002;
  if (RDB[DATA_PHOTON_EMIN] < Elim)
    Die(FUNCTION_NAME, "Compton scattering model doesn't work below %f MeV", Elim);

  Nss = 0;

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_CP_FNAME));

  if ((fp = fopen(fname, "r")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file %s for reading", fname);

  while (fgets(linebuf, Nbuf, fp)) {
    if ((sscanf(linebuf, "#S %s %*s\n", tmpstr) == 1)
        && (strcmp(strZ, tmpstr) == 0)) {

      /* Read and store number of subshells */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#N", strlen("#N")))
          Die(FUNCTION_NAME, "Can't read compton profile data: #N missing");
        buf0 = &linebuf[strlen("#N")];
        Nss = strtod(buf0, &buf0) - 2;
        WDB[loc0 + PHOTON_DIST_N_UOCCUP] = Nss;
      }
      else
        Die(FUNCTION_NAME, "Can't read compton profile data");

      if (Nss == 0)
        Die(FUNCTION_NAME, "Couldn't find number of subshells. WTF?");

      /* Read and store number of electrons in subshells */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#UOCCUP", strlen("#UOCCUP")))
          Die(FUNCTION_NAME, "Can't read compton profile data: #UOCCUP missing");
        buf0 = &linebuf[strlen("#UOCCUP")];

        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_PTR_UOCCUP] = (double)ptr0;

        for (i = 0; i < Nss; i++)
          WDB[ptr0++] = strtod(buf0, &buf0);
      }
      else
        Die(FUNCTION_NAME, "Can't read compton profile data");


      /* Read and store the binding energies */
      /* NOTE: These binding energies differ from the ENDF binding energies! */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#UBIND", strlen("#UBIND")))
          Die(FUNCTION_NAME, "Can't read compton profile data: #UBIND missing");
        buf0 = &linebuf[strlen("#UBIND")];

        ptr0 = ReallocMem(DATA_ARRAY, Nss);
        WDB[loc0 + PHOTON_DIST_PTR_UI] = (double)ptr0;

        for (i = 0; i < Nss; i++)
          WDB[ptr0++] = strtod(buf0, &buf0) / 1.0e6; /* from eV to MeV*/
      }
      else
        Die(FUNCTION_NAME, "Can't read compton profile data");


      /* Skip the header */
      if (fgets(linebuf, Nbuf, fp)) {
        if (strncmp(linebuf, "#L", strlen("#L")))
          Die(FUNCTION_NAME, "Can't read compton profile data: #L missing");
      }
      else
        Die(FUNCTION_NAME, "Can't read compton profile data");


      /* Allocate memory for cpmatrixd and pzarrd */
      cpdmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
      for (i = 0; i < Nss; i++)
        cpdmatrix[i] = (double *)Mem(MEM_ALLOC, Npzd, sizeof(double));
      pzdarr = (double *)Mem(MEM_ALLOC, Npzd, sizeof(double));

      /* Read data to pzarr and cpmatrix */
      for (j = 0; j < Npzd; j++) {
        if (fgets(linebuf, Nbuf, fp)) {
          buf0 = linebuf;
          pzdarr[j] = strtod(buf0, &buf0);
          dummy = strtod(buf0, &buf0); /* Skip the sum of the compton profiles*/
          for (i = 0; i < Nss; i++)
            cpdmatrix[i][j] = strtod(buf0, &buf0);
        }
        else
          Die(FUNCTION_NAME, "Can't read compton profile data");
      }


      /***** Interpolation and extrapolation of the Compton profile data *****/

      Npz = Npzd + Npzincl;

      cpmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
      for (i = 0; i < Nss; i++)
        cpmatrix[i] = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));
      pzarr = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));

      /* Include interpolation and extrapolation points to pzarr  */
      for (i = 0; i < Npzd; i++)
        pzarr[i] = pzdarr[i];

      for (i = 0; i < Npzincl; i++)
        pzarr[Npzd + i] = pzincl[i];

      /* Sort the array */
      SortArray(pzarr, Npz);

      /* Check for duplicates */
      for (i = 0; i < Npz - 1; i++)
        if (pzarr[i] == pzarr[i+1])
          Die(FUNCTION_NAME, "Compton profile data pz contains duplicate values");

      for (j = 0; j < Npz; j++) {

        if (pzarr[j] < pzdarr[Npzd-1]) {
          /* Interpolate log-lin */

          if (pzarr[j] == pzarr[Npz - 1])
            idx = Npz - 2;
          else if ((idx = SearchArray(pzdarr, pzarr[j], Npzd)) == -1)
            Die(FUNCTION_NAME, "Compton profile data pz not found: %.5E", pzarr[j]);

          for (i = 0; i < Nss; i++) {
            if (pzdarr[idx] == 0.0)
              cpmatrix[i][j] = cpdmatrix[i][idx];
            else
              cpmatrix[i][j] = ENDFInterp(4, pzarr[j], pzdarr[idx], pzdarr[idx+1], cpdmatrix[i][idx],
                             cpdmatrix[i][idx+1]);
          }
        }
        else {
          /* Extrapolate log-lin */
          for (i = 0; i < Nss; i++) {
            cpmatrix[i][j] = cpdmatrix[i][Npzd-1] *
                pow(cpdmatrix[i][Npzd-1]/cpdmatrix[i][Npzd-2],
                (pzarr[j] - pzdarr[Npzd-1])/(pzdarr[Npzd-1] - pzdarr[Npzd-2]));
          }
        }
      }

      /* Check that all the values are positive */
      for (j = 0; j < Npz; j++) {
        for (i = 0; i < Nss; i++) {
          if (cpmatrix[i][j] < 0.0)
            Die(FUNCTION_NAME, "Negative value in the Compton profile data:\n cpmatrix[%ld][%ld] = %.f", i, j, cpmatrix[i][j]);
        }
      }

      /* Find the index of pzminabs if it's included in pzincl*/
      pzminidx = -1;
      for (i = 0; i < Npzincl; i++) {
        if (pzincl[i] == pzminabs) {
          for (j = 0; j < Npz; j++) {
            if (pzarr[j] == pzminabs) {
              pzminidx = j;
              break;
            }
          }
          break;
        }
      }

      if (pzminidx == -1)
        Die(FUNCTION_NAME, "pzminabs not found");

      /***********************************************************************/

      /***** Normalize and integrate Compton profiles ************************/
      
      cpintmatrix = (double **)Mem(MEM_ALLOC, Nss, sizeof(double*));
      for (i = 0; i < Nss; i++)
        cpintmatrix[i] = (double *)Mem(MEM_ALLOC, Npz, sizeof(double));
      exta = (double *)Mem(MEM_ALLOC, Nss, sizeof(double));

      for (i = 0; i < Nss; i++) {

        /* Numerical integration using trapezoidal rule */
        cpintmatrix[i][0] = 0.0;
        for (j = 1; j < Npz; j++)
          cpintmatrix[i][j] = cpintmatrix[i][j-1] + 0.5*(pzarr[j] - pzarr[j-1])
                              *(cpmatrix[i][j] + cpmatrix[i][j-1]);

        /* Extrapolation coefficient */
        exta[i] = log(cpdmatrix[i][Npzd-1]/cpdmatrix[i][Npzd-2])
            /(pzdarr[Npzd-1] - pzdarr[Npzd-2]);

        /* Check that exta is negative */
        if (exta[i] >= 0)
          Die(FUNCTION_NAME, "Positive extrapolation coefficient a in exp(a*pz)");

        /* Integral from -infty to infty */
        cpintinf = 2.0*(cpintmatrix[i][Npz-1] - cpmatrix[i][Npz-1]/exta[i]);

        /* Check the integral */
        if ((cpintinf < 0.95) || (cpintinf > 1.06))
          Die(FUNCTION_NAME, "Compton profile integral out of limits");

        /* Normalize Compton profiles and the integrals */
        for (j = 0; j < Npz; j++) {
          cpmatrix[i][j] = cpmatrix[i][j]/cpintinf;
          cpintmatrix[i][j] = cpintmatrix[i][j]/cpintinf;
        }

      }

      /***********************************************************************/

      /* Store array size */
      WDB[loc0 + PHOTON_DIST_N_INC_CP] = Npz;

      /* Store pz */
      ptr0 = ReallocMem(DATA_ARRAY, Npz);
      WDB[loc0 + PHOTON_DIST_PTR_INC_CPPZ] = (double)ptr0;
      for (j = 0; j < Npz; j++)
        WDB[ptr0++] = pzarr[j];

      /* Store maximum pz */
      WDB[loc0 + PHOTON_DIST_INC_PZMAX] = pzarr[Npz-1];

      /* Store the index of pzminabs */
      WDB[loc0 + PHOTON_DIST_INC_PZMINIDX] = pzminidx;

      /* Store compton profiles */
      ptr0 = ReallocMem(DATA_ARRAY, Nss);
      WDB[loc0 + PHOTON_DIST_PTR_INC_CP] = (double)ptr0;

      for (i = 0; i < Nss; i++) {
        ptr1 = ReallocMem(DATA_ARRAY, Npz);
        WDB[ptr0++] = (double)ptr1;
        for (j = 0; j < Npz; j++)
          WDB[ptr1++] = cpmatrix[i][j];
      }

      /* Store integrated Compton profiles (=cdf) */
      ptr0 = ReallocMem(DATA_ARRAY, Nss);
      WDB[loc0 + PHOTON_DIST_PTR_INC_CPINT] = (double)ptr0;

      for (i = 0; i < Nss; i++) {
        ptr1 = ReallocMem(DATA_ARRAY, Npz);
        WDB[ptr0++] = (double)ptr1;

        for (j = 0; j < Npz; j++)
          WDB[ptr1++] = cpintmatrix[i][j];
      }

      /* Store the extrapolation coefficients */
      ptr0 = ReallocMem(DATA_ARRAY, Nss);
      WDB[loc0 + PHOTON_DIST_PTR_INC_EXTA] = (double)ptr0;
      for (i = 0; i < Nss; i++)
        WDB[ptr0++] = exta[i];

      /* Store electron shell number cdf */
      ptr0 = ReallocMem(DATA_ARRAY, Nss + 1);
      WDB[loc0 + PHOTON_DIST_PTR_INC_ELNCDF] = (double)ptr0;
      WDB[ptr0] = 0.0;
      for (i = 0; i < Nss; i++)
        WDB[ptr0+1+i] = WDB[ptr0+i] + WDB[(long)WDB[loc0 + PHOTON_DIST_PTR_UOCCUP] + i];

      /* Check the last element */
      if ((long)RDB[nuc + NUCLIDE_Z] != (long)WDB[ptr0+Nss])
        Die(FUNCTION_NAME, "incorrect electron number cdf");

      /* Free memory */
      for (i = 0; i < Nss; i++) {
        Mem(MEM_FREE, cpdmatrix[i]);
        Mem(MEM_FREE, cpmatrix[i]);
        Mem(MEM_FREE, cpintmatrix[i]);
      }
      Mem(MEM_FREE, cpdmatrix);
      Mem(MEM_FREE, cpmatrix);
      Mem(MEM_FREE, cpintmatrix);
      Mem(MEM_FREE, pzdarr);
      Mem(MEM_FREE, pzarr);
      Mem(MEM_FREE, exta);

      /* Exit loop */
      break;

    }
  }
  fclose(fp);


  /***********************************************************************/

}
/*****************************************************************************/
