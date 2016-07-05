/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processrelaxation.c                            */
/*                                                                           */
/* Created:       2014/11/23 (TKa)                                           */
/* Last modified: 2015/07/03 (TKa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Processes atomic relaxation data                             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessRelaxation:"


/*****************************************************************************/

void ProcessRelaxation() {

  long i, j, elementfound, Nss, subi, ntr, eln, maxdesignator, ptrd2i, ptrntr,
      ptrebi, ptreln, ptretr, ptrftr, ptrsubi, ptrsubj, ptrsubk, ptrsubj2,
      ptrsubk2, ptretr2, ptrftr2, ptrelncdf, nuc, loc0;
  double ebi, Etrtmp, ftrtmp, ftrsum;
  const long Nbuf = 500;    /* NOTE: same as linebuf length */
  char fname[MAX_STR], linebuf[500], strZ[5], tmpstr[100];
  FILE *fp;

  /* Set file name */
  strcpy(fname, GetText(DATA_PHOTON_DATA_DIR));
  strcat(fname, GetText(DATA_PHOTON_RELAX_FNAME));

  /* Open file */
  if (!(fp = fopen(fname, "r")))
    Die(FUNCTION_NAME, "Unable to open file for reading");

  /* Loop over nuclides */
  nuc = RDB[DATA_PTR_NUC0];
  while (nuc > VALID_PTR) {

    if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON) {

      /* Rewind the file because the data might not be ordered by Z */
      rewind(fp);

      /* Set Z string */
      sprintf(strZ, "%ld", (long)RDB[nuc + NUCLIDE_Z]);

      elementfound = 0;

      while (fgets(linebuf, Nbuf, fp)) {

        if ((sscanf(linebuf, "Element %s\n", tmpstr) == 1) &&
           (strcmp(tmpstr, strZ) == 0)) {

          elementfound = 1;

          /* Allocate memory for atomic relaxation block */
          loc0 = NewItem(nuc + NUCLIDE_PTR_RELAX, RELAX_BLOCK_SIZE);

          if (!fgets(linebuf, Nbuf, fp))
            Die(FUNCTION_NAME, "Can't read relaxation data");

          /* Read and store number of subshells */
          if (sscanf(linebuf, "NSS %ld\n", &Nss) == 0)
            Die(FUNCTION_NAME, "Can't read relaxation data NSS");
          CheckValue(FUNCTION_NAME, "Nss", "", Nss, 0, 30);
          WDB[loc0 + RELAX_NSS] = Nss;

          /* Create a map from subshell designator to subshell index */
          /* (Designator from ENDF, index 0,1,2 etc) */
          maxdesignator = 599 - 534 + 1;
          ptrd2i = ReallocMem(DATA_ARRAY, maxdesignator + 1);
          WDB[loc0 + RELAX_D2IMAP] = ptrd2i;

          /* Initialize map */
          for (i = 0; i < maxdesignator + 1; i++)
            WDB[ptrd2i + i] = -1.0;

          /* Set maximum number of vacancies */
          /*WDB[loc0 + PHOTON_DIST_AR_MAXVAC] = RDB[nuc + NUCLIDE_Z];*/

          /* Read and init subshell data */
          if (Nss > 0) {

            /* Subshell data array size */
            ptrntr = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_NTR] = (double)ptrntr;

            /* Binding energy array */
            ptrebi = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_EBI] = (double)ptrebi;

            /* Number of electrons in subshell */
            ptreln = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_ELN] = (double)ptreln;

            /* Transition energy */
            ptretr = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_ETR] = (double)ptretr;

            /* Transition cumulative probability function (discrete) */
            ptrftr = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_FTR] = (double)ptrftr;

            /* Subshell designator array */
            ptrsubi = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_SUBI] = (double)ptrsubi;

            /* Secondary subshell designator array */
            ptrsubj = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_SUBJ] = (double)ptrsubj;

            /* Tertiary subshell designator array */
            ptrsubk = ReallocMem(DATA_ARRAY, Nss);
            WDB[loc0 + RELAX_SUBK] = (double)ptrsubk;

            /* Electron number CDF */
            ptrelncdf = ReallocMem(DATA_ARRAY, Nss+1);
            WDB[loc0 + RELAX_ELNCDF] = (double)ptrelncdf;
            WDB[ptrelncdf] = 0.0;

            /* Loop shells */
            for (i = 0; i < Nss; i++) {

              /* Read and store SUBI*/
              if (fgets(linebuf, Nbuf, fp)) {
                if (sscanf(linebuf, "SUBI %ld\n", &subi) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data SUBI");

                /* Check subi*/
                CheckValue(FUNCTION_NAME, "subi", "", subi, 1, maxdesignator);

                WDB[ptrsubi++] = subi;

                /* Map SUBI to index */
                WDB[ptrd2i + subi] = i;
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Read number of transitions */
              if (fgets(linebuf, Nbuf, fp)) {
                if (sscanf(linebuf, "NTR %ld\n", &ntr) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data NTR");
                CheckValue(FUNCTION_NAME, "ntr", "", ntr, 0, INFTY);
                WDB[ptrntr++] = ntr;
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Read electron binding energy */
              if (fgets(linebuf, Nbuf, fp)) {
                if (sscanf(linebuf, "EBI %lf\n", &ebi) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data EBI");
                CheckValue(FUNCTION_NAME, "ebi", "", ebi, 0, INFTY);
                WDB[ptrebi++] = ebi/1.0e6; /* ev to MeV */
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Read number of electrons in subshell */
              if (fgets(linebuf, Nbuf, fp)) {
                if (sscanf(linebuf, "ELN %ld\n", &eln) == 0)
                  Die(FUNCTION_NAME, "Can't read relaxation data ELN");
                CheckValue(FUNCTION_NAME, "eln", "", eln, 2, 6);
                WDB[ptreln++] = eln;
                WDB[ptrelncdf+i+1] = RDB[ptrelncdf+i] + eln;
              }
              else
                Die(FUNCTION_NAME, "Can't read relaxation data");

              /* Secondary subshell designator array */
              ptrsubj2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptrsubj++] = (double)ptrsubj2;

              /* Tertiary subshell designator array */
              ptrsubk2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptrsubk++] = (double)ptrsubk2;

              /* Transition energy */
              ptretr2 = ReallocMem(DATA_ARRAY, ntr);
              WDB[ptretr++] = (double)ptretr2;

              /* Transition cumulative probability function (discrete) */
              ptrftr2 = ReallocMem(DATA_ARRAY, ntr + 1);
              WDB[ptrftr++] = (double)ptrftr2;

              /* The first element of the cdf array is zero */
              WDB[ptrftr2++] = 0.0;

              /* Read transitions */
              ftrsum = 0.0;
              for (j = 0; j < ntr; j++) {

                if (fgets(linebuf, Nbuf, fp) != NULL) {
                  if (sscanf(linebuf, "%lf %lf %lf %lf\n", &WDB[ptrsubj2++],
                      &WDB[ptrsubk2++], &Etrtmp, &ftrtmp) != 4)
                    Die(FUNCTION_NAME, "Can't read relaxation data, line:\n %s", linebuf);

                  /* TODO: CheckValue */

                  WDB[ptretr2++] = Etrtmp/1.0e6; /* From ev to MeV */

                  ftrsum += ftrtmp;
                  WDB[ptrftr2+j] = ftrsum;
                }
                else
                  Die(FUNCTION_NAME, "Can't read relaxation data");
              }

              /* Normalize ftr */
              for (j = 0; j < ntr; j++)
                WDB[ptrftr2+j] /= ftrsum;

            }
          }

          /* Exit loop */
          break;
        }
      }

      if (!elementfound)
        Die(FUNCTION_NAME, "Atomic relaxation data not found for element Z=%s", strZ);
    }

    nuc = NextItem(nuc);
  }

  /* Close the file */
  fclose(fp);

  /*************************************************************************/

}
