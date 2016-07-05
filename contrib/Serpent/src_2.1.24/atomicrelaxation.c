/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : atomicrelaxation.c                             */
/*                                                                           */
/* Created:       2014/07/12 (TKa)                                           */
/* Last modified: 2015/05/03 (TKa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Handles atomic relaxtion process                             */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AtomicRelaxation:"

/*****************************************************************************/

/* TODO:
 *  - return locally deposited energy or not?
 *  - variables ebi and ntr are not needed
 * */


double AtomicRelaxation(long rea, long part, long ss0, double x, double y,
                   double z, double wgt, double t, long id) {

  long ptd, ss, ssj, ssk, subj, subk, tr, newp, ssqueuefront, ssqueueback,
      ptretr, ptrftr, ptrsubj, ptrsubk, Nssd;
  double ebi0, etr, Edeploc, Etrtot, u, v, w;
  const double *d2imap, *ntrd, *ebid, *elnd;
  long ssqueue[111];  /* NOTE: the length comes from the maximum number of electrons  */
  long vacancy[PHOTON_NSS_MAX] = {0};   /* NOTE: the length comes from the maximum number of subshells */

  /* Check the initial vacancy subshell index ss0 */
  if (ss0 < 0) {
    Die(FUNCTION_NAME, "Negative primary subshell index");
  }

  /* Pointer to photon reaction data */
  ptd = (long)RDB[rea + REACTION_PTR_PHOTON_DIST];
  CheckPointer(FUNCTION_NAME, "(ptd)", DATA_ARRAY, ptd);

  /* Subshell data */
  Nssd = (long)RDB[ptd + PHOTON_DIST_N_AR_NSS];
  CheckValue(FUNCTION_NAME, "Nssd", "", Nssd, 1, 40);

  if (ss0 >= Nssd) {
    /* Subshell data not available */
    /* TODO: Energy deposited in this case? */
    return -1;
  }

  d2imap = &RDB[(long)RDB[ptd + PHOTON_DIST_PTR_AR_D2IMAP]];
  CheckPointer(FUNCTION_NAME, "(AR_D2IMAP)", DATA_ARRAY, (long)RDB[ptd + PHOTON_DIST_PTR_AR_D2IMAP]);

  ntrd = &RDB[(long)RDB[ptd + PHOTON_DIST_PTR_AR_NTR]];
  CheckPointer(FUNCTION_NAME, "(AR_NTR)", DATA_ARRAY, (long)RDB[ptd + PHOTON_DIST_PTR_AR_NTR]);

  ebid = &RDB[(long)RDB[ptd + PHOTON_DIST_PTR_EBI]];
  CheckPointer(FUNCTION_NAME, "(EBI)", DATA_ARRAY, (long)RDB[ptd + PHOTON_DIST_PTR_EBI]);

  elnd = &RDB[(long)RDB[ptd + PHOTON_DIST_PTR_AR_ELN]];
  CheckPointer(FUNCTION_NAME, "(AR_ELN)", DATA_ARRAY, (long)RDB[ptd + PHOTON_DIST_PTR_AR_ELN]);

  /* Pointers to subshell data */
  ptretr = (long)RDB[ptd + PHOTON_DIST_PTR_AR_ETR];
  CheckPointer(FUNCTION_NAME, "(AR_ETR)", DATA_ARRAY, ptretr);

  ptrftr = (long)RDB[ptd + PHOTON_DIST_PTR_AR_FTR];
  CheckPointer(FUNCTION_NAME, "(AR_FTR)", DATA_ARRAY, ptrftr);

  ptrsubj = (long)RDB[ptd + PHOTON_DIST_PTR_AR_SUBJ];
  CheckPointer(FUNCTION_NAME, "(AR_SUBJ)", DATA_ARRAY, ptrsubj);

  ptrsubk = (long)RDB[ptd + PHOTON_DIST_PTR_AR_SUBK];
  CheckPointer(FUNCTION_NAME, "(AR_SUBK)", DATA_ARRAY, ptrsubk);


  /* Binding energy of the initial vacancy subshell */
  ebi0 = ebid[ss0];

  /* Check the lower energy limit */
  if (ebi0 < RDB[DATA_PHOTON_EMIN])
    return ebi0;

  /* Set the initial vacancy */
  vacancy[ss0] = 1;

  /* Intialize subshell queue */
  ssqueuefront = 0;
  ssqueueback = 0;
  ssqueue[ssqueueback] = ss0;

  /* Total transported energy */
  Etrtot = 0.0;

  while (ssqueuefront <= ssqueueback) {

    ss = ssqueue[ssqueuefront++];

    /* Check the subshell index */
    if (ss == -1) {
      /* Data not available for the subshell (ebi below DATA_PHOTON_EMIN) */
      /* NOTE: Vacancy array is not changed */
      continue;
    }

    /* Decrease primary subshell vacancy */
    vacancy[ss]--;

    /* Check the binding energy */
    if (ebid[ss] < RDB[DATA_PHOTON_EMIN])
      continue;

    while (1) {

      /* Sample transition */
      tr = SearchArray(&RDB[(long)RDB[ptrftr + ss]], RandF(id), ntrd[ss] + 1);

      CheckValue(FUNCTION_NAME, "tr", "", tr, 0, ntrd[ss] - 1);

      subj = RDB[(long)RDB[ptrsubj + ss] + tr];
      subk = RDB[(long)RDB[ptrsubk + ss] + tr];

      ssj = d2imap[subj];
      ssk = d2imap[subk];

      /* Check vacancies */
      /* NOTE: Vacancies are checked only from available subshells */
      if ((ssj != -1) && (vacancy[ssj] == elnd[ssj]))
        continue;
      else if ((ssk != -1) && (vacancy[ssk] == elnd[ssk]))
        continue;
      else
        break;
    }

    /* Transition energy */
    etr = RDB[(long)RDB[ptretr + ss] + tr];

    if (ssj >= 0) {
      /* Increase secondary subshell vacancy */
      vacancy[ssj]++;

      /* Put secondary subshell to queue */
      ssqueue[++ssqueueback] = ssj;

      CheckValue(FUNCTION_NAME, "ssqueueback", "Maximum number of vacancies exceeded", ssqueueback, ssqueuefront, 111);
    }

    if (subk == 0) {
      /* Radiative transition */

      if (etr > RDB[DATA_PHOTON_EMIN]) {

        /* Create a new photon */
        newp = DuplicateParticle(part, id);

        /* Put variables */
        WDB[newp + PARTICLE_X] = x;
        WDB[newp + PARTICLE_Y] = y;
        WDB[newp + PARTICLE_Z] = z;
        WDB[newp + PARTICLE_WGT] = wgt; /*TODO: Onko oikein? */
        WDB[newp + PARTICLE_T] = t;

        /* Put energy */
        WDB[newp + PARTICLE_E] = etr;

        /* Sample direction isotropically */
        IsotropicDirection(&u, &v, &w, id);

        /* Put direction cosines */
        WDB[newp + PARTICLE_U] = u;
        WDB[newp + PARTICLE_V] = v;
        WDB[newp + PARTICLE_W] = w;

        /* Put photon in queue */
        ToQue(newp, id);

        /* Transported energy */
        Etrtot += etr;

      }

    }
    else {
      /* Non-radiative transition */

      if (ssk >= 0) {
        /* Increase tertiary subshell vacancy */
        vacancy[ssk]++;

        /* Put tertiary subshell to queue */
        ssqueue[++ssqueueback] = ssk;

        CheckValue(FUNCTION_NAME, "ssqueueback", "Maximum number of vacancies exceeded", ssqueueback, ssqueuefront, 111);
      }

      if (etr > RDB[DATA_PHOTON_EMIN]) {

        /* Sample electron direction isotropically */
        IsotropicDirection(&u, &v, &w, id);

        /* Use TTB-approximation for the electron */
        TTB(rea, part, etr, x, y, z, u, v, w, wgt, t, 0, id);

        /* Transported energy */
        Etrtot += etr;

      }
    }
  }

  /* Locally deposit energy */
  Edeploc = ebi0 - Etrtot;

  return Edeploc;

}

/*****************************************************************************/
