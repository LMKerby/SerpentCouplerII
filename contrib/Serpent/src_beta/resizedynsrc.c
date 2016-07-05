/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : resizedynsrc.c                                 */
/*                                                                           */
/* Created:       2015/18/05 (VVa)                                           */
/* Last modified: 2016/02/17 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Resizes live source population to PRECDET_N_LIVE at time     */
/*              interval boundaries                                          */
/*                                                                           */
/* Comments:   -Number to normalize is calculated based on live weight and   */
/*              weight to emit on the upcoming interval.                     */
/*             -Live weight is calculated in countdynsrc.c                   */
/*             -Weight to emit is calculated in decayprecdet.c               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ResizeDynSrc:"

/*****************************************************************************/

void ResizeDynSrc()
{
  long loc0, ptr, part, id, N, m;
  long nbatch, nsrc, np;
  long Nlive, Nemit;
  double wgt, wgt0, P, mul;
  double Wlive, Wemit, avewgt;
  

  /***************************************************************************/

  /* Get pointer to precursor detector */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) < VALID_PTR)
    return;

#ifdef DNPRINT
  fprintf(out, "resizedynsrc.c -->\n");  
#endif

  /* Get live weight at interval boundary */

  Wlive = RDB[loc0 + PRECDET_W_LIVE];

  /* Get weight to emit over next interval */

  Wemit = RDB[loc0 + PRECDET_W_EMIT];

#ifdef DNPRINT
  fprintf(out, "Wlive %E\n", Wlive);
  fprintf(out, "Wemit %E\n", Wemit);
#endif
  
  /* Get wanted criticality population */

  nbatch = (long)RDB[DATA_SRC_POP];

  /* Calculate number of live neutrons and neutrons to emit */

  Nlive = (long)round(nbatch*Wlive/(Wlive + Wemit));
  Nemit = (long)round(nbatch*Wemit/(Wlive + Wemit));
  
#ifdef DNPRINT
  fprintf(out, "Nlive %ld\n", Nlive);
  fprintf(out, "Nemit %ld\n", Nemit);
#endif
  
  /* Store numbers */

  WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
  WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

  /* Calculate average weight */

  avewgt = (Wemit + Wlive)/NormCoef(PARTICLE_TYPE_NEUTRON)/RDB[DATA_SRC_POP];

  /* Store average weight */

  WDB[loc0 + PRECDET_W_AVE] = avewgt;

  /* Reset total weight and source size */

  wgt0 = 0.0;
  nsrc = 0;

  /* Loop over threads */
  
  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to dummy */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      while((ptr = NextItem(ptr)) > VALID_PTR)
	{
	  /* Check type */
	  
	  if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
	    Die(FUNCTION_NAME, "Invalid particle type");
	  
	  /* Add to total weight and source size */
	  
	  wgt0 = wgt0 + RDB[ptr + PARTICLE_WGT];
	  nsrc = nsrc + 1;

	}
    }

#ifdef DNPRINT
  fprintf(out, "%ld particles in banks\n", nsrc);
#endif

  /***** Population control **************************************************/

  /* Reset weight and number of particles */

  wgt = wgt0;
  np = nsrc;

  /* Compare population size to number of live neutrons */

  if (nsrc < Nlive)
    {
      /* Calculate multiplication */

      P = ((double)Nlive)/((double)nsrc);
      mul = (long)P;
      P = P - (double)mul;

      /* Loop over threads */
  
      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	{
	  /* Get particles from bank */

	  /* Get pointer to bank */

	  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Get pointer to last item */
	  /* We have to loop backwards since new neutrons are  */
	  /* added to the end of the bank and we don't want to */
	  /* multiply them */

	  ptr = LastItem(ptr);
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  while(ptr > VALID_PTR)
	    {
	      /* Check type */
	  
	      if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
		{
		  if ((long)RDB[ptr + PARTICLE_TYPE] == PARTICLE_TYPE_DUMMY)
		    break;
		  else
		    Die(FUNCTION_NAME, "Invalid particle type");
		}

	      /* Sample multiplication */

	      if (RandF(0) < P)
		N = mul;
	      else
		N = mul - 1;

	      /* Loop over multiplication */

	      for (m = 0; m < N; m++)
		{
		  /* Duplicate neutron */
	  
		  part = DuplicateParticle(ptr, id);

		  /* Add to weight and population size */

		  wgt = wgt + RDB[ptr + PARTICLE_WGT];
		  np = np + 1;

		  /* Put particle in source*/

		  ToBank(part, id);
		}

	      /* Next particle */

	      ptr = PrevItem(ptr);

	    }

	}
    }
  else if (nsrc > Nlive)
    {
      /* Calculate probability */

      P = 1.0 - ((double)Nlive)/((double)nsrc);

      /* Loop over threads */
  
      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	{
	  /* Get particles from bank */

	  /* Get pointer to bank */

	  ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Get pointer to dummy */

	  ptr = FirstItem(ptr);
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Get pointer to first item after dummy */

	  ptr = NextItem(ptr);

	  while(ptr > VALID_PTR)
	    {
	      /* Check type */
	  
	      if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
		Die(FUNCTION_NAME, "Invalid particle type");

	      /* Get pointer and next */
	  
	      part = ptr;
	      ptr = NextItem(ptr);

	      /* Sample removal */

	      if (RandF(0) < P)
		{
		  /* Remove particle */

		  RemoveItem(part);

		  /* subtract from weight and population size */
	      
		  wgt = wgt - RDB[part + PARTICLE_WGT];
		  np = np - 1;

		  /* Put particle back to stack */

		  ToStack(part, id);

		}
	    }
	}
    }

  /* Check */

  if (np < 1)
    {
      /* All initial neutrons on next interval will be emitted */

      /* Calculate number of live neutrons and neutrons to emit */

      Nlive = 0;
      Nemit = nbatch;

#ifdef DNPRINT
      fprintf(out, "Nlive %ld\n", Nlive);
      fprintf(out, "Nemit %ld\n", Nemit);
#endif

      /* Store numbers */

      WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
      WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

    }
  else if (np > nbatch)
    {
      /* All initial neutrons on next interval will be live */

      /* Calculate number of live neutrons and neutrons to emit */

      Nlive = np;
      Nemit = 0;

#ifdef DNPRINT
      fprintf(out, "Nlive %ld\n", Nlive);
      fprintf(out, "Nemit %ld\n", Nemit);
#endif

      /* Store numbers */

      WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
      WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

    }
  else
    {

      /* Calculate eventual number of live neutrons and neutrons to emit */

      Nlive = np;
      Nemit = nbatch - np;

#ifdef DNPRINT
      fprintf(out, "Nlive %ld\n", Nlive);
      fprintf(out, "Nemit %ld\n", Nemit);
#endif

      /* Store numbers */

      WDB[loc0 + PRECDET_N_EMIT] = (double)Nemit;
      WDB[loc0 + PRECDET_N_LIVE] = (double)Nlive;

    }

  /***************************************************************************/

  /***** Normalize source ****************************************************/

  /* Normalize weights */

  /* Loop over threads */
  
  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      /* Get particles from bank */

      /* Get pointer to bank */

      ptr = (long)RDB[OMPPtr(DATA_PART_PTR_BANK, id)];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to dummy */

      ptr = FirstItem(ptr);
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      /* Get pointer to first item after dummy */

      while((ptr = NextItem(ptr)) > VALID_PTR)
	{
	  /* Check type */
	  
	  if ((long)RDB[ptr + PARTICLE_TYPE] != PARTICLE_TYPE_NEUTRON)
	    Die(FUNCTION_NAME, "Invalid particle type");

	  /* Normalize */

	  WDB[ptr + PARTICLE_WGT] = RDB[ptr + PARTICLE_WGT]*wgt0/wgt;

	}
    }

#ifdef DNPRINT
  fprintf(out, "<-- resizedynsrc.c\n\n");  
#endif

}

/*****************************************************************************/
