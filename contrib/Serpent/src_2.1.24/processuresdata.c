/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processuresdata.c                              */
/*                                                                           */
/* Created:       2011/01/08 (JLe)                                           */
/* Last modified: 2015/04/24 (TVi)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes unresolved resonance probability table data        */
/*                                                                           */
/* Comments: - Serpent 1 hakee ensimmäisen energiapisteen etukäteen          */
/*             generoidusta listasta, mikä voi nopeuttaa sämpläystä. Samaa   */
/*             voisi harkita tähän.                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessUresData:"

/*****************************************************************************/

void ProcessUresData(long nuc) 
{
  long L, N, M, INT, IFF, loc1, neg, n, m, i, j, rea, rea0, ptr, urs, nuc1;
  long pte, ptp, ace, JXS[32];
  double *E0, f, *XSS;

  /* Check ures flag */

  if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED))
    return;

  /* Check TMS */

  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS)
    {
      /* Reset ures flag */

      ResetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_USED);

      /* Exit */

      return;
    }

  /* Pointer to ACE data */

  ace = (long)RDB[nuc + NUCLIDE_PTR_ACE];
  
  /* Read data to JXS array */
  
  ptr = (long)ACE[ace + ACE_PTR_JXS];
  
  for (n = 0; n < 32; n++)
    JXS[n] = (long)ACE[ptr++];
  
  /* Get pointer to XSS array */
  
  XSS = &ACE[(long)ACE[ace + ACE_PTR_XSS]];

  /* Get pointer to UNR block */

  if ((L = JXS[22] - 1) < 0)
    Die(FUNCTION_NAME, "Ures flag set but ACE block not found");

  /* Reset number of negative points */
  
  neg = 0.0;

  /***************************************************************************/

  /***** Read pointers and flags *********************************************/

  /* Read number of incident energies */

  N = (long)XSS[L];

  /* Read number of probabilities */

  M = (long)XSS[L + 1];

  /* Read interpolation parameter (2 = lin-lin, 5 = log-log) */

  INT = (long)XSS[L + 2];

  /* Check interpolation scheme, (2 = lin-lin, 5 = log-log) */

  if ((INT != 2) && (INT != 5))
    Die(FUNCTION_NAME, "Unrecognized interpolation scheme in URES tables. INT = %ld (nuclide %s)", INT, GetText(nuc + NUCLIDE_PTR_NAME));
  
  /* XSS[L + 3] is ILF (Inelastic competition flag), XSS[L + 4] is IOA */
  /* (other absorption flag). Neither is used here. */
  
  /* Other factors flag (0 = cross sections, 1 = factors of smooth) */

  IFF = (long)XSS[L + 5];

  /* Get pointer to incident energies */

  E0 = &XSS[L + 6];


  if (1 != 2)
    {
      /* Check points and discard data if negative values are found. This  */
      /* is the case for some natural elements and the values seem to make */
      /* sense without the minus sign. (27.9.2009) */

      for (n = 0; n < N; n++)
	if (E0[n] < 0.0)
	  {
	    /* Print warning */
	    
	    Warn(FUNCTION_NAME, 
		 "Overlapping ures regions, data discarded (%s)", 
		 GetText(nuc + NUCLIDE_PTR_NAME));

	    /* Reset ures options */

	    ResetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_AVAIL);
	    ResetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_URES_USED);

	    /* Update counters */

	    WDB[DATA_URES_AVAIL] = RDB[DATA_URES_AVAIL] - 1.0;
	    WDB[DATA_URES_USED] = RDB[DATA_URES_USED] - 1.0;
	    
	    /* Exit subroutine */
	    
	    return;
	  } 
    }
  else
    {
      /* Or take absolute values (28.9.2009). PURR sets negative energies */
      /* for materials consisting of multiple isotopes with overlapping   */
      /* ures regions. */
      
      for (n = 0; n < N; n++)
	E0[n] = fabs(E0[n]);
    }

  /* Set energy limits for nuclide */

  WDB[nuc + NUCLIDE_URES_EMIN] = E0[0];
  WDB[nuc + NUCLIDE_URES_EMAX] = E0[N - 1];

  /* Set interpolation scheme for nuclide */
	  
  WDB[nuc + NUCLIDE_URES_INT] = (double)INT;

  /* Check limits */

  CheckValue(FUNCTION_NAME, "Emin", "", E0[0], 1E-6, 20.0);
  CheckValue(FUNCTION_NAME, "Emax", "", E0[N - 1], E0[0], 20.0);

  /* Compare to global boundaries */
	      
  if (RDB[nuc + NUCLIDE_URES_EMIN] < RDB[DATA_URES_EMIN])
    WDB[DATA_URES_EMIN] = RDB[nuc + NUCLIDE_URES_EMIN];
  
  if (RDB[nuc + NUCLIDE_URES_EMAX] > RDB[DATA_URES_EMAX])
    WDB[DATA_URES_EMAX] = RDB[nuc + NUCLIDE_URES_EMAX];

  /* Allocate memory for random number (used for all reaction modes) */
  /* TVi 2015-01-23: Must be used also for all nuclides */

  /* Allocate NUCLIDE_PTR_URES_RND and set pointers to all nuclides sharing */
  /* the same ZAI */

  if ((long)RDB[nuc + NUCLIDE_PTR_URES_RND] < VALID_PTR)
    {
      /* Allocate memory */
      
      AllocValuePair(nuc + NUCLIDE_PTR_URES_RND);

      /* Loop over nuclides */

      nuc1 = (long)RDB[DATA_PTR_NUC0];
      while (nuc1 > VALID_PTR)
	{
	  /* Check ZAI and copy pointer */
	  
	  if ((long)RDB[nuc + NUCLIDE_ZAI] == (long)RDB[nuc1 + NUCLIDE_ZAI])
	    WDB[nuc1 + NUCLIDE_PTR_URES_RND] = RDB[nuc + NUCLIDE_PTR_URES_RND];
	  
	  /* Next */

	  nuc1 = NextItem(nuc1);
	}
    } 

  /* Update pointer */

  L = L + 6 + N;

  /***************************************************************************/
      
  /***** Store reaction data *************************************************/

  /* Store energy grid */

  pte = MakeEnergyGrid(N, 0, 0, -1, E0, EG_INTERP_MODE_LIN);

  /* Allocate memory for probabilities */

  ptp = ReallocMem(DATA_ARRAY, N*M);

  /* Store points */

  i = 0;
  for (n = 0; n < N; n++)
    {
      for (m = 0; m < M; m++)
	{
	  WDB[ptp + i++] = XSS[L + n*6*M + m];
	  
	  /* Check ascending order */

	  if (m > 0)
	    if (XSS[L + n*6*M + m] < XSS[L + n*6*M + m - 1])
	      Die(FUNCTION_NAME, "Probabilities not in ascending order");
	}
      
      /* Check last point */

      if (XSS[L + n*6*M + M - 1] != 1.0)
	Die(FUNCTION_NAME, "Probability distribution not normalized");
    }

  /* Loop over reaction modes */
      
  for (j = 0; j < 4; j++)
    {
      /* Get pointer to reaction data */
      
      rea = -1;
      
      if(j == 0)
	rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
      else if(j == 1)
	rea = (long)RDB[nuc + NUCLIDE_PTR_ELAXS];
      else if (j == 2)
	rea = (long)RDB[nuc + NUCLIDE_PTR_FISSXS];
      else if (j == 3)      
	rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];

      /* Check that channel exists */

      if (rea > VALID_PTR)
	{
	  /* Set minimum and maximum energy */

	  WDB[rea + REACTION_URES_EMIN] = RDB[nuc + NUCLIDE_URES_EMIN];
	  WDB[rea + REACTION_URES_EMAX] = RDB[nuc + NUCLIDE_URES_EMAX];

	  /* Set interpolation scheme for reaction */
	  
	  WDB[rea + REACTION_URES_INT] = RDB[nuc + NUCLIDE_URES_INT];

	  /* Allocate memory for ures data block */

	  urs = NewItem(rea + REACTION_PTR_URES, URES_BLOCK_SIZE);

	  /* Allocate memory for previous factor */

	  AllocValuePair(urs + URES_PTR_PREV_FACT);

	  /* Set pointer to grid */

	  WDB[urs + URES_PTR_EGRID] = (double)pte;

	  /* Set number of probabilities and pointer to grid */

	  WDB[urs + URES_NP] = (double)M;
	  WDB[urs + URES_PTR_PROB] = (double)ptp;

	  /* Link pointer to random number */

	  WDB[urs + URES_PTR_RND] = RDB[nuc + NUCLIDE_PTR_URES_RND];

	  /* Store interpolation flag */

	  WDB[urs + URES_IFF] = (double)IFF;

	  /* Allocate memory for factors */

	  ptr = ReallocMem(DATA_ARRAY, N*M);

	  /* Set pointer */
	  
	  WDB[urs + URES_PTR_FACT] = (double)ptr;

	  /* Store points */

	  i = 0;
	  for (n = 0; n < N; n++)
	    for (m = 0; m < M; m++)
	      {
		/* Get value */
		
		f = XSS[L + n*6*M + m + (j + 1)*M];
		
		/* Check */
		
		if ((f < -10.0) || (f > MAX_XS))
		  Die(FUNCTION_NAME, "Invalid table value %1.5E", f);
		else if (f < 0.0)
		  {
		    /* Add warning message counter */
		    
		    neg++;
		  }
		
		/* Add point */
		
		WDB[ptr + i++] = f;
	      }
	  
	  /* Allocate memory for maximum factors */
	      
	  loc1 = ReallocMem(DATA_ARRAY, N);
	  
	  /* Set pointer */
	  
	  WDB[urs + URES_PTR_MAXF] = (double)loc1;

	  /* Get maximum factors */
	  
	  for (n = 0; n < N; n++)
	    for (m = 0; m < M; m++)
	      {
		/* Value from first interpolation table */
		
		f = XSS[L + n*6*M + m + (j + 1)*M];
		
		/* Compare to maximum */

		if (f > RDB[loc1 + n])
		  WDB[loc1 + n] = f;
		
		/* Compare to maximum in next table */
		
		if ((n < N - 1) && (f > RDB[loc1 + n + 1]))
		  WDB[loc1 + n + 1] = f;

		/* Value from second interpolation table */
		
		if (n < N - 1)
		  {
		    f = XSS[L + (n + 1)*6*M + m + (j + 1)*M];
		
		    /* Compare to maximum */
		    
		    if (f > RDB[loc1 + n])
		      WDB[loc1 + n] = f;
		    
		    /* Compare to maximum in next table */
		    
		    if ((n < N - 1) && (f > RDB[loc1 + n + 1]))
		      WDB[loc1 + n + 1] = f;
		  }
	      }
	}  
    }

  /***************************************************************************/

  /***** Sum of absorptions **************************************************/

  /* Pointer to absorption cross section */

  rea = (long)RDB[nuc + NUCLIDE_PTR_SUM_ABSXS];
  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

  /* Pointer to (n,gamma) */

  rea0 = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);

  /* Copy data */
  
  WDB[rea + REACTION_PTR_URES] = WDB[rea0 + REACTION_PTR_URES];      
  WDB[rea + REACTION_URES_EMIN] = WDB[rea0 + REACTION_URES_EMIN];
  WDB[rea + REACTION_URES_EMAX] = WDB[rea0 + REACTION_URES_EMAX];
  WDB[rea + REACTION_PTR_URES_MAX] = WDB[rea0 + REACTION_PTR_URES_MAX];
  WDB[rea + REACTION_URES_MAX_N0] = WDB[rea0 + REACTION_URES_MAX_N0];
  WDB[rea + REACTION_URES_MAX_NP] = WDB[rea0 + REACTION_URES_MAX_NP];
  WDB[rea + REACTION_URES_INT] = WDB[rea0 + REACTION_URES_INT];

  /* Check for non-(n,gamma) capture reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check implicit capture, ty and mt */

      if (((long)RDB[DATA_OPT_IMPL_CAPT] == YES) &&
	  ((long)RDB[rea + REACTION_TY] == 0) && 
	  ((long)RDB[rea + REACTION_MT] != 102))
	Warn(FUNCTION_NAME, "Non-(n,gamma) reaction mode");
						 
      /* Next reaction */

      rea = NextItem(rea);
    }

  /***************************************************************************/

  /***** Copy data to channels with energy-dependent branching ***************/

  /* Check pointer to data */

  if ((long)RDB[nuc + NUCLIDE_BRA_TYPE] != BRA_TYPE_ENE)
    return;

  /* Get pointer to (n,gamma) to isomeric state */

  if ((rea = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS_ISO]) < VALID_PTR)
    return;

  /* Pointer to (n,gamma) to ground state */

  rea0 = (long)RDB[nuc + NUCLIDE_PTR_NGAMMAXS];
  CheckPointer(FUNCTION_NAME, "(rea0)", DATA_ARRAY, rea0);
  
  /* Copy data */
  
  WDB[rea + REACTION_PTR_URES] = WDB[rea0 + REACTION_PTR_URES];      
  WDB[rea + REACTION_URES_EMIN] = WDB[rea0 + REACTION_URES_EMIN];
  WDB[rea + REACTION_URES_EMAX] = WDB[rea0 + REACTION_URES_EMAX];
  WDB[rea + REACTION_PTR_URES_MAX] = WDB[rea0 + REACTION_PTR_URES_MAX];
  WDB[rea + REACTION_URES_MAX_N0] = WDB[rea0 + REACTION_URES_MAX_N0];
  WDB[rea + REACTION_URES_MAX_NP] = WDB[rea0 + REACTION_URES_MAX_NP];
  WDB[rea + REACTION_URES_INT] = WDB[rea0 + REACTION_URES_INT];
  
  /***************************************************************************/
}

/*****************************************************************************/



