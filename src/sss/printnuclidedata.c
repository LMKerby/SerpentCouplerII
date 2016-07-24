#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printnuclidedata.c                             */
/*                                                                           */
/* Created:       2011/01/24 (JLe)                                           */
/* Last modified: 2014/02/21 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Prints misc. nuclide data in standard output                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintNuclideData:"

/*****************************************************************************/

void PrintNuclideData(long nuc, long all)
{
  long ZAI, ace, ptr;
  double t, mem; 
  static char tmpstr[MAX_STR];
 
  /* Check nuclide pointer */

  if (nuc > VALID_PTR)
    {
      /***** Print nuclide-wise data *****************************************/
      
      /* Check if nuclide has transport data */

      if (((long)RDB[DATA_BURN_DECAY_CALC] == NO) &&
	  !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TRANSPORT_DATA)
	  &&
	  !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DOSIMETRY_DATA) 
	  && 
	  !((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_PHOTON_DATA))
	return;

      /* Get ZAI */
      
      ZAI = RDB[nuc + NUCLIDE_ZAI];
      
      /* Use temporary string to print name */
      
      sprintf(tmpstr, "%s", ZAItoIso(ZAI,2)); 

      if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
	fprintf(out, "Element %10s -- %s (photon library)", 
		GetText(nuc + NUCLIDE_PTR_NAME), 
		ZAItoIso((long)RDB[nuc + NUCLIDE_Z],1));
      else if (RDB[nuc + NUCLIDE_TEMP] > 0.0)
	fprintf(out, "Nuclide %10s -- %s at %1.0fK (%s)", 
		GetText(nuc + NUCLIDE_PTR_NAME), tmpstr,
		RDB[nuc + NUCLIDE_TEMP], ZAItoIso(ZAI,1));
      else
	fprintf(out, "Nuclide %10s -- %s (%s)", 
		GetText(nuc + NUCLIDE_PTR_NAME), tmpstr,
		ZAItoIso(ZAI,1));

      /* Check mode */

      if (all == NO)
	{
	  fprintf(out, "\n");

	  return;
	}
      else if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_PHOTON)
	{
	  if (RDB[DATA_N_TRANSPORT_NUCLIDES] > 0.0)
	    fprintf(out, "\n\n");
	  else
	    fprintf(out, "\n");

	  return;
	}
      else
	fprintf(out, ":\n");

      /* Print number of reactions */
      
      if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] == 0)
	fprintf(out, " - %ld reaction channels\n", 
		(long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA]);
      else if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] == 1)
	fprintf(out, 
		" - %ld reaction channels + 1 branch\n", 
		(long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA]);
      else
	fprintf(out, 
		" - %ld reaction channels + %ld branches\n", 
		(long)RDB[nuc + NUCLIDE_N_TRANSPORT_REA],
		(long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH]);
      
      if ((long)RDB[nuc + NUCLIDE_N_SPECIAL_REA] > 0)
	fprintf(out, " - %ld special cross sections\n", 
		(long)RDB[nuc + NUCLIDE_N_SPECIAL_REA]);
      
      if ((long)RDB[DATA_BURNUP_CALCULATION_MODE] == YES)
	{
	  if ((long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA] == 0)
	    fprintf(out, " - No transmutation reactions\n");
	  else if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] == 0)
	    fprintf(out, " - %ld transmutation reactions\n", 
		    (long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA]);
	  else if ((long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH] == 1)
	    fprintf(out, " - %ld transmutation reactions + 1 branch\n", 
		    (long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA]);
	  else
	    fprintf(out, " - %ld transmutation reactions + %ld branches\n", 
		    (long)RDB[nuc + NUCLIDE_N_TRANSMUTATION_REA],
		    (long)RDB[nuc + NUCLIDE_N_TRANSPORT_BRANCH]);
	  
	  /* Check if decay data is read */
	  
	  if ((long)RDB[DATA_PTR_DECDATA_FNAME_LIST] > 0)
	    {	  
	      /* Check decay constant */
	      
	      if ((ace = (long)RDB[nuc + NUCLIDE_PTR_DECAY_ACE]) 
		  < VALID_PTR)
		fprintf(out, " - Nuclide has no decay data\n");
	      else if (ACE[ace + ACE_LAMBDA] == 0.0)
		fprintf(out, " - Nuclide is stable\n");
	      else if (RDB[nuc + NUCLIDE_LAMBDA] == 0.0)
		fprintf(out, 
			" - Nuclide is stable beyond cutoff (Half-life > %s)\n",
			TimeIntervalStr(RDB[DATA_DEP_HALF_LIFE_CUTOFF]));
	      else
		{
		  if ((long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH] == 0)
		    fprintf(out, " - %ld decay modes\n", 
			    (long)RDB[nuc + NUCLIDE_N_DECAY_REA]);
		  else if ((long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH] == 1)
		    fprintf(out, 
			    " - %ld decay modes + 1 branch\n", 
			    (long)RDB[nuc + NUCLIDE_N_DECAY_REA]);
		  else
		    fprintf(out, 
			    " - %ld reaction decay modes + %ld branches\n", 
			    (long)RDB[nuc + NUCLIDE_N_DECAY_REA],
			    (long)RDB[nuc + NUCLIDE_N_DECAY_BRANCH]);
		  
		  /* Convert to half-life */
		  
		  t = log(2.0)/RDB[nuc + NUCLIDE_LAMBDA];
		  
		  /* Print */
		  
		  if (t < 60.0)
		    fprintf(out, " - Half-life = %s \n", 
			    TimeIntervalStr(t));
		  else
		    fprintf(out, " - Half-life = %1.2E seconds (%s)\n", t, 
			    TimeIntervalStr(t));
		}
	    }
	  
	  /* Check if NFY data is used */
	  
	  if ((long)RDB[DATA_PTR_NFYDATA_FNAME_LIST] > 0)
	    {
	      /* Pointer to data */
	      
	      if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_NFY_DATA]) 
		  < VALID_PTR)
		fprintf(out, " - Nuclide has no NFY data\n");
	      else
		{
		  if ((long)RDB[ptr + FISSION_YIELD_PARENT_ZAI] !=
		      (long)RDB[nuc + NUCLIDE_ZAI])
		    fprintf(out, 
			    " - Using NFY data from %s\n",
			    ZAItoIso((long)RDB[ptr + 
						FISSION_YIELD_PARENT_ZAI], 1));
		  else
		    fprintf(out, " - NFY data included\n");
		}
	    }
	  
	  /* Check if SFY data is used */
	  
	  if ((long)RDB[DATA_PTR_SFYDATA_FNAME_LIST] > 0)
	    {
	      /* Pointer to data */
	      
	      if ((ptr = (long)RDB[nuc + NUCLIDE_PTR_SFY_DATA]) 
		  < VALID_PTR)
		fprintf(out, " - Nuclide has no SFY data\n");
	      else
		fprintf(out, " - SFY data included\n");
	    }
	}
      
      if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_AVAIL))
	fprintf(out, " - No ures ptable data available\n"); 
      else if (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] &
		 NUCLIDE_FLAG_URES_USED))
	fprintf(out, " - Ures ptable data available but not used\n"); 
      else
	fprintf(out, " - Ures ptable data available and used\n"); 
      
      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_SAB_DATA)
	fprintf(out, " - Nuclide has S(a,b) data\n"); 
      else
	fprintf(out, " - No S(a,b) data\n"); 
      
      if ((mem = RDB[nuc + NUCLIDE_MEMSIZE]) > 0.0)
	{
	  if ((mem  = mem/KILO)< KILO)
	    fprintf(out, " - %1.2f kb of memory allocated for data\n", mem);
	  else if ((mem  = mem/KILO) < KILO)
	    fprintf(out, " - %1.2f Mb of memory allocated for data\n", mem);
	  else
	    fprintf(out, " - %1.2f Gb of memory allocated for data\n", 
		    mem/KILO);
	}
      
      fprintf(out, "\n");

      /***********************************************************************/
    }
  else
    {
      /***** Print summary ***************************************************/
      
      fprintf(out, "\n");

      fprintf(out, "SUMMARY -- %ld nuclides included in calculation:\n\n",
	      (long)RDB[DATA_N_TRANSPORT_NUCLIDES] + 
	      (long)RDB[DATA_N_DOSIMETRY_NUCLIDES] +
	      (long)RDB[DATA_N_DECAY_NUCLIDES] +
	      (long)RDB[DATA_N_PHOTON_NUCLIDES]);

      if (RDB[DATA_N_TRANSPORT_NUCLIDES] > 0.0)
	fprintf(out, " - %ld transport nuclides\n", 
		(long)RDB[DATA_N_TRANSPORT_NUCLIDES]);
      
      if (RDB[DATA_N_DOSIMETRY_NUCLIDES] > 0.0)
	fprintf(out, " - %ld dosimetry nuclides\n", 
		(long)RDB[DATA_N_DOSIMETRY_NUCLIDES]);
      
      if (RDB[DATA_N_DECAY_NUCLIDES] > 0.0)
	fprintf(out, " - %ld decay nuclides (not listed above)\n", 
		(long)RDB[DATA_N_DECAY_NUCLIDES]);
      
      if (RDB[DATA_N_PHOTON_NUCLIDES] > 0.0)
	fprintf(out, " - %ld elements with photon interaction data\n", 
		(long)RDB[DATA_N_PHOTON_NUCLIDES]);

      if ((long)RDB[DATA_N_TRANSPORT_NUCLIDES] > 0)
	fprintf(out, " - Neutron energy cut-offs at %1.2E and %1.1f MeV\n", 
		RDB[DATA_NEUTRON_EMIN], RDB[DATA_NEUTRON_EMAX]);

      if (RDB[DATA_N_PHOTON_NUCLIDES] > 0.0)
	fprintf(out, " - Photon energy cut-offs at %1.2E and %1.1f MeV\n", 
		RDB[DATA_PHOTON_EMIN], RDB[DATA_PHOTON_EMAX]);

      if (RDB[DATA_N_TRANSPORT_REA] > 0.0)
	fprintf(out, " - %ld transport reactions\n", 
		(long)RDB[DATA_N_TRANSPORT_REA]);

      if (RDB[DATA_N_SPECIAL_REA] > 0.0)
	fprintf(out, " - %ld special reactions\n", 
		(long)RDB[DATA_N_SPECIAL_REA]);
     
      if (RDB[DATA_N_DECAY_REA] > 0.0)
	{      
	  if ((long)RDB[DATA_N_DECAY_BRANCH] == 0)
	    fprintf(out, " - %ld decay reactions\n", 
		    (long)RDB[DATA_N_DECAY_REA]);
	  else if ((long)RDB[DATA_N_DECAY_BRANCH] == 1)
	    fprintf(out, 
		    " - %ld decay reactions + 1 branch\n", 
		    (long)RDB[DATA_N_DECAY_REA]);
	  else
	    fprintf(out, 
		    " - %ld decay reactions + %ld branches\n", 
		    (long)RDB[DATA_N_DECAY_REA],
		    (long)RDB[DATA_N_DECAY_BRANCH]);
	}

      if (RDB[DATA_N_TRANSMUTATION_REA] > 0.0)
	{      
	  if ((long)RDB[DATA_N_TRANSPORT_BRANCH] == 0)
	    fprintf(out, " - %ld transmutation reactions\n", 
		    (long)RDB[DATA_N_TRANSMUTATION_REA]);
	  else if ((long)RDB[DATA_N_TRANSPORT_BRANCH] == 1)
	    fprintf(out, 
		    " - %ld transmutation reactions + 1 branch\n", 
		    (long)RDB[DATA_N_TRANSMUTATION_REA]);
	  else
	    fprintf(out, 
		    " - %ld transmutation reactions + %ld branches\n", 
		    (long)RDB[DATA_N_TRANSMUTATION_REA],
		    (long)RDB[DATA_N_TRANSPORT_BRANCH]);
	}

      /* Calculate total memory */

      mem = 0;
      nuc = (long)RDB[DATA_PTR_NUC0];
      
      while (nuc > VALID_PTR)
	{
	  /* Add to memory */
	  
	  mem = mem + RDB[nuc + NUCLIDE_MEMSIZE];

	  /* Next nuclide */

	  nuc = NextItem(nuc);
	}

      /* Print */

      if ((mem  = mem/KILO) < KILO)
	fprintf(out, " - %1.2f kb of memory allocated for data\n", mem);
      else if ((mem  = mem/KILO) < KILO)
	fprintf(out, " - %1.2f Mb of memory allocated for data\n", mem);
      else
	fprintf(out, " - %1.2f Gb of memory allocated for data\n", 
		mem/KILO);

      /* List ures nuclides */
      
      if ((long)RDB[DATA_URES_USED] > 0)
	{
	  if ((long)RDB[DATA_URES_USED] == 1)
	    fprintf(out, " - One nuclide with ures ptable data:\n\n");
	  else
	    fprintf(out, " - %ld nuclides with ures ptable table data:\n\n",
		    (long)RDB[DATA_URES_USED]);
	  
	  nuc = (long)RDB[DATA_PTR_NUC0];
	  while (nuc > VALID_PTR)
	    {
	      if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_URES_USED)
		fprintf(out, "%12s (from %1.2E to %1.2E MeV)\n", 
			GetText(nuc + NUCLIDE_PTR_NAME),
			RDB[nuc + NUCLIDE_URES_EMIN], 
			RDB[nuc + NUCLIDE_URES_EMAX]);
	      
	      /* Next nuclide */
	      
	      nuc = NextItem(nuc);
	    }
	}

      fprintf(out, "\n");

      /***********************************************************************/
    }
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
