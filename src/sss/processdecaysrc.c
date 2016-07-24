#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdecaysrc.c                              */
/*                                                                           */
/* Created:       2015/05/23 (JLe)                                           */
/* Last modified: 2015/06/13 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes decay source into a sampling list                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessDecaySrc:"

/*****************************************************************************/

void ProcessDecaySrc()
{  
  long mat, iso, nuc, ptr, src;
  double vol, adens, lambda, tot, prev, I;

  fprintf(out, "Processing decay source...\n");
  
  /* Check total source rate */
  
  if (RDB[DATA_TOT_PHOTON_SRC_RATE] == 0.0)
    Error(0, "Total photon source rate is zero in decay source mode");

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get volume */

      vol = RDB[mat + MATERIAL_VOLUME];

      /* Check total photon source rate */
	  
      if ((vol == 0.0) || ((RDB[mat + MATERIAL_PHOTON_SRC_RATE]/
			    RDB[DATA_TOT_PHOTON_SRC_RATE]) < 1E-19))
	{
	  /* Next material */
	  
	  mat = NextItem(mat);

	  /* Cycle loop */

	  continue;
	}

      /* Check that pointer is not defined (NOTE: Tää on lähinnä sitä   */
      /* varten että jos tätä rutiinia joskus kutsutaan silmukasta niin */
      /* muistia ei varata turhaan. */

      if ((long)RDB[mat + MATERIAL_PTR_DECAY_SRC] > VALID_PTR)
	Die(FUNCTION_NAME, "Pointer to decay source already exists");

      /* Avoid compiler warning */

      src = -1;

      /* Reset total */

      tot = 0.0;

      /* Loop over composition */

      iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP];
      while (iso > VALID_PTR)
	{
	  /* Get atomic density */
	  
	  adens = RDB[iso + COMPOSITION_ADENS]*1E+24;

	  /* Get pointer to nuclide data */

	  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	  /* Get decay constant */

	  lambda = RDB[nuc + NUCLIDE_LAMBDA];

	  /* Get total intensity */

	  I = RDB[nuc + NUCLIDE_SPEC_PHOTON_I];

	  /* Check intensity */

	  if (I*lambda*adens*vol/RDB[mat + MATERIAL_PHOTON_SRC_RATE] < 1E-18)
	    {
	      /* Skip nuclide */

	      iso = NextItem(iso);

	      /* Cycle loop */

	      continue;
	    }

	  /* Create entry */

	  src = NewItem(mat + MATERIAL_PTR_DECAY_SRC, SRC_DECCAY_BLOCK_SIZE);
	      
	  /* Put nuclide pointer */

	  WDB[src + SRC_DECCAY_PTR_NUCLIDE] = (double)nuc;
	      
	  /* Put emission rate */

	  WDB[src + SRC_DECCAY_I] = I*lambda*adens*vol;

	  /* Reset weight */

	  WDB[src + SRC_DECCAY_WGT] = 1.0;

	  /* Add to total */

	  tot = tot + RDB[src + SRC_DECCAY_I];

	  /* Put pointer to line spectra */

	  ptr = (long)RDB[nuc + NUCLIDE_PTR_PHOTON_LINE_SPEC];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  WDB[src + SRC_DECCAY_PTR_SPEC] = (double)ptr;
	  
	  /* Next nuclide */

	  iso = NextItem(iso);
	}
     
      /* Check if source was created */

      if (src > VALID_PTR)
	{
	  /* Close list */

	  CloseList(src);

	  /* Sort */

	  SortList(src, SRC_DECCAY_I, SORT_MODE_DESCEND);

	  /* Check total */

	  if (tot == 0.0)
	    Die(FUNCTION_NAME, "WTF?");
	  
	  /* Calculate cumulative probabilities */

	  prev = 0.0;
	  
	  src = (long)RDB[mat + MATERIAL_PTR_DECAY_SRC];
	  while (src > VALID_PTR)
	    {
	      /* Calculate cumulative probability */
	      
	      WDB[src + SRC_DECCAY_CUM_P] = 
		(prev + RDB[src + SRC_DECCAY_I])/tot;
	      
	      /* Update previous */
	      
	      prev = prev + RDB[src + SRC_DECCAY_I];

	      /* Next */
	      
	      src = NextItem(src);
	    }

	  /* Allocate memory for sampled stats */

	  ptr = NewStat("SAMPLED_PHOTON_SRC", 1, 1);
	  WDB[mat + MATERIAL_SAMPLED_PHOTON_SRC] = (double)ptr;  
	}

      /* Next material */

      mat = NextItem(mat);
    } 

  fprintf(out, "OK.\n\n");
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
