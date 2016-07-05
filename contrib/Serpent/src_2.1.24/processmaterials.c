/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processmaterials.c                             */
/*                                                                           */
/* Created:       2010/12/28 (JLe)                                           */
/* Last modified: 2015/06/03 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes material compositions, sets nuclide and reaction   */
/*              lists.                                                       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessMaterials:"

/*****************************************************************************/

void ProcessMaterials()
{
  long mat, ptr, nmat, nbumat, i, TMS, iso, nuc;
  double T;

  /* Check burnup step */

  if ((long)RDB[DATA_BURN_STEP] > 0)
    Die(FUNCTION_NAME, "Should not be here");

  /* Pointer to material list */
  
  if ((mat = (long)RDB[DATA_PTR_M0]) < VALID_PTR)
    Error(0, "No material definitions");

  /***************************************************************************/

  /***** Finalize material compositions **************************************/
  
  fprintf(out, "Normalizing compositions and processing mixtures...\n");

  /* Calculate fractions */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      IsotopeFractions(mat);
      mat = NextItem(mat);
    }

  /* Process mixtures */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      ProcessMixture(mat, 0);
      mat = NextItem(mat);
    }

  fprintf(out, "OK.\n\n");

  /* Replace isotopic with atomic data in photon transport calculation */
  /* (NOTE: tätä kutsutaan myöhemmin myös tuolla alempana). */
  
  ReplacePhotonData();
  
  /* Identify TFB materials */

  IdentifyTFBMat();

  /***************************************************************************/
  
  /***** Sort composition for better memory management ***********************/

  /* Check if macroscopic cross sections are pregenerated (if */
  /* number of materials is large, the calculation hangs here) */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES)
    {
      /* Use MPI task numbers to remember initial order */
  
      i = 0;
      
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
	{
	  /* Check divided and burn flags */
	  
	  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
	    WDB[mat + MATERIAL_MPI_ID] = 1E+12;
	  else if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
	    WDB[mat + MATERIAL_MPI_ID] = 1E+11;
	  else
	    WDB[mat + MATERIAL_MPI_ID] = (double)(i++);
	  
	  /* Next material */
	  
	  mat = NextItem(mat);
	}
      
      /* Use OpenMP thread number for sort */
      
      i = 0;
      
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
	{
	  /* Check divided and burn flags */
	  
	  if ((long)RDB[mat + MATERIAL_DIV_TYPE] == MAT_DIV_TYPE_PARENT)
	    WDB[mat + MATERIAL_OMP_ID] = 1E+12;
	  else if (!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
	    WDB[mat + MATERIAL_OMP_ID] = 1E+11;
	  else
	    {
	      /* Set id */
	  
	      WDB[mat + MATERIAL_OMP_ID] = (double)(i++);
	      
	      /* Check id */
	      
	      if (i == (long)RDB[DATA_OMP_MAX_THREADS])
		i = 0;
	    }
	  
	  /* Next material */
	  
	  mat = NextItem(mat);
	}
      
      /* Sort */
      
      mat = (long)RDB[DATA_PTR_M0];
      SortList(mat, MATERIAL_OMP_ID, SORT_MODE_ASCEND);
    }

  /***************************************************************************/

  /***** Allocate memory and process *****************************************/

  /* Process compositions of burnable materials */

  BurnMatCompositions();
  
  /* Put composition to divided materials */

  PutCompositions();

  /* Override material compositions */
  
  ReadMaterialComp();

  /* Process burnable materials */

  ProcessBurnMat();

  /* Calculate masses (to get the value on output) */

  CalculateMasses();

  /* Close composition lists */
  
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    { 
      ptr = (long)RDB[mat + MATERIAL_PTR_COMP];
      CloseList(ptr);

      /* Tämä on lista alkuperäiseen nuklidikoostumukseen joka korvataan */
      /* alkuainekohtaisella koostumuksella fotonitransportmoodissa jos  */
      /* inputti on annettu neutronidatana (JLe 3.6.2015 / 2.1.24). */

      if ((ptr = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) > VALID_PTR)
	CloseList(ptr);

      mat = NextItem(mat);
    }

  /* Check temperatures and TMS flags */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Get temperature */
	  
      T = RDB[mat + MATERIAL_DOPPLER_TEMP];

      /* Get TMS flag */

      if (RDB[mat + MATERIAL_TMS_MODE] != TMS_MODE_NONE)
	TMS = NUCLIDE_FLAG_TMS;
      else
	TMS = 0;
      
      /* Loop over composition */
      
      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
      while (iso > VALID_PTR)
	{
	  /* Pointer to nuclide */

	  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	  /* Skip lost */

	  if (nuc == (long)RDB[DATA_PTR_NUCLIDE_LOST])
	    {
	      /* Next */

	      iso = NextItem(iso);

	      /* Cycle loop */

	      continue;
	    }
	  
	  /* Compare temperature */

	  if ((T > 0) && (T != RDB[nuc + NUCLIDE_TEMP]))
	    Die(FUNCTION_NAME, "Error in temperature: %s %s %E %E",
		GetText(mat + MATERIAL_PTR_NAME), 
		GetText(nuc + NUCLIDE_PTR_NAME), T, RDB[nuc + NUCLIDE_TEMP]);
	  
	  /* Check TMS flag */

	  if (TMS != ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_TMS))
	    Die(FUNCTION_NAME, "Error in TMS flag");

	  /* Check TMS limits */

	  if (TMS != NO)
	    {
	      /* Check minimum */

	      if (RDB[nuc + NUCLIDE_TMS_MIN_TEMP] > 
		  RDB[mat + MATERIAL_TMS_TMIN])
		Die(FUNCTION_NAME, "Error in TMS Tmin: %s %s %E %E",
		    GetText(mat + MATERIAL_PTR_NAME), 
		    GetText(nuc + NUCLIDE_PTR_NAME),
		    RDB[nuc + NUCLIDE_TMS_MIN_TEMP], 
		    RDB[mat + MATERIAL_TMS_TMIN]);

	      /* Check maximum */

	      if (RDB[nuc + NUCLIDE_TMS_MAX_TEMP] < 
		  RDB[mat + MATERIAL_TMS_TMAX])
		Die(FUNCTION_NAME, "Error in TMS Tmax: %s %s %E %E",
		    GetText(mat + MATERIAL_PTR_NAME), 
		    GetText(nuc + NUCLIDE_PTR_NAME),
		    RDB[nuc + NUCLIDE_TMS_MAX_TEMP], 
		    RDB[mat + MATERIAL_TMS_TMAX]);
	    }
	  
	  /* Next */

	  iso = NextItem(iso);
	}

      /* Next material */

      mat = NextItem(mat);
    }

  /* Allocate memory for reaction lists and macroscopic cross sections */

  AllocMacroXS();

  /* Sort composition to get initial order */

  if ((long)RDB[DATA_OPTI_RECONSTRUCT_MACROXS] == YES)
    {
      mat = (long)RDB[DATA_PTR_M0];
      SortList(mat, MATERIAL_MPI_ID, SORT_MODE_ASCEND);
    }

  /* Re-read compositions from restart file (JLe: tämä pitää lukea uudestaan  */
  /* siitä syystä että alialuejako tehdään vasta PutCompositions():issa, jota */
  /* kutsutaan tuolla ylempänä). */

  ReadRestartFile(RESTART_OVERRIDE);

  /* Replace isotopic with atomic data in photon transport calculation */
  /* (JLe: Tätä joudutaan kutsumaan uudestaan että myös alialueiden    */
  /* koostumukset menee oikein). */

  ReplacePhotonData();

  /* Check if decay source is used */

  if ((long)RDB[DATA_USE_DECAY_SRC] == YES)
    {
      /* Calculate activities (Aktiivisuudet lasketaan uudestaan */
      /* transportcycle.c:n lopussa.) */

      CalculateActivities();

      /* Process decay source */
      
      ProcessDecaySrc();
      
      /* Print gamma source (spectra for testing) */
      
      PrintGammaSpectra();
    }

  /* Process photon attenuation data if photon transport mode */

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    ProcessPhotonAtt();

  /***************************************************************************/
  
  /***** Assign MPI numbers to materials *************************************/

  /* NOTE: tää pitää miettiä uudestaan niin että toi järjestys menee */
  /* jotenkin fiksusti */

  /* Set MPI id's */

  i = 0;

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check burn-flag */

      if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
	{
	  /* Set id */

	  WDB[mat + MATERIAL_MPI_ID] = (double)(i++);
	  
	  /* Check id */
	  
	  if (i == mpitasks)
	    i = 0;
	}
      else
	WDB[mat + MATERIAL_MPI_ID] = -1.0;
            
      /* Reset OpenMP thread number */

      WDB[mat + MATERIAL_OMP_ID] = -1.0;

      /* Next material */

      mat = NextItem(mat);
    }
  
  /***************************************************************************/
  
  /***** Count materials and set indexes for printout ************************/

  /* Reset counters */

  nmat = 0;
  nbumat = 0;

  /* Loop over materials */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Add total counter */

      nmat++;

      /* Put index */

      WDB[mat + MATERIAL_PROC_IDX] = (double)nmat;

      /* Check burn flag */
      
      if (((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT) &&
	  ((long)RDB[mat + MATERIAL_DIV_TYPE] != MAT_DIV_TYPE_PARENT))
	{
	  /* Add burn counter */

	  nbumat++;

	  /* Put index */

	  WDB[mat + MATERIAL_BURN_IDX] = (double)nmat;
	}
      
      /* Next material */
      
      mat = NextItem(mat);
    }

  /* Put counters */

  WDB[DATA_N_MATERIALS] = (double)nmat;
  WDB[DATA_N_BURN_MATERIALS] = (double)nbumat;

  /***************************************************************************/

  /***** Calculate memory size and print summary *****************************/

  /* Calculated divided material total memory */

  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check pointer to parent */

      if ((ptr = (long)RDB[mat + MATERIAL_DIV_PTR_PARENT]) > VALID_PTR)
	WDB[ptr + MATERIAL_TOT_DIV_MEM_SIZE] =
	  RDB[ptr + MATERIAL_TOT_DIV_MEM_SIZE] + RDB[mat + MATERIAL_MEM_SIZE];

      /* Next material */

      mat = NextItem(mat);
    }

  /* Print material data */
  
  PrintMaterialData();

  /***************************************************************************/
}

/*****************************************************************************/
