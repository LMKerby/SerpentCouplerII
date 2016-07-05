/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readrestartfile.c                              */
/*                                                                           */
/* Created:       2014/01/23 (JLe)                                           */
/* Last modified: 2015/10/01 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Reads material compositions from a reastart file before      */
/*              running the simulation.                                      */
/*                                                                           */
/* Comments: RESTART_CHECK    -- Check compatibility with materials          */
/*           RESTART_OVERRIDE -- Override material compositions (used in     */
/*                               burnup mode)                                */
/*           RESTART_REPLACE  -- Replace initial compositions (used in       */
/*                               transport mode)                             */
/*                                                                           */
/*           - Replace mode works only with photon data (neutron transport   */
/*             calculation requires full nuclide names to be stored and a    */
/*             revised file format)                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadRestartFile:"

/*****************************************************************************/

void ReadRestartFile(long mode)
{
  long sz, n, nnuc, zai, mat, iso, nuc, *ptr, ok, idx, i;
  double pt, bu, days, adens, mdens, mbu, d0, closest;
  char tmpstr[MAX_STR], fname[MAX_STR], ZAI[MAX_STR];
  FILE *fp;

  /* Check if file is read */

  if ((long)RDB[DATA_READ_RESTART_FILE] == NO)
    return;

  if (mode == RESTART_CHECK)
    fprintf(out, "Checking restart file...\n");
  else
    fprintf(out, "Reading material compositions from restart file...\n");

  /* File name */

  if ((long)RDB[DATA_RESTART_READ_PTR_FNAME] > VALID_PTR)
    sprintf(fname, "%s", GetText(DATA_RESTART_READ_PTR_FNAME));
  else
    sprintf(fname, "%s.wrk", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for reading */

  if ((fp = fopen(fname, "r")) == NULL)
    Error(0, "Restart file \"%s\" does not exist", fname);

  /***************************************************************************/

  /***** Find point **********************************************************/
	  
  /* Check if index is given */

  if ((idx = (long)RDB[DATA_RESTART_READ_IDX]) > 0)
    {
      /* Reset point */

      pt = -1;
    }
  else
    {
      /* Get point */
      
      if ((pt = RDB[DATA_RESTART_READ_POINT]) == 0.0)
	Error(0, "Restart to zero burnup");
    }

  /* Reset closest or point if index is given */

  closest = INFTY;
  d0 = 0.0;
  i = 0;
     
  /* Read loop */
 
  while ((sz = fread(&n, sizeof(long), 1, fp)) > 0)
    {
      /* Read name */

      if ((sz = fread(tmpstr, sizeof(char), n, fp)) == 0)
	Error(0, "Error in restart file");
      
      /* Read nominal burnup and time */

      if ((sz = fread(&bu, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");
      
      if ((sz = fread(&days, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      /* Read number of nuclides, atomic and mass density and burnup */

      if ((sz = fread(&nnuc, sizeof(long), 1, fp)) == 0)
	Error(0, "Error in restart file");

      if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      if ((sz = fread(&mdens, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      if ((sz = fread(&mbu, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      /* Check values */

      CheckValue(FUNCTION_NAME, "nnuc", "", nnuc, 1, 10000);
      CheckValue(FUNCTION_NAME, "adens", "", adens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mdens", "", mdens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mbu", "", mbu, 0.0, 1000.0);
      
      /* Check index was given */

      if ((idx > 0) && (idx == i))
	{
	  /* Set point */

	  pt = -days;
	  closest = days;

	  /* Break loop */

	  break;
	}

      /* Check closest */

      if (pt > 0.0)
	{
	  if (fabs(bu - pt) < fabs(closest - pt))
	    closest = bu;
	}
      else if (pt < 0.0)
	{
	  if (fabs(days + pt) < fabs(closest + pt))
	    closest = days;
	}

      /* Skip nuclide data */
      
      n = (2*nnuc)*sizeof(double);
      fseek(fp, n, SEEK_CUR);

      /* Update counter */

      if (d0 != days)
	{
	  i++;
	  d0 = days;
	}
    }

  /* Rewind file */

  rewind(fp);

  /* Check point */
  
  if (((pt > 0.0) && (fabs(closest/pt - 1.0) > 0.001)) ||
      ((pt < 0.0) && (fabs(-closest/pt - 1.0) > 0.001)) ||
      ((idx > 0) && (idx != i)))
    {
      /* Point not found, print points */

      fprintf(stdout, 
	      "\nRestart file \"%s\" contains the following burnup points:\n\n",
	      fname);

      /* Reset days */

      d0 = -1;
      i = 0;

      /* Read loop */

      while ((sz = fread(&n, sizeof(long), 1, fp)) > 0)
	{
	  /* Read data */

	  if ((sz = fread(tmpstr, sizeof(char), n, fp)) == 0)
	    Die(FUNCTION_NAME, "Error in restart file");
	  if ((sz = fread(&bu, sizeof(double), 1, fp)) == 0)
	    Die(FUNCTION_NAME, "Error in restart file");
	  if ((sz = fread(&days, sizeof(double), 1, fp)) == 0)
	    Die(FUNCTION_NAME, "Error in restart file");
	  if ((sz = fread(&nnuc, sizeof(long), 1, fp)) == 0)
	    Die(FUNCTION_NAME, "Error in restart file");
	  if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
	    Die(FUNCTION_NAME, "Error in restart file");
	  if ((sz = fread(&mdens, sizeof(double), 1, fp)) == 0)
	    Die(FUNCTION_NAME, "Error in restart file");
	  if ((sz = fread(&mbu, sizeof(double), 1, fp)) == 0)
	    Die(FUNCTION_NAME, "Error in restart file");

	  /* Print */

	  if (days != d0)
	    {
	      /* Print burnup and time */

	      fprintf(stdout, "%2ld : BU = %7.3f MWd/kgu, time = %1.5E days", 
		      i++, bu, days);

	      /* Identify closest point */
	      
	      if ((idx == 0) && (((pt > 0) && (bu == closest)) || 
				 ((pt < 0) && (days == closest))))
		fprintf(stdout, " (closest point)\n");
	      else
		fprintf(stdout, "\n");
	    }

	  /* Set days */

	  d0 = days;

	  /* Skip material data */
	  
	  n = (2*nnuc)*sizeof(double);
	  fseek(fp, n, SEEK_CUR);
	}

      if (idx > 0)
	Error(0, "Burnup point %ld not found in restart file", idx);
      else if (pt > 0.0)
	Error(0, "Burnup point %1.3f MWd/kgU not found in restart file", pt);
      else
	Error(0, "Burnup point %1.3f days not found in restart file", -pt);
    }
  else if ((pt > 0.0) && (fabs(closest/pt - 1.0) > 1E-6))
    Note(0, "No exact match in restart file, using %1.2f MWd/kgU", closest);
  else if ((pt < 0.0) && (fabs(-closest/pt - 1.0) > 1E-6))
    Note(0, "No exact match in restart file, using %s", 
	 TimeIntervalStr(closest*60*60*24));

  /* Put point to closest value */

  if (pt > 0)
    pt = closest;
  else
    pt = -closest;

  /***************************************************************************/

  /***** Read data ***********************************************************/
  
  /* Reset ok flag */

  ok = NO;
    
  /* Read loop */

  while ((sz = fread(&n, sizeof(long), 1, fp)) > 0)
    {
      /* Read name */

      if ((sz = fread(tmpstr, sizeof(char), n, fp)) == 0)
	Error(0, "Error in restart file");
      
      /* Put string terminator */

      tmpstr[n] = '\0';

      /* Read nominal burnup and time */

      if ((sz = fread(&bu, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");
      
      if ((sz = fread(&days, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      /* Read number of nuclides, atomic and mass density and burnup */

      if ((sz = fread(&nnuc, sizeof(long), 1, fp)) == 0)
	Error(0, "Error in restart file");

      if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      if ((sz = fread(&mdens, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      if ((sz = fread(&mbu, sizeof(double), 1, fp)) == 0)
	Error(0, "Error in restart file");

      /* Check values */

      CheckValue(FUNCTION_NAME, "nnuc", "", nnuc, 1, 10000);
      CheckValue(FUNCTION_NAME, "adens", "", adens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mdens", "", mdens, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "mbu", "", mbu, 0.0, 1000.0);

      /* Check point */

      if (((pt > 0.0) && (fabs(bu/pt - 1.0) > 1E-12)) ||
	  ((pt < 0.0) && (fabs(-days/pt - 1.0) > 1E-12)))
	{
	  /* Skip material */
      
	  n = (2*nnuc)*sizeof(double);
	  fseek(fp, n, SEEK_CUR);
	}
      else
	{
	  /* Set flag */

	  ok = YES;

	  /* Check check mode */

	  if (mode == RESTART_CHECK)
	    {
	      /* Skip material */
      
	      n = (2*nnuc)*sizeof(double);
	      fseek(fp, n, SEEK_CUR);

	      /* Cycle loop */
	      
	      continue;
	    }

	  /* Set burnup and irradiation time */

	  WDB[DATA_BURN_CUM_BURNUP] = bu;
	  WDB[DATA_BURN_CUM_BURNTIME] = days*24.0*60.0*60.0;

	  /* Find material */

	  mat = (long)RDB[DATA_PTR_M0];
	  while (mat > VALID_PTR)
	    {
	      /* Compare */

	      if (!strcmp(tmpstr, GetText(mat + MATERIAL_PTR_NAME)))
		break;
	      
	      /* Next material */

	      mat = NextItem(mat);
	    }
	
	  /* Check mode */

	  if (mode == RESTART_REPLACE)
	    {
	      /***************************************************************/

	      /***** Replace entire material data ****************************/

	      /* Check material pointer and create new if not found */

	      if (mat < VALID_PTR)
		mat = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);

	      /* Put parameter name and file name */

	      WDB[mat + PARAM_PTR_NAME] = (double)PutText("mat");
	      WDB[mat + PARAM_PTR_FNAME] = (double)PutText(fname);

	      /* Put material name and atomic density */

	      WDB[mat + MATERIAL_PTR_NAME] = (double)PutText(tmpstr);
	      WDB[mat + MATERIAL_ADENS] = adens;

	      /* Reset pointer to composition */

	      WDB[mat + MATERIAL_PTR_COMP] = -1;

	      /* Loop over composition */

	      for (n = 0; n < nnuc; n++)
		{
		  /* Read ZAI and atomic density */
		  
		  if ((sz = fread(&zai, sizeof(long), 1, fp)) == 0)
		    Error(0, "Error in restart file");
		  
		  if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
		    Error(0, "Error in restart file");

		  /* Skip lost at beginning */

		  if (n > 0)
		    {
		      /* Create new item */
		  
		      iso = NewItem(mat + MATERIAL_PTR_COMP, 
				    COMPOSITION_BLOCK_SIZE);
		      
		      /* Put name */
		      
		      sprintf(ZAI, "%ld", zai);
		      WDB[iso + COMPOSITION_PTR_NUCLIDE] = PutText(ZAI);
		      
		      /* Put atomic denisty */
		      
		      WDB[iso + COMPOSITION_ADENS] = adens;
		    }	
		}	  

	      /***************************************************************/
	    }
	  else
	    {
	      /***************************************************************/

	      /***** Replace material composition ****************************/

	      /* Check pointer */

	      if (mat > VALID_PTR) 
		{
		  /* Put material-wise values */
		  
		  WDB[mat + MATERIAL_ADENS] = adens;
		  WDB[mat + MATERIAL_MDENS] = mdens;
		  WDB[mat + MATERIAL_BURNUP] = mbu;
		  
		  /* Get pointer to composition list */
		  
		  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) 
		      < VALID_PTR)
		    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

		  CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

		  /* Allocate memory for temporary nuclide pointers */

		  n = ListSize(iso);
		  ptr = (long *)Mem(MEM_ALLOC, n, sizeof(double));
		  
		  /* Store negative ZAI's in composition vector */
		  
		  n = 0;

		  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) 
		      < VALID_PTR)
		    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

		  while (iso > VALID_PTR)
		    {
		      /* Pointer to nuclide */
		      
		      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
		      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
		      
		      /* Store pointer */
		      
		      ptr[n++] = nuc;
		      
		      /* Replace pointer with ZAI */
		      
		      WDB[iso + COMPOSITION_PTR_NUCLIDE] = 
			RDB[nuc + NUCLIDE_ZAI];
		      
		      /* Reset density */
		      
		      WDB[iso + COMPOSITION_ADENS] = 0.0;
		      
		      /* Pointer to next */
		      
		      iso = NextItem(iso);		  
		    }
		}
	      else
		{
		  /* Not found, Skip material */
		  
		  n = (2*nnuc)*sizeof(double);
		  fseek(fp, n, SEEK_CUR);
		  
		  /* Reset nuclide count to avoid loop */
		  
		  nnuc = -1;
		  iso = -1;
		  ptr = NULL;
		}
	      
	      /* Read data */
	      
	      for (n = 0; n < nnuc; n++)
		{
		  /* Read ZAI and atomic density */
		  
		  if ((sz = fread(&zai, sizeof(long), 1, fp)) == 0)
		    Error(0, "Error in restart file");
		  
		  if ((sz = fread(&adens, sizeof(double), 1, fp)) == 0)
		    Error(0, "Error in restart file");
		  
		  /* Pointer to composition vector */
		  
		  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) 
		      < VALID_PTR)
		    iso = (long)RDB[mat + MATERIAL_PTR_COMP];

		  CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);
		  
		  /* Find nuclide */
		  
		  if ((iso = SeekList(iso, COMPOSITION_PTR_NUCLIDE, zai, 
				      SORT_MODE_ASCEND)) > VALID_PTR)
		    {
		      /* Put density */
		      
		      WDB[iso + COMPOSITION_ADENS] = adens;
		    }
		  else if (zai > 0)
		    Note(0, "Nuclide %ld not found in composition", zai);
		}
	      
	      /* Restore pointers */
	      
	      if (ptr != NULL)
		{
		  /* Loop over composition */
		  
		  n = 0;

		  if ((iso = (long)RDB[mat + MATERIAL_PTR_ORIG_NUC_COMP]) 
		      < VALID_PTR)
		    iso = (long)RDB[mat + MATERIAL_PTR_COMP];
		 
		  while (iso > VALID_PTR)
		    {
		      /* Put pointer */
		      
		      WDB[iso + COMPOSITION_PTR_NUCLIDE] = (double)ptr[n++];
		      
		      /* Pointer to next */
		      
		      iso = NextItem(iso);		  
		    }
		  
		  /* Free memory */
		  
		  Mem(MEM_FREE, ptr);
		}
	    
	      /***************************************************************/
	    }
	}
    }

  /* Check ok flag */

  if (ok == NO)
    Die(FUNCTION_NAME, "WTF?");

  /***************************************************************************/

  /* Calculate burnups, etc */

  SumDivCompositions();

  /* Close file */

  fclose(fp);

  /* Exit OK. */

  fprintf(out, "OK.\n\n");
}

/*****************************************************************************/
