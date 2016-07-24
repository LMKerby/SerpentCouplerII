#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readmaterialcomp.c                             */
/*                                                                           */
/* Created:       2013/06/06 (JLe)                                           */
/* Last modified: 2015/07/17 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Reads material compositions from an input file and overrides */
/*              the densities.                                               */
/*                                                                           */
/* Comments: - Tää on poistuva feature                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadMaterialComp:"

/*****************************************************************************/

void ReadMaterialComp()
{
  char *input, fname[MAX_STR], word[MAX_STR], pname[MAX_STR], **params;
  long loc0, loc1, mat, mat0, iso, iso0, nuc, i0, i, np, j, n, line, r, g, b;
  double val, sum;
  FILE *fp;

  /* Get pointer to file list */

  if ((loc0 = (long)RDB[DATA_PTR_COMP_FILE]) < VALID_PTR)
    return;

  fprintf(out, "Overriding initial material compositions...\n");

  /* Reset previous pointer */

  mat0 = -1;

  /* Reset counters for line number calculation */

  WDB[DATA_LINE_NUM_N0] = 0.0;
  WDB[DATA_LINE_NUM_NL0] = 1.0;

  /* Loop over list */

  while (RDB[loc0] > VALID_PTR)
    {
      /* Get file name */

      sprintf(fname, "%s", GetText(loc0));

      /* Check that file exists */

      if ((fp = fopen(fname, "r")) != NULL)
	fclose(fp);
      else
	{
	  /* File not found */
	  
	  Error(0, "Material composition file \"%s\" does not exist", fname);
	}

      /* Read input file */
  
      input = ReadTextFile(fname);

      /* Avoid compiler warning */
      
      params = NULL;

      /* Loop over file */
      
      i0 = 0;
      while ((i = NextWord(&input[i0], word)) > 0)
	{
	  /* update pointer */
	  
	  i0 = i0 + i;

	  /* Get line number for error messages */
	  
	  line = GetLineNumber(input, i0);

	  /* Look for material definition */

	  if (!strcasecmp(word, "mat"))
	    {
	      /* Copy parameter name */

	      strcpy (pname, word);

	      /* Read parameters */

	      params = GetParams(word, input, &np, &i0, 4, 4*MAX_ISOTOPES + 8, 
				 fname);

	      /* Read data */

	      j = 0;

	      /* Find material (try starting from previous) */

	      mat = mat0;
	      while (mat > VALID_PTR)
		{
		  /* Compare */

		  if (!strcmp(params[j], GetText(mat + MATERIAL_PTR_NAME)))
		    break;

		  /* Next */

		  mat = NextItem(mat);
		}

	      /* Find material (start from beginning) */

	      if (mat < VALID_PTR)
		{
		  mat = (long)RDB[DATA_PTR_M0];
		  while (mat > VALID_PTR)
		    {
		      /* Compare */
		      
		      if (!strcmp(params[j], GetText(mat + MATERIAL_PTR_NAME)))
			break;
		      
		      /* Next */
		      
		      mat = NextItem(mat);
		    }
		}

	      /* Check */

	      if (mat < VALID_PTR)
		Error(-1, pname, fname, line, "Material %s is not defined",
		      params[j]);
	      else
		j++;

	      /* Remember previous */

	      mat0 = NextItem(mat);
	      
	      /* Material density */

	      if (!strcmp(params[j], "sum"))
		{
		  /* Set value to -inf to calculate sum from composition */
		  
		  WDB[mat + MATERIAL_ADENS] = -INFTY;
		  
		  j++;
		}
	      else
		{
		  /* Read value */
		  
		  WDB[mat + MATERIAL_ADENS] = 
		    TestParam(pname, fname, line, params[j++], PTYPE_REAL,
			      -1000.0, 1000.0);
		}

	      /* Reset sum */

	      sum = 0.0;

	      /* Reset previous pointer */

	      iso0 = -1;

	      /* Loop over parameters */
	  
	      while (j < np)
		{
		  /* Check parameter */
		  
		  if (!strcmp(params[j], "tmp"))
		    {
		      /***** Temperature for Doppler-breadening **************/
		  
		      j++;
		      
		      /* Get temperature */
		      
		      WDB[mat + MATERIAL_DOPPLER_TEMP] = 
			TestParam(pname, fname, line, params[j++], 
				  PTYPE_REAL, 0.0, 100000.0);

		      /* Set option */
		      
		      WDB[DATA_USE_DOPPLER_PREPROCESSOR] = (double)YES;
		      
		      /*******************************************************/
		    }
		  if (!strcmp(params[j], "tms") || !strcmp(params[j], "ettm"))
		    {
		      /***** Temperature for TMS *****************************/
		  
		      j++;
		      
		      /* Get temperature */
		      
		      WDB[mat + MATERIAL_TMS_TMIN] = 
			TestParam(pname, fname, line, params[j++], 
				  PTYPE_REAL, 0.0, 100000.0);
		      
		      /* Copy to maximum */
		      
		      WDB[mat + MATERIAL_TMS_TMAX] = 
			RDB[mat + MATERIAL_TMS_TMIN];
		      
		      /* Set mode */
		      
		      WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;
		      WDB[mat + MATERIAL_TMS_MODE] = (double)YES;
		      
		      /*******************************************************/
		    }
		  else if (!strcmp(params[j], "rgb"))
		    {
		      /***** Material colour *********************************/

		      j++;
		      
		      /* Get r, b and g */
		      
		      r = TestParam(pname, fname, line, params[j++], 
				    PTYPE_INT, 0, 255);
		      
		      g = TestParam(pname, fname, line, params[j++], 
				    PTYPE_INT, 0, 255);

		      b = TestParam(pname, fname, line, params[j++], 
				    PTYPE_INT, 0, 255);

		      /* Set color */
		      
		      WDB[mat + MATERIAL_RGB] = b + 1000.0*g + 1000000.0*r;
		      
		      /*******************************************************/
		    }
		  else if (!strcmp(params[j], "vol"))
		    {
		      /***** Material volume *********************************/
		      
		      j++;
		      
		      /* Get volume */
		      
		      WDB[mat + MATERIAL_VOLUME_GIVEN] =
			TestParam(pname, fname, line, params[j++], PTYPE_REAL, 
				  0.0, INFTY);
		      
		      /*******************************************************/
		    }
		  else if (!strcmp(params[j], "fix"))
		    {
		      /***** Default library ID and temperature***************/

		      j++;
		  
		      /* Get default ID and temperature */
		      
		      WDB[mat + MATERIAL_DEFAULT_PTR_LIB_ID] = 
			PutText(params[j++]);
		      WDB[mat + MATERIAL_DEFAULT_TMP] = 
			TestParam(pname, fname, line, params[j++], PTYPE_REAL, 
				  0.0, INFTY);
		      
		      /*******************************************************/
		    }
		  else if (!strcmp(params[j], "mass"))
		    {
		      /***** Material mass ***********************************/

		      j++;
		      
		      /* Get mass */
		      
		      WDB[mat + MATERIAL_MASS_GIVEN] = 
			TestParam(pname, fname, line, params[j++], PTYPE_REAL, 
				  0.0, INFTY);

		      /*******************************************************/
		    }
		  else if (!strcmp(params[j], "burn"))
		    {
		      /***** Burnable material *******************************/
		      
		      j++;
		      
		      /* Set burn flag */
		      
		      SetOption(mat + MATERIAL_OPTIONS, OPT_BURN_MAT);
		      
		      /* Set burn sort flag and materials flag */
		      
		      WDB[mat + MATERIAL_BURN_SORT_FLAG] = 1.0;
		      WDB[DATA_BURN_MATERIALS_FLAG] = (double)YES;
		      
		      /* Get number of rings */
		      
		      WDB[mat + MATERIAL_BURN_RINGS] = 
			(double)TestParam(pname, fname, line, params[j++], 
					  PTYPE_INT, 0, 10000000);
		  
		      /*******************************************************/
		    }
		  else if (!strcmp(params[j], "moder"))
		    {
		      /***** Thermal scattering data *************************/

		      j++;
		      
		      /* Check number of parameters */
		      
		      if (j > np - 3)
			Error(mat, "Invalid number of parameters");
		      
		      /* Create new item (use the same structure as with */
		      /* the therm card) */
		      
		      WDB[mat + MATERIAL_PTR_SAB] = NULLPTR;
		      loc1 = NewItem(mat + MATERIAL_PTR_SAB, 
				     THERM_BLOCK_SIZE);

		      /* Read name */
		      
		      WDB[loc1 + THERM_PTR_ALIAS] = 
			(double)PutText(params[j++]);
		      
		      /* Read ZA */
		      
		      WDB[loc1 + THERM_ZA] =  
			(double)TestParam(pname, fname, line, params[j++], 
					  PTYPE_INT, 1001, 120000);
		      
		      /*******************************************************/
		    }
		  else if (!strcmp(params[j], "tft"))
		    {
		      /***** Minimum and maximum temperatures for TMS ********/

		      j++;
		      
		      /* Get minimum temperature */
		      
		      WDB[mat + MATERIAL_TMS_TMIN] = 
			TestParam(pname, fname, line, params[j++], 
				  PTYPE_REAL, 0.0, 100000.0);
		      
		      /* Get maximum temperature */
		      
		      WDB[mat + MATERIAL_TMS_TMAX] = 
			TestParam(pname, fname, line, params[j++], 
				  PTYPE_REAL, RDB[mat + MATERIAL_TMS_TMIN],
				  100000.0);
		      
		      /* Set mode */
		      
		      WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;
		      WDB[mat + MATERIAL_TMS_MODE] = (double)YES;
		      
		      /*******************************************************/
		    }
		  else 
		    {
		      /***** Composition *************************************/

		      /* Find nuclide in composition (start from previous) */

		      iso = iso0;
		      while (iso > VALID_PTR)
			{
			  /* Pointer to nuclide data */

			  nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
			  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY,nuc);

			  /* Compare */
			  
			  if (!strcmp(GetText(nuc + NUCLIDE_PTR_NAME),
				      params[j]))
			    break;

			  /* Next */

			  iso = NextItem(iso);
			}

		      /* Find nuclide in composition (start from beginning) */

		      if (iso < VALID_PTR)
			{			  
			  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
			  while (iso > VALID_PTR)
			    {
			      /* Pointer to nuclide data */
			      
			      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
			      CheckPointer(FUNCTION_NAME, "(nuc)", 
					   DATA_ARRAY,nuc);

			      /* Compare */
			  
			      if (!strcmp(GetText(nuc + NUCLIDE_PTR_NAME),
					  params[j]))
				break;
			      
			      /* Next */
			      
			      iso = NextItem(iso);
			    }
			}

		      /* Check pointer */
		      
		      if (iso < VALID_PTR)
			Error(-1, pname, fname, line, 
			      "Material %s has no nuclide %s in composition",
			      GetText(mat + MATERIAL_PTR_NAME), params[j]);
		      else
			j++;

		      /* Remember pointer */

		      iso0 = iso;
			      
		      /* Read fraction */
		  
		      val = TestParam(pname, fname, line, params[j++], 
				      PTYPE_REAL, -100.0, 1E+25);
		  
		      /* Put value */

		      WDB[iso + COMPOSITION_ADENS] = val;

		      /* Add to sum */
		      
		      sum = sum + val;
		      
		      /*******************************************************/
		    }
		}
	      
	      /* Set density if sum */

	      if (RDB[mat + MATERIAL_ADENS] == -INFTY)
		WDB[mat + MATERIAL_ADENS] = sum;

	      /* Calculate normalized fractions */
	      
	      IsotopeFractions(mat);
  	    }

	  /* Free parameter list */
      
	  if (np > 0)
	    for (n = 0; n < np + 1; n++)
	      Mem(MEM_FREE, params[n]);
	}

      /* Free memory */
  
      Mem(MEM_FREE, input);

      /* Next file */

      loc0++;
    }

  /* This must be called to get the divided compositions into material */
  /* structures */

  SumDivCompositions();
  
  fprintf(out, "OK.\n\n");
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
