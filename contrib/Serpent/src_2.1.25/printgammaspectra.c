/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printgammaspectra.c                            */
/*                                                                           */
/* Created:       2015/05/18 (JLe)                                           */
/* Last modified: 2015/05/22 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Prints gamma spectra for different materials                 */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintGammaSpectra:"

/*****************************************************************************/

void PrintGammaSpectra()
{
  long mat, nuc, src, ptr;
  char outfile[MAX_STR];
  double cum;
  FILE *fp;

  fprintf(out, "Printing photon spectra in file...\n");

  /* File name */

  sprintf(outfile, "%s_gsrc.m", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file */

  if ((fp = fopen(outfile, "w")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Print columns */

  fprintf(fp, "\n%% Column 1: Nuclide ZAI\n");
  fprintf(fp, "%% Column 2: Specific intensity (photons per decay)\n");
  fprintf(fp, "%% Column 3: Total emission rate (photons/s)\n");
  fprintf(fp, "%% Column 4: Cumulative fraction of material total\n");
  fprintf(fp, "%% Column 5: Emission line energy\n");
  fprintf(fp, "%% Column 6: Relative intensity (photons per decay)\n");
  fprintf(fp, "%% Column 7: Cumulative fraction of nuclide total\n");

  /* Newline */

  fprintf(fp, "\n");

  /* Loop over materials */
      
  mat = (long)RDB[DATA_PTR_M0];
  while (mat > VALID_PTR)
    {
      /* Check decay source */

      if ((src = (long)RDB[mat + MATERIAL_PTR_DECAY_SRC]) < VALID_PTR)
	{
	  /* Next material */
	  
	  mat = NextItem(mat);

	  /* Cycle loop */

	  continue;
	}

      /* Print material name */

      fprintf(fp, "mat_%s = [\n", GetText(mat + MATERIAL_PTR_NAME));

      /* Loop over source */

      while (src > VALID_PTR)
	{
	  /* Pointer to nuclide */

	  nuc = (long)RDB[src + SRC_DECCAY_PTR_NUCLIDE];
	  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	  /* Pointer to line spectra */

	  ptr = (long)RDB[src + SRC_DECCAY_PTR_SPEC];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  fprintf(fp, "\n%% --- %s :\n\n", 
		  ZAItoIso((long)RDB[nuc + NUCLIDE_ZAI], 1));
				
	  /* Reset sum */

	  cum = 0.0;
	  
	  /* Loop over spectra */

	  while (ptr > VALID_PTR)
	    {
	      /* Print ZAI, specific and total intensity and fraction */
	      
	      fprintf(fp, "%6ld %11.5E %11.5E %11.5E ", 
		      (long)RDB[nuc + NUCLIDE_ZAI],
		      RDB[nuc + NUCLIDE_SPEC_PHOTON_I], 
		      RDB[src + SRC_DECCAY_I], 
		      RDB[src + SRC_DECCAY_CUM_P]);

	      /* Add to cumulative */

	      cum = cum + RDB[ptr + PHOTON_LINE_SPEC_RI];

	      /* Print energy, relative intensity and total */

	      fprintf(fp, "%11.5E %11.5E %11.5E", RDB[ptr + PHOTON_LINE_SPEC_E],
		      RDB[ptr + PHOTON_LINE_SPEC_RI],
		      cum/RDB[nuc + NUCLIDE_SPEC_PHOTON_I]);

	      /* Newline */

	      fprintf(fp, "\n");

	      /* Next line */

	      ptr = NextItem(ptr);
	    }

	  /* Next */

	  src = NextItem(src);
	}

      fprintf(fp, "];\n\n");
      
      /* Print volume and total */

      fprintf(fp, "mat_%s_vol = %11.5E;\n", GetText(mat + MATERIAL_PTR_NAME),
	      RDB[mat + MATERIAL_VOLUME]);

      fprintf(fp, "mat_%s_tot = %11.5E;\n\n", GetText(mat + MATERIAL_PTR_NAME),
	      RDB[mat + MATERIAL_PHOTON_SRC_RATE]);

      /* Next material */

      mat = NextItem(mat);
    }      

  /* Print total */

  fprintf(fp, "tot = %11.5E;\n\n", RDB[DATA_TOT_PHOTON_SRC_RATE]);

  /* Close file */

  fclose(fp);

  fprintf(out, "OK.\n\n");
}
 
/*****************************************************************************/
