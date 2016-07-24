#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : retrievecomposition.c                          */
/*                                                                           */
/* Created:       2012/08/24 (JLe)                                           */
/* Last modified: 2014/03/07 (JLe)                                           */
/* Version:       2.1.19                                                     */
/*                                                                           */
/* Description: Retrieves material composition from a binary work file       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "RetrieveComposition:"

/*****************************************************************************/

void RetrieveComposition(long mat0, long idx0)
{
  long nnuc, iso, mat, idx;
  char tmpstr[MAX_STR];
  double val;
  FILE *fp;
  size_t sz;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat0)", DATA_ARRAY, mat0);

  /* File name */

  sprintf(tmpstr, "%s.wrk", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for writing (vai append?) */

  fp = fopen(tmpstr, "r");
  
  /* Check pointer */

  if (fp == NULL)
    Die(FUNCTION_NAME, "Unable to open file for reading");

  /* Read loop */

  while ((sz = fread(&mat, sizeof(long), 1, fp)) > 0)
    {
      /* Read index */

      sz = fread(&idx, sizeof(long), 1, fp);

      /* Read number of nuclides */

      sz  = fread(&nnuc, sizeof(long), 1, fp); 

      /* Compare pointer and index */

      if ((mat == mat0) && (idx == idx0))
	{
	  /* Get atomic density, mass density and burnup */

	  sz  = fread(&val, sizeof(double), 1, fp);
	  WDB[mat + MATERIAL_ADENS] = val;

	  sz  = fread(&val, sizeof(double), 1, fp);
	  WDB[mat + MATERIAL_MDENS] = val;

	  sz  = fread(&val, sizeof(double), 1, fp);
	  WDB[mat + MATERIAL_BURNUP] = val;

	  /* Get pointer to composition */

	  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
	  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

	  /* Compare size */

	  if (nnuc != ListSize(iso))
	    Die(FUNCTION_NAME, "Mismatch in size");

	  /* Loop over composition */

	  while(iso > VALID_PTR)
	    {
	      /* Get atomic density */

	      sz  = fread(&val, sizeof(double), 1, fp);
	      WDB[iso + COMPOSITION_ADENS] = val;

	      /* Next */

	      iso = NextItem(iso);
	    }

	  /* Break loop */

	  break;
	}
      
      /* Skip material */
      
      sz = (nnuc + 3)*sizeof(double);
      fseek(fp, sz, SEEK_CUR);
    }

  /* Close file and exit */

  fclose(fp);  
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
