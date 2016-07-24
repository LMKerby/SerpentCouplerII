#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storecomposition.c                             */
/*                                                                           */
/* Created:       2012/08/24 (JLe)                                           */
/* Last modified: 2014/01/23 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Stores material composition in a binary work file            */
/*                                                                           */
/* Comments: - Tän vois nimetä writerestartfile.c:ksi                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreComposition:"

/*****************************************************************************/

void StoreComposition(long mat, double bu, double days)
{
  long nnuc, iso, nuc, n;
  char tmpstr[MAX_STR];
  double val;
  FILE *fp;

  /* Check if file is written */

  if ((long)RDB[DATA_WRITE_RESTART_FILE] == NO)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* File name */

  if ((long)RDB[DATA_RESTART_WRITE_PTR_FNAME] > VALID_PTR)
    sprintf(tmpstr, "%s", GetText(DATA_RESTART_WRITE_PTR_FNAME));
  else
    sprintf(tmpstr, "%s.wrk", GetText(DATA_PTR_INPUT_FNAME));

  /* Open file for writing (vai append?) */

  fp = fopen(tmpstr, "a");
  
  /* Check pointer */

  if (fp == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Get material name */

  sprintf(tmpstr, "%s", GetText(mat + MATERIAL_PTR_NAME));

  /* Write length of name */

  n = strlen(tmpstr);
  fwrite(&n, sizeof(long), 1, fp);

  /* Write material name, nominal burnup and burn time */

  fwrite(&tmpstr, sizeof(char), n, fp);
  fwrite(&bu, sizeof(double), 1, fp);
  fwrite(&days, sizeof(double), 1, fp);
  
  /* Get number of nuclides */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  nnuc = ListSize(iso);

  /* Write number of nuclides */

  fwrite(&nnuc, sizeof(long), 1, fp);

  /* Write atomic density, mass density and burnup */

  val = RDB[mat + MATERIAL_ADENS];
  fwrite(&val, sizeof(double), 1, fp);

  val = RDB[mat + MATERIAL_MDENS];
  fwrite(&val, sizeof(double), 1, fp);

  val = RDB[mat + MATERIAL_BURNUP];
  fwrite(&val, sizeof(double), 1, fp);

  /* Loop over composition */

  iso = (long)RDB[mat + MATERIAL_PTR_COMP];
  while (iso > VALID_PTR)
    {
      /* Pointer to nuclide */

      nuc =(long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

      /* Write ZAI */

      n = RDB[nuc + NUCLIDE_ZAI];
      fwrite(&n, sizeof(long), 1, fp);

      /* Write atomic density */
      
      val = RDB[iso + COMPOSITION_ADENS];
      fwrite(&val, sizeof(double), 1, fp);
	  
      /* Next nuclide in composition */

      iso = NextItem(iso);
    }

  /* Close file and exit */

  fclose(fp);  
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
