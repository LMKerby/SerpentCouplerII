/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : checkrealistsum.c                              */
/*                                                                           */
/* Created:       2011/01/05 (JLe)                                           */
/* Last modified: 2012/10/04 (JLe)                                           */
/* Version:       2.1.9                                                      */
/*                                                                           */
/* Description: Check that reaction list sum matches total                   */
/*                                                                           */
/* Comments: - This is used to debug SampleReaction()                        */
/*           - Kutsu SampleMu:hun on väliaikainen ratkaisu                   */
/*           - Tästä tulee joskus errori kun E == 20 ja moodi on 0 tai 1     */
/*             (majorantti on nolla mutta reaktioiden summa ei)              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "CheckRealistSum:"

/*****************************************************************************/

void CheckReaListSum(long mat, double E, long kill, long id)
{
  long ptr, loc0, rea, n, nuc, lst1, rls1;
  double totxs, adens, sum, xs, Emin, Emax;

  /* Check material pointer */

  CheckPointer(FUNCTION_NAME, "(mat)", DATA_ARRAY, mat);

  /* Check energy */

  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);

  if (kill == YES)
    {
      fprintf(err, "\n\n%s Fatal error encountered...\n\n", FUNCTION_NAME); 
      fprintf(err, "mat = %s E = %E\n", GetText(mat + MATERIAL_PTR_NAME), E);
    }

  /* Get total cross section */

  ptr = (long)RDB[mat + MATERIAL_PTR_TOTXS];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  totxs = MacroXS(ptr, E, id);

  /* Check implicit capture mode */
  
  if (RDB[DATA_OPT_IMPL_CAPT] == YES)
    {
      /* Subtract absorption */
      
      ptr = (long)RDB[mat + MATERIAL_PTR_ABSXS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      totxs = totxs - MacroXS(ptr, E, id);
    }

  if (kill == YES)
    fprintf(err, "totxs    = %E (ures [%E %E])\n", totxs, 
	    RDB[ptr + REACTION_URES_EMIN], RDB[ptr + REACTION_URES_EMAX]);

  /* Get majorant */

  ptr = (long)RDB[DATA_PTR_MAJORANT];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  Die(FUNCTION_NAME, "Ton ptr:n paikalle pitää laittaa tyyppi..");

  xs = MajorantXS(ptr, E, id);

  if (kill == YES)
    fprintf(err, "majorant = %E (ures [%E %E])\n\n", xs, 
	    RDB[ptr + REACTION_URES_EMIN], RDB[ptr + REACTION_URES_EMAX]);

  /* Reset sum */

  sum = 0.0;

  if ((loc0 = (long)RDB[mat + MATERIAL_PTR_TOT_REA_LIST]) > VALID_PTR)
    {
      /***** Sum over nuclides and reactions *********************************/

      if (kill == YES)
	fprintf(err, "Calculating sum over nuclides and reactions:\n\n");
      
      /* Pointer to material total (toimii vaan neutroneille) */

      ptr = (long)RDB[mat + MATERIAL_PTR_TOTXS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      
      /* Get pointer to partial list */
      
      lst1 = (long)RDB[ptr + REACTION_PTR_PARTIAL_LIST];
      CheckPointer(FUNCTION_NAME, "(lst1)", DATA_ARRAY, lst1);

      /* Reset reaction pointer */

      rea = -1;
      
      /* Loop over reactions */
      
      while ((rls1 = NextReaction(lst1, &rea, &adens, &Emin, &Emax, id)) > 
	     VALID_PTR)
	{
	  /* Check reaction pointer */
	  
	  CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);
	  
	  /* Pointer to nuclide */
	  
	  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];
	  CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);

	  /* Pointer to partial reaction list */
      
	  ptr = (long)RDB[nuc + NUCLIDE_PTR_SAMPLE_REA_LIST];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Loop over partials */
	  
	  while (ptr > VALID_PTR)
	    {
	      /* Pointer to reaction data */
	      
	      rea = (long)RDB[ptr + RLS_DATA_PTR_REA];
	      CheckPointer(FUNCTION_NAME, "(rea)", DATA_ARRAY, rea);

	      /* Get cross section */

	      xs = MicroXS(rea, E, id);

	      /* Sample mu */
	      
	      if (xs > 0.0)
		for (n = 0; n < 2; n++)
		  SampleMu(rea, E, id);

	      /* Check energy */

	      if (E < RDB[rea + REACTION_EMIN])
		{
		  if (xs > 0.0)
		    Die(FUNCTION_NAME,
			"nonzero xs %E below boundary (%s %ld %E %E)",
			xs, GetText(nuc + NUCLIDE_PTR_NAME),
			(long)RDB[rea + REACTION_MT], E, 
			RDB[rea + REACTION_EMIN]);
		}
	      else
		{
		  /* Add partial cross section */

		  sum = sum + adens*xs;
		  
		  if (kill == YES)
		    {
		      fprintf(err, "nuc = %10s adens = %E mt = %4ld ",
			      GetText(nuc + NUCLIDE_PTR_NAME), adens,
			      (long)RDB[rea + REACTION_MT]);

		      if ((E > RDB[rea + REACTION_URES_EMIN]) &&
			  (E < RDB[rea + REACTION_URES_EMAX]))
			fprintf(err, "(u) ");
		      else
			fprintf(err, "    ");

		      fprintf(err, "Emin = %E ", RDB[rea + REACTION_EMIN]);
		      fprintf(err, "xs = %E sum = %E\n", adens*xs, sum);
		    }
		}

	      /* Next */

	      ptr = NextItem(ptr);
	    }
	}
    
      /***********************************************************************/
    }
  else
    Die(FUNCTION_NAME, "Material %s has no partial or total list",
	GetText(mat + MATERIAL_PTR_NAME));
  
  if (kill == YES)
    {
      fprintf(err, "\nOptimization mode = %ld\n", 
	      (long)RDB[DATA_OPTI_MODE]);
      fprintf(err, "Reconstruction tolerance = %E\n", RDB[DATA_ERG_TOL]);
      fprintf(err, "totxs    = %E  sum = %E diff = %E\n", totxs, sum,
	      sum/totxs - 1.0);

      Die(FUNCTION_NAME, "Ton rea:n paikalle pitää laittaa tyyppi..");
  
      if ((rea = (long)RDB[DATA_PTR_MAJORANT]) > VALID_PTR)
	fprintf(err, "majorant = %E diff = %E\n", 
		MajorantXS(rea, E, id), totxs/MajorantXS(rea, E, id) - 1.0);

      Die(FUNCTION_NAME, 
	  "Error in total or majorant, or reaction sampling failed");
    }

  /* Tähän kosahtaa helposti opt = 3 -moodissa */

  /* Compare sum and total */

  if ((fabs(sum - totxs) > 1E-10) && (fabs(sum/totxs - 1.0) > 1E-10))
    CheckReaListSum(mat, E, YES, id);
    
  /* Get majorant cross section and check */

  Die(FUNCTION_NAME, "Ton rea:n paikalle pitää laittaa tyyppi..");

  if ((rea = (long)RDB[DATA_PTR_MAJORANT]) > VALID_PTR)
    if (totxs/MajorantXS(rea, E, id) - 1.0 > 1E-6)
      CheckReaListSum(mat, E, YES, id);
}

/*****************************************************************************/
