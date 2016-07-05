/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : linkreations.c                                 */
/*                                                                           */
/* Created:       2011/03/04 (JLe)                                           */
/* Last modified: 2014/06/07 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Links reactions to sources and detectors                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "LinkReactions:"

/*****************************************************************************/

void LinkReactions()
{
  long det, src, ptr, mat, rea, nuc, iso, mt, rfs;

  /***************************************************************************/

  /***** Detectors ***********************************************************/

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];  
  while (det > VALID_PTR)
    {
      /* Loop over reaction bins */

      ptr = (long)RDB[det + DET_PTR_RBINS];
      while (ptr > VALID_PTR)
	{
	  /* Get material pointer */
	  
	  mat = (long)RDB[ptr + DET_RBIN_PTR_MAT];

	  /* Get mt */

	  mt = (long)RDB[ptr + DET_RBIN_MT];

	  /* Check for TLE type */

	  if ((mat < VALID_PTR) && (mt != 0) &&
	      ((long)RDB[det + DET_PTR_SBINS] > VALID_PTR))
	    Error(det, "void response material not allowed with detector type");

	  /* Reset reaction pointer */

	  rea = -1;

	  /* Check mt */

	  if (mt == MT_MACRO_TOTXS)
	    {
	      if (mat > VALID_PTR)
		{
		  if ((long)RDB[det + DET_PARTICLE] == PARTICLE_TYPE_NEUTRON)
		    rea = (long)RDB[mat + MATERIAL_PTR_TOTXS];
		  else
		    rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
		}
	    }
	  else if (mt == MT_MACRO_ABSXS)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_ABSXS];
	    }
	  else if (mt == MT_MACRO_ELAXS)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_ELAXS];
	    }
	  else if (mt == MT_MACRO_INLPRODXS)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_INLPXS];
	    }
	  else if (mt == MT_MACRO_FISSXS)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_FISSXS];
	    }
	  else if (mt == MT_MACRO_FISSE)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_FISSE];
	    }
	  else if (mt == MT_MACRO_NSF)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_NSF];
	    }
	  else if (mt == MT_MACRO_RECOILE)
	    {
	      rea = -1;
	    }
	  else if (mt == MT_SOURCE_RATE)
	    {
	      rea = -1;
	    }
	  else if (mt == MT_MACRO_MAJORANT)
	    {
	      rea = -1;
	    }
	  else if (mt == MT_NEUTRON_DENSITY)
	    {
	      rea = -1;
	    }
	  else if (mt == MT_USER_DEFINED)
	    {
	      rea = -1;
	    }
	  else if ((mt == MT_PHOTON_DOSE) || ((mt > -249) && (mt < -200)))
	    {
	      rea = -1;
	    }
	  else if (mt == MT_MACRO_TOTPHOTXS)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_TOTPHOTXS];
	    }
	  else if (mt == MT_MACRO_TMP_MAJORANTXS)
            {
              if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_TMP_MAJORANTXS];
            }
	  else if (mt == MT_MACRO_HEATPHOTXS)
	    {
	      if (mat > VALID_PTR)
		rea = (long)RDB[mat + MATERIAL_PTR_HEATPHOTXS];
	    }
	  else if (mt < 0)
	    Error(det, "MT %ld not allowed in response function", mt);
	  else if (mt > 0)
	    {
	      /* Check material pointer */

	      if (mat < VALID_PTR)
		Error(det, "Response function with mt %ld must be associated with a material", mt);

	      /* Pointer to composition */

	      iso = (long)RDB[mat + MATERIAL_PTR_COMP];
	      CheckPointer(FUNCTION_NAME, "(iso)", DATA_ARRAY, iso);

	      /* Check number of nuclides */

	      if (NextItem(iso) > VALID_PTR)
		Error(det,  "Material %s used in response must consist of single nuclide", GetText(mat + MATERIAL_PTR_NAME));

	      /* Pointer to nuclide */

	      nuc = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];
	      CheckPointer(FUNCTION_NAME, "(nuc)", DATA_ARRAY, nuc);
	    
	      /* Check total */

	      if (mt == 1)
		rea = (long)RDB[nuc + NUCLIDE_PTR_TOTXS];
	      else
		{
		  /* Check for 102 to isomeric state */

		  if (mt == 1021)
		    {
		      /* Check branch data type */
		      /*
		      if ((long)RDB[nuc + NUCLIDE_BRA_TYPE] != BRA_TYPE_ENE)
			Error(det, "mt 1021 allowed only with energy-dependent branching ratios");
		      */
		      /* Set mt and state */
		      
		      mt = 102;
		      rfs = 1;
		    }
		  else
		    {
		      /* Reset state */

		      rfs = 0;
		    }

		  /* Find matching mt and rfs */
	      
		  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
		  while (rea > VALID_PTR)
		    {
		      /* Compare */
		      
		      if (((long)RDB[rea + REACTION_MT] == mt) && 
			  ((long)RDB[rea + REACTION_RFS] == rfs))
			break;
			  
		      /* Pointer to next */
		      
		      rea = NextItem(rea);
		    }

		  /* Check pointer */
	      
		  if (rea < VALID_PTR)
		    {
		      fprintf(stdout, 
			      "\nNuclide %s has the following reactions:\n\n",
			      GetText(nuc + NUCLIDE_PTR_NAME));
			  
		      /* Loop over reactions */
			  
		      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
		      while (rea > VALID_PTR)
			{
			  /* Check type */
			  
			  if (((long)RDB[rea + REACTION_TYPE] == 
			       REACTION_TYPE_PARTIAL) ||
			      ((long)RDB[rea + REACTION_TYPE] == 
			       REACTION_TYPE_SPECIAL))
			    {
			      if ((long)RDB[rea + REACTION_RFS] == 0)
				fprintf(stdout, "MT %-4ld : %s\n",
					(long)RDB[rea + REACTION_MT],
					ReactionMT((long)RDB[rea + REACTION_MT]));
			      else
				fprintf(stdout, "MT %-4ld : %s to isomeric state\n", 					(long)(10*RDB[rea + REACTION_MT] + 1),
					ReactionMT((long)RDB[rea + REACTION_MT]));
			    }
			      
			  /* Next reaction */
			      
			  rea = NextItem(rea);
			}
			  
		      Error(det, 
			    "Reaction mt %ld not found for response function",
			    (long)RDB[ptr + DET_RBIN_MT]);
		    }
		}
	    }
	
	  /* Check pointer */

	  if ((mt != 0) && (mt != MT_PHOTON_DOSE) && (rea < VALID_PTR) && 
	      (mat > VALID_PTR))
	    {
	      if (mt == MT_MACRO_RECOILE)
		Error(det, "Material entry must be void with mt %ld", mt);
	      else if ((mt != MT_MACRO_FISSXS) && (mt != MT_MACRO_FISSE))
		Error(det, "Reaction mt %ld not found for response function", 
		      mt);
	      else
		Note(det, "Reaction mt %ld not found for response function", 
		      mt);
	    }
	  
	  /* Set pointer */

	  WDB[ptr + DET_RBIN_PTR_REA] = (double)rea;

	  /* Next bin */

	  ptr = NextItem(ptr);
	}
      
      /* Next detector */

      det = NextItem(det);
    }

  /***************************************************************************/

  /***** Sources *************************************************************/

  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];  
  while (src > VALID_PTR)
    {
      /* Pointer to nuclide */

      if ((nuc = (long)RDB[src + SRC_PTR_XSDATA]) > VALID_PTR)
	{
	  /* Get mt */

	  mt = (long)RDB[src + SRC_PTR_REA];
	  
	  /* Find matching mt */
	  
	  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
	  if ((rea = SeekList(rea, REACTION_MT, mt, NO)) < VALID_PTR)
	    {
	      fprintf(stdout, 
		      "\nNuclide %s has the following reactions:\n\n",
		      GetText(nuc + NUCLIDE_PTR_NAME));
		  
	      /* Loop over reactions */
	      
	      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
	      while (rea > VALID_PTR)
		{
		  /* Check type */
		  
		  if (((long)RDB[rea + REACTION_TYPE] == 
		       REACTION_TYPE_PARTIAL) ||
		      ((long)RDB[rea + REACTION_TYPE] == 
		       REACTION_TYPE_SPECIAL))
		    fprintf(stdout, "MT %-3ld : %s\n",
			    (long)RDB[rea + REACTION_MT],
			    ReactionMT((long)RDB[rea + REACTION_MT]));
		  
		  /* Next reaction */
		  
		  rea = NextItem(rea);
		}
	      
	      Error(src, "Reaction mt %ld not found for source distribution", 
		    mt);
	    }
	  
	  /* Put pointer */

	  WDB[src + SRC_PTR_REA] = (double)rea;

	  /* Check that reaction has distribution data */

	  if ((long)RDB[rea + REACTION_PTR_ERG] < VALID_PTR)
	    Error(src, "Reaction mt %ld has no distribution data", mt);
	}
      
      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/
}  

/*****************************************************************************/
