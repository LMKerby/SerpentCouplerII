/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processdetectors.c                             */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2015/05/25 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Processes detector definitions                               */
/*                                                                           */
/* Comments: - Detector cells and universes are processed in                 */
/*             creategeometry.c                                              */
/*           - TODO: se automaattinen volume-juttu                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"
#include "photon_attenuation.h"

#define FUNCTION_NAME "ProcessDetectors:"

/*****************************************************************************/

void ProcessDetectors()
{
  long det, ene, ptr, mat, uni, lat, cell, surf, tme, tot, n1, n2, msh1, msh2;
  long ebins, ubins, cbins, mbins, lbins, rbins, zbins, ybins, xbins, tbins;
  long mt, ne, n, loc0, loc1, loc2, umsh, sflag, idx, i0;
  double sum;
  char str[MAX_STR];

  /* Avoid compiler warning (photon_attenuation.h sisältää dataa jota */
  /* käytetään useammassa aliohjelmassa) */

  n = idx0[0][0];
  sum = dat0[0][0];

  /* Adjust coincident points */

  ne = (long)idx1[48][0] + (long)idx1[48][1];
  
  for (n = 1; n < ne; n++)
    if (dat1[n][1] == dat1[n - 1][1])
      dat1[n][1] = dat1[n][1] + 1E-11;

  /***************************************************************************/

  /***** Structure for criticality source detector ***************************/

  /* Check if file name is given */

  if ((long)RDB[DATA_PTR_CRIT_SRC_DET]  > VALID_PTR)
    {
      /* Allocate memory for detector structure (not included in list) */

      det = ReallocMem(DATA_ARRAY, DET_BLOCK_SIZE);

      /* Put file name */

      WDB[det + DET_WRITE_PTR_FILE] = RDB[DATA_PTR_CRIT_SRC_DET];
      
      /* Allocate memory for buffer */
      
      ptr = ReallocMem(DATA_ARRAY, SRC_FILE_BUF_SIZE*SRC_BUF_BLOCK_SIZE);

      /* Put pointer and buffer size */
      
      WDB[det + DET_WRITE_PTR_BUF] = (double)ptr;
      WDB[det + DET_WRITE_BUF_SZ] = (double)SRC_FILE_BUF_SIZE;

      /* Reset index */

      WDB[det + DET_WRITE_BUF_IDX] = 0.0;

      /* Put fraction */

      WDB[det + DET_WRITE_PROB] = 1.0;

      /* Remove old file */
      
      remove(GetText(det + DET_WRITE_PTR_FILE));

      /* Put pointer */

      WDB[DATA_PTR_CRIT_SRC_DET] = (double)det;
    }

  /* Exit if detectors are not defined */

  if ((long)RDB[DATA_PTR_DET0] < VALID_PTR)
    return;

  /***************************************************************************/

  /***** Link pointers and allocate memory ***********************************/

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];  
  while (det > VALID_PTR)
    {
      /* Number of mesh bins */
  
      if ((msh1 = (long)RDB[det + DET_PTR_MESH]) > VALID_PTR)
	{
	  xbins = (long)RDB[msh1 + MESH_N0];
	  ybins = (long)RDB[msh1 + MESH_N1];
	  zbins = (long)RDB[msh1 + MESH_N2];
	}
      else
	{
	  xbins = 1;
	  ybins = 1;
	  zbins = 1;
	}
      
      /***********************************************************************/
      
      /***** Link energy grids to detectors **********************************/

      /* Check pointer */
      
      if ((long)RDB[det + DET_PTR_EGRID] > VALID_PTR)
	{
	  /* Find grid */

	  ene = RDB[DATA_PTR_ENE0];
	  if ((ene = SeekListStr(ene, ENE_PTR_NAME, 
				 GetText(det + DET_PTR_EGRID))) < VALID_PTR)
	    Error(det, "Energy grid %s in detector %s is not defined", 
		  GetText(det + DET_PTR_EGRID), 
		  GetText(det + DET_PTR_NAME));
	  
	  /* Check if same grid is used for group constant generation */

	  if (ene == (long)RDB[DATA_ERG_FG_PTR_PREDEF])
	    Note(det, "Energy grid %s is used for both group constant generation\n and detector calculation (detector %s)", GetText(ene + ENE_PTR_NAME),
		 GetText(det + DET_PTR_NAME));

	  /* Set pointer */

	  WDB[det + DET_PTR_EGRID] = RDB[ene + ENE_PTR_GRID];
	  
	  /* Set number of bins */

	  ebins = (long)RDB[ene + ENE_NB];

	  /* Check */
	  
	  if (ebins < 1)
	    Error(ene, "Error in bin structure");
	}
      else
	ebins = 1;

      /***********************************************************************/

      /***** Link time bins to detectors *************************************/

      /* Check pointer */
      
      if ((long)RDB[det + DET_PTR_TME] > VALID_PTR)
	{
	  /* Find bins */
	  
	  tme = RDB[DATA_PTR_TME0];
	  if ((tme = SeekListStr(tme, TME_PTR_NAME, 
				 GetText(det + DET_PTR_TME))) < VALID_PTR)
	    Error(det, "Time binning %s in detector %s is not defined", 
		  GetText(det + DET_PTR_TME), 
		  GetText(det + DET_PTR_NAME));
	  
	  /* Set pointer */
	  
	  WDB[det + DET_PTR_TME] = RDB[tme + TME_PTR_BINS];
	  
	  /* Set number of bins */

	  tbins = (long)RDB[tme + TME_NB];

	  /* Stop tracks at outer boundary */

	  WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;
	}
      else
	tbins = 1;

      /***********************************************************************/
      
      /***** Link materials to response functions ****************************/

      /* Check pointer */

      if ((ptr = (long)RDB[det + DET_PTR_RBINS]) < VALID_PTR)
	{
	  /* Create response bin for flux */
	  
	  ptr = NewItem(det + DET_PTR_RBINS, DET_RBIN_BLOCK_SIZE);
	      
	  /* Put mt and material pointer */
	  
	  WDB[ptr + DET_RBIN_MT] = 0.0;
	  WDB[ptr + DET_RBIN_PTR_MAT] = NULLPTR;
	  
	  /* Set void mode */
	  
	  WDB[ptr + DET_RBIN_VOID_MODE] = (double)YES;
	  
	  /* Set counter */

	  rbins = 1;
	}
      else
	{
	  /* Reset counter */

	  rbins = 0;
	  
	  /* Loop over reaction bins */
 
	  while (ptr > VALID_PTR)
	    {
	      /* Get mt */

	      mt = (long)RDB[ptr + DET_RBIN_MT];

	      /* Put conversion flag for photon dose rate */

	      if (mt == MT_PHOTON_DOSE)
		WDB[ptr + DET_RBIN_CONVERT] = (double)DET_CONVERT_GDOSE;

	      /* Check type */

	      if (mt == MT_USER_DEFINED)
		{
		  /* Get number of values */

		  ne = (long)RDB[ptr + DET_RBIN_USDEF_NP];
		  
		  /* Get pointer to energy array */

		  loc0 = (long)RDB[ptr + DET_RBIN_USDEF_PTR_E];
		  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

		  /* Check that energy values are in ascending order */

		  for (n = 1; n < ne ; n++)
		    if (RDB[loc0 + n - 1] >= RDB[loc0 + n])
		      Error(det, 
			    "User-defined energies not in ascending order");

		  /* Get pointer to values */

		  loc0 = (long)RDB[ptr + DET_RBIN_USDEF_PTR_F];
		  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

		  /* Check values for log-interpolations */

		  if (((long)RDB[ptr + DET_RBIN_USDEF_INTT] == 3) ||
		      ((long)RDB[ptr + DET_RBIN_USDEF_INTT] == 5))
		    for (n = 0; n < ne ; n++)
		      if (RDB[loc0 + n] <= 0)
			Error(det, 
			      "Value %1.5E not allowed with log interpolation",
			      RDB[loc0 + n]);
		}
	      else if ((mt > -249) && (mt < -200))
		{
		  /* Photon attenuation, get index */

		  idx = -mt - 201;

		  /* Get number of energy points and index to first */

		  ne = (long)idx1[idx][1];
		  i0 = (long)idx1[idx][0];

		  /* Allcate memory for energy array */

		  loc0 = ReallocMem(DATA_ARRAY, ne);
		  WDB[ptr + DET_RBIN_USDEF_PTR_E] = (double)loc0;

		  /* Read data */

		  for (n = 0; n < ne; n++)
		    WDB[loc0++] = dat1[i0 + n][1];

		  /* Allcate memory for data array */

		  loc0 = ReallocMem(DATA_ARRAY, ne);
		  WDB[ptr + DET_RBIN_USDEF_PTR_F] = (double)loc0;

		  /* Read data */

		  for (n = 0; n < ne; n++)
		    WDB[loc0++] = dat1[i0 + n][3];

		  /* Put number of points and interpolation scheme */

		  WDB[ptr + DET_RBIN_USDEF_NP] = (double)ne;
		  WDB[ptr + DET_RBIN_USDEF_INTT] = 5.0;

		  /* Set null material pointer */
		  
		  WDB[ptr + DET_RBIN_PTR_MAT] = NULLPTR;

		  /* Set void mode */

		  WDB[ptr + DET_RBIN_VOID_MODE] = (double)YES;

		  /* Override mt */

		  WDB[ptr + DET_RBIN_MT] = (double)MT_USER_DEFINED;

		  /* Put conversion option */

		  WDB[ptr + DET_RBIN_CONVERT] = (double)DET_CONVERT_GDOSE;
		}
	      else if ((!strcmp(GetText(ptr + DET_RBIN_PTR_MAT), "void")) ||
		       (mt == MT_NEUTRON_DENSITY))
		{
		  /* Set null pointer */
		  
		  WDB[ptr + DET_RBIN_PTR_MAT] = NULLPTR;

		  /* Set void mode */

		  WDB[ptr + DET_RBIN_VOID_MODE] = (double)YES;
		}
	      else
		{
		  /* Find material */

		  mat = RDB[DATA_PTR_M0];
		  if ((mat = SeekListStr(mat,MATERIAL_PTR_NAME, 
					 GetText(ptr + DET_RBIN_PTR_MAT))) 
		      < VALID_PTR)
		    Error(det,
			  "Material %s in detector %s response is not defined", 
			  GetText(ptr + DET_RBIN_PTR_MAT), 
			  GetText(det + DET_PTR_NAME));
		  
		  /* Set pointer */
		  
		  WDB[ptr + DET_RBIN_PTR_MAT] = (double)mat;
		  
		  /* Reset void mode */

		  WDB[ptr + DET_RBIN_VOID_MODE] = (double)NO;

		  /* Set used-flags */
		  
		  SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
		}

	      /* Update number of bins */
	      
	      rbins++;
	      
	      /* Next bin */
	      
	      ptr = NextItem(ptr);
	    }
	}

      /***********************************************************************/

      /***** Link universes to detectors *************************************/

      /* Check pointer */

      if ((ptr = (long)RDB[det + DET_PTR_UBINS]) < VALID_PTR)
	ubins = 1;
      else
	{
	  /* Reset counter */

	  ubins = 0;

	  /* Loop over universe bins */

	  while (ptr > VALID_PTR)
	    {
	      /* Find universe */
	      
	      uni = RDB[DATA_PTR_U0];
	      if ((uni = SeekListStr(uni, UNIVERSE_PTR_NAME, 
				     GetText(ptr + DET_UBIN_PTR_UNI))) 
		  < VALID_PTR)
		Error(det,"Universe %s in detector %s is not defined", 
		      GetText(ptr + DET_UBIN_PTR_UNI), 
		      GetText(det + DET_PTR_NAME));
	      
	      /* Set pointer */
	      
	      WDB[ptr + DET_UBIN_PTR_UNI] = (double)uni;
	      
	      /* Update number of bins */
	      
	      ubins++;
	      
	      /* Next bin */
	      
	      ptr = NextItem(ptr);
	    }
	}
      
      /***********************************************************************/

      /***** Link lattices to detectors **************************************/

      /* Check pointer */

      if ((ptr = (long)RDB[det + DET_PTR_LBINS]) < VALID_PTR)
	lbins = 1;
      else
	{
	  /* Find lattice */
	      
	  lat = RDB[DATA_PTR_L0];
	  if ((lat = SeekListStr(lat, LAT_PTR_NAME, 
				 GetText(ptr + DET_LBIN_PTR_LAT))) 
	      < VALID_PTR)
	    Error(det,"Lattice %s in detector %s is not defined", 
		  GetText(ptr + DET_LBIN_PTR_LAT), 
		  GetText(det + DET_PTR_NAME));
	  
	  /* Set pointer */
	      
	  WDB[ptr + DET_LBIN_PTR_LAT] = (double)lat;
	  
	  /* Set number of bins */
	  
	  lbins = (long)RDB[lat + LAT_NTOT];
	}
      
      /***********************************************************************/

      /***** Link materials to detectors *************************************/

      /* Check pointer */

      if ((ptr = (long)RDB[det + DET_PTR_MBINS]) < VALID_PTR)
	mbins = 1;
      else
	{
	  /* Reset counter */

	  mbins = 0;

	  /* Loop over bins */
      
	  while (ptr > VALID_PTR)
	    {
	      /* Find material */
	      
	      mat = RDB[DATA_PTR_M0];
	      if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, 
				     GetText(ptr + DET_MBIN_PTR_MAT))) 
		  < VALID_PTR)
		{
		  /* Check for scoring in fissile regions only */

		  if (!strcmp(GetText(ptr + DET_MBIN_PTR_MAT), "fiss"))
		    {
		      /* Set flag */

		      WDB[det + DET_SCORE_FISS_REG_ONLY] = (long)YES;
		      
		      /* Reset material pointer */

		      WDB[det + DET_PTR_MBINS] = NULLPTR;

		      /* Set number of bins */

		      mbins = 1;

		      /* Break loop */

		      break;
		    }
		  else
		    Error(det,"Material %s in detector %s is not defined", 
			  GetText(ptr + DET_MBIN_PTR_MAT), 
			  GetText(det + DET_PTR_NAME));
		}
	      else
		{
		  /* Set pointer */
	      
		  WDB[ptr + DET_MBIN_PTR_MAT] = (double)mat;
		  
		  /* Add to material structure */

		  loc0 = NewItem(mat + MATERIAL_PTR_DETBIN, DETBIN_BLOCK_SIZE);

		  /* Put pointer and bin */

		  WDB[loc0 + DETBIN_PTR_DET] = det;
		  WDB[loc0 + DETBIN_BIN] = mbins;
		  
		  /* Update number of bins */
		  
		  mbins++;
		}

	      /* Next bin */
	      
	      ptr = NextItem(ptr);
	    }
	}
      
      /***********************************************************************/

      /***** Link cells to detectors *****************************************/

      /* Check pointer */

      if ((ptr = (long)RDB[det + DET_PTR_CBINS]) < VALID_PTR)
	cbins = 1;
      else
	{
	  /* Check pointer to UMSH geometry */

	  if ((long)RDB[ptr + DET_CBIN_UMSH_PTR_UMSH] > VALID_PTR)
	    {
	      /* Check other cells */

	      if (NextItem(ptr) > VALID_PTR)
		Error(det, "Multiple cell bins not allowed with UMSH type");

	      /* Find UMSH geometry */
	      
	      umsh = RDB[DATA_PTR_UMSH0];
	      if ((umsh = SeekListStr(umsh, UMSH_PTR_NAME, 
				      GetText(ptr + DET_CBIN_UMSH_PTR_UMSH))) 
		  < VALID_PTR)
		Error(det,"UMSH geometry %s in detector %s is not defined", 
		      GetText(ptr + DET_CBIN_UMSH_PTR_UMSH), 
		      GetText(det + DET_PTR_NAME));
	      
	      /* Set pointer */
	      
	      WDB[ptr + DET_CBIN_UMSH_PTR_UMSH] = (double)umsh;
	  
	      /* Get pointers to cells and bins */
	      
	      loc0 = (long)RDB[ptr + DET_CBIN_UMSH_PTR_CELLS];
	      CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);
	      
	      loc1 = (long)RDB[ptr + DET_CBIN_UMSH_PTR_BINS];
	      CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
	      
	      /* Get pointer to tet cells */
	      
	      ptr = (long)RDB[umsh + UMSH_PTR_IFC];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	      loc2 = (long)RDB[ptr + IFC_PTR_TET_MSH];
	      CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);
	      
	      /* Reset number of bins */
	      
	      cbins = -1;

	      /* Loop over data */
	      
	      while ((n1 = (long)RDB[loc0++]) > 0)
		{
		  /* Check cell index */
		  
		  if (n1 > (long)RDB[umsh + UMSH_N_ORIG_CELLS])
		    Error(det, "Invalid cell index %ld\n", n1);
		  
		  /* Get bin index */
		  
		  n2 = (long)RDB[loc1++];
		  
		  /* Get pointer to tet cell data */
		  
		  ptr = ListPtr(loc2, n1 - 1);
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  
		  /* Put stat index */
		  
		  WDB[ptr + IFC_TET_MSH_STAT_IDX] = (double)(n2 - 1);
		  
		  /* Compare to bins */
		  
		  if (n2 > cbins)
		    cbins = n2;
		}
	      
	      /* Check number of bins */

	      if ((cbins < 1) || 
		  (cbins > (long)RDB[umsh + UMSH_N_ORIG_CELLS]))
		Die(FUNCTION_NAME, "Something wrong here");
	    }
	  else
	    {
	      /* Reset counter and super-imposed flag */

	      cbins = 0;
	      sflag = -1;
	      
	      /* Loop over cell bins */
	      
	      while (ptr > VALID_PTR)
		{
		  /* Check type */

		  if ((long)RDB[ptr + DET_CBIN_UMSH_PTR_UMSH] > VALID_PTR)
		    Error(det, "Multiple cell bins not allowed with UMSH type");
		  
		  /* Find cell */
		  
		  cell = RDB[DATA_PTR_C0];
		  if ((cell = SeekListStr(cell, CELL_PTR_NAME, 
					  GetText(ptr + DET_CBIN_PTR_CELL))) 
		      < VALID_PTR)
		    Error(det,"Cell %s in detector %s is not defined", 
			  GetText(ptr + DET_CBIN_PTR_CELL), 
			  GetText(det + DET_PTR_NAME));
		  
		  /* Set pointer */
		  
		  WDB[ptr + DET_CBIN_PTR_CELL] = (double)cell;
		  
		  /* Pointer to universe */
		  
		  uni = (long)RDB[cell + CELL_PTR_UNI];
		  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
		  
		  /* Check universe type and set super-imposed flag */
		  
		  if ((long)RDB[uni + UNIVERSE_TYPE] == UNIVERSE_TYPE_SUPER)
		    WDB[ptr + DET_CBIN_SUPER_CELL] = (double)YES;
		  else
		    WDB[ptr + DET_CBIN_SUPER_CELL] = (double)NO;
		  
		  /* Check flag */

		  if (sflag < 0)
		    sflag = (long)RDB[ptr + DET_CBIN_SUPER_CELL];
		  else if (sflag != (long)RDB[ptr + DET_CBIN_SUPER_CELL])
		    Error(det, 
			  "Super-imposed and physical cell boms not allowed in same detector");

		  /* Add to material structure */

		  loc0 = NewItem(cell + CELL_PTR_DETBIN, DETBIN_BLOCK_SIZE);
		  
		  /* Put pointer and bin */

		  WDB[loc0 + DETBIN_PTR_DET] = det;
		  WDB[loc0 + DETBIN_BIN] = cbins;

		  /* Update number of bins */
		  
		  cbins++;
		  
		  /* Next bin */
		  
		  ptr = NextItem(ptr);
		}
	    }
	}
      
      /***********************************************************************/

      /***** Link dividers and multipliers ***********************************/

      /* Check pointer */

      if ((long)RDB[det + DET_PTR_MUL] > VALID_PTR)
	{
	  /* Find detector */
	  
	  ptr = RDB[DATA_PTR_DET0];
	  if ((ptr = SeekListStr(ptr, DET_PTR_NAME,
				 GetText(det + DET_PTR_MUL)))
	      < VALID_PTR)
	    Error(det, "Multiplier / divider detector %s is not defined", 
		  GetText(det + DET_PTR_MUL));
	  
	  /* Set pointer */
	  
	  WDB[det + DET_PTR_MUL] = (double)ptr;
	}
      
      /***********************************************************************/

      /***** Link adjoints ***************************************************/

      /* Check pointer */

      if ((long)RDB[det + DET_PTR_ADJOINT] > VALID_PTR)
	{
	  /* Find detector */
	  
	  ptr = RDB[DATA_PTR_DET0];
	  if ((ptr = SeekListStr(ptr, DET_PTR_NAME,
				 GetText(det + DET_PTR_ADJOINT)))
	      < VALID_PTR)
	    Error(det,"Adjoint detector %s is not defined", 
		  GetText(det + DET_PTR_ADJOINT));
	  
	  /* Set pointer */
	  
	  WDB[det + DET_PTR_ADJOINT] = (double)ptr;

	  /* Switch history lists on */
	  
	  WDB[DATA_HIST_LIST_SIZE] = fabs(RDB[DATA_HIST_LIST_SIZE]);
	}
      
      /***********************************************************************/

      /***** Link surfaces to detectors **************************************/

      /* Check pointer */

      if ((ptr = (long)RDB[det + DET_PTR_SBINS]) > VALID_PTR)
	{
	  /* Find surface */
	      
	  surf = RDB[DATA_PTR_S0];
	  if ((surf = SeekListStr(surf, SURFACE_PTR_NAME, 
				  GetText(ptr + DET_SBIN_PTR_SURF))) 
	      < VALID_PTR)
	    Error(det,"Surface %s in detector %s is not defined", 
		  GetText(ptr + DET_SBIN_PTR_SURF), 
		  GetText(det + DET_PTR_NAME));
	  
	  /* Set pointer */
	  
	  WDB[ptr + DET_SBIN_PTR_SURF] = (double)surf;

	  /* Check invalid options (Tää ei tarkista kaikkea) */

	  if ((long)RDB[det + DET_PTR_UBINS] > VALID_PTR)
	    Error(det, "Universe bins not allowed with current detector");
	  if ((long)RDB[det + DET_PTR_LBINS] > VALID_PTR)
	    Error(det, "Lattice bins not allowed with current detector");
	  if ((long)RDB[det + DET_PTR_MBINS] > VALID_PTR)
	    Error(det, "Material bins not allowed with current detector");
	  if ((long)RDB[det + DET_PTR_CBINS] > VALID_PTR)
	    Error(det, "Cell bins not allowed with current detector");

	  /* Stop tracks at outer boundary */

	  WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;
	}

      /***********************************************************************/
      
      /***** Detector source files *******************************************/

      /* Check if file name is given */

      if ((long)RDB[det + DET_WRITE_PTR_FILE] > VALID_PTR)
	{
	  /* Allocate memory for buffer */

	  ptr = ReallocMem(DATA_ARRAY, SRC_FILE_BUF_SIZE*SRC_BUF_BLOCK_SIZE);

	  /* Put pointer and buffer size */

	  WDB[det + DET_WRITE_PTR_BUF] = (double)ptr;
	  WDB[det + DET_WRITE_BUF_SZ] = (double)SRC_FILE_BUF_SIZE;

	  /* Reset index */

	  WDB[det + DET_WRITE_BUF_IDX] = 0.0;

	  /* Remove old file */

	  remove(GetText(det + DET_WRITE_PTR_FILE));
	}
      
      /***********************************************************************/

      /***** Normalize direction vector **************************************/

      /* Get square sum */

      if ((sum = RDB[det + DET_DIRVEC_U]*RDB[det + DET_DIRVEC_U] +
	   RDB[det + DET_DIRVEC_V]*RDB[det + DET_DIRVEC_V] +
	   RDB[det + DET_DIRVEC_W]*RDB[det + DET_DIRVEC_W]) > 0.0)
	{
	  /* Take square root */

	  sum = sqrt(sum);

	  /* Normalize components */

	  WDB[det + DET_DIRVEC_U] = RDB[det + DET_DIRVEC_U]/sum;
	  WDB[det + DET_DIRVEC_V] = RDB[det + DET_DIRVEC_V]/sum;
	  WDB[det + DET_DIRVEC_W] = RDB[det + DET_DIRVEC_W]/sum;
	}

      /***********************************************************************/

      /***** Allocate memory for statistics **********************************/
    
      /* Set final bin sizes */

      WDB[det + DET_N_EBINS] = (double)ebins;
      WDB[det + DET_N_UBINS] = (double)ubins;
      WDB[det + DET_N_CBINS] = (double)cbins;
      WDB[det + DET_N_MBINS] = (double)mbins;
      WDB[det + DET_N_LBINS] = (double)lbins;
      WDB[det + DET_N_RBINS] = (double)rbins;
      WDB[det + DET_N_TBINS] = (double)tbins;

      /* Calculate total size */

      tot = ebins*ubins*cbins*mbins*lbins*rbins*zbins*ybins*xbins*tbins;

      /* Check maximum size */

      if (tot > 10000000)
	Error(det, "Total number of bins exceeds maximum");

      /* Check zero */

      if (tot < 1)
	Die(FUNCTION_NAME, "Error in bin sizes");

      /* Put total number fo bins */

      WDB[det + DET_N_TOT_BINS] = (double)tot;

      /* Put name */

      sprintf(str, "DET_%s", GetText(det + DET_PTR_NAME));

      /* Allocate memory for results for non-adjoint detectors */

      if ((long)RDB[det + DET_PTR_ADJOINT] < VALID_PTR)
	{
	  /* Allocate memory */
	  
	  ptr = NewStat(str, 1, tot);
	
	  /* Put pointer */
	  
	  WDB[det + DET_PTR_STAT] = (double)ptr;

	  /* Allocate memory for history */

	  if((long)RDB[det + DET_WRITE_HIS] == 1)
	    AllocStatHistory(ptr);
	}

      /***********************************************************************/

      /***** Set particle type ***********************************************/

      /* Check if not set */

      if ((long)RDB[det + DET_PARTICLE] == 0)
	{
	  /* Check combined mode */

	  if (((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES) &&
	      ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES))
	    Error(det, "Particle type must be set in combined neutron / photon mode");
	  else if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
	    WDB[det + DET_PARTICLE] = (double)PARTICLE_TYPE_NEUTRON;
	  else if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
	    WDB[det + DET_PARTICLE] = (double)PARTICLE_TYPE_GAMMA;
	  else
	    Die(FUNCTION_NAME, "What the hell am I doing here?");
	}

      /***********************************************************************/

      /* Next detector */

      det = NextItem(det);
    }
  
  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];  
  while (det > VALID_PTR)
    {
      /* Check adjoint type */
      
      if ((ptr = (long)RDB[det + DET_PTR_ADJOINT]) > VALID_PTR)
	{
	  /* Get sizes */

	  n1 = (long)RDB[det + DET_N_TOT_BINS];
	  n2 = (long)RDB[ptr + DET_N_TOT_BINS];
  
	  /* Allocate memory */

	  ptr = NewStat(str, 2, n1, n2);

	  /* Put pointer */
	  
	  WDB[det + DET_PTR_STAT] = (double)ptr;
	}
      
      /* Next detector */

      det = NextItem(det);
    }

  /***************************************************************************/

  /***** Checks **************************************************************/

  /* Loop over detectors */

  det = (long)RDB[DATA_PTR_DET0];  
  while (det > VALID_PTR)
    {
      /***********************************************************************/

      /***** Check multiplier / divider detectors ****************************/

      /* Check that energy bin structure is defined when needed */

      if ((((long)RDB[det + DET_TYPE] == -2) || 
	   ((long)RDB[det + DET_TYPE] == -3)) && 
	  ((long)RDB[det + DET_PTR_EGRID] < VALID_PTR))
	Error(det, "Energy binning must be provided for type %ld detector",
	      (long)RDB[det + DET_TYPE]);

      /* Check multipliers / dividers */
	  
      if ((ptr = (long)RDB[det + DET_PTR_MUL]) > VALID_PTR)
	{
	  /* Check number of values */

	  if ((long)RDB[ptr + DET_N_TOT_BINS] > 1)
	    {
	      /* Check each bin */

	      if (RDB[det + DET_N_EBINS] != RDB[ptr + DET_N_EBINS])
		Error(det, "Mismatch in energy bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      if (RDB[det + DET_N_UBINS] != RDB[ptr + DET_N_UBINS])
		Error(det, "Mismatch in universe bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      if (RDB[det + DET_N_CBINS] != RDB[ptr + DET_N_CBINS])
		Error(det, "Mismatch in cell bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      if (RDB[det + DET_N_MBINS] != RDB[ptr + DET_N_MBINS])
		Error(det, "Mismatch in material bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      if (RDB[det + DET_N_LBINS] != RDB[ptr + DET_N_LBINS])
		Error(det, "Mismatch in lattice bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      if (RDB[det + DET_N_RBINS] != RDB[ptr + DET_N_RBINS])
		Error(det, "Mismatch in reaction bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      if (RDB[det + DET_N_TBINS] != RDB[ptr + DET_N_TBINS])
		Error(det, "Mismatch in time bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      
	      /* Mesh bins */

	      msh1 = (long)RDB[det + DET_PTR_MESH];
	      msh2 = (long)RDB[ptr + DET_PTR_MESH];

	      if (msh1*msh2 < 0)
		Error(det, "Mismatch in mesh bins of detector %s",
		      GetText(ptr + DET_PTR_NAME));
	      else if ((msh1 > VALID_PTR) && (msh2 > VALID_PTR))
		{
		  if (RDB[msh1 + MESH_N0] != RDB[msh2 + MESH_N0])
		    Error(det, "Mismatch in mesh bins of detector %s",
			  GetText(ptr + DET_PTR_NAME));
		  if (RDB[msh1 + MESH_N1] != RDB[msh2 + MESH_N1])
		    Error(det, "Mismatch in mesh bins of detector %s",
			  GetText(ptr + DET_PTR_NAME));
		  if (RDB[msh1 + MESH_N2] != RDB[msh2 + MESH_N2])
		    Error(det, "Mismatch in mesh bins of detector %s",
			  GetText(ptr + DET_PTR_NAME));
		}
	    }
	}
      
      /***********************************************************************/

      /* Next detector */

      det = NextItem(det);
    }

  /***************************************************************************/
}  

/*****************************************************************************/
