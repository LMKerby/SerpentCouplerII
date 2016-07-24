#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processsources.c                               */
/*                                                                           */
/* Created:       2011/03/02 (JLe)                                           */
/* Last modified: 2016/02/01 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Processes source definitions                                 */
/*                                                                           */
/* Comments: - From serpent 1.1.15                                           */
/*                                                                           */
/*           - Source cells and universes are processed in creategeometry.c  */
/*                                                                           */
/*           - Ton radioactive decay sourcen kanssa ei pitäisi sallia muita  */
/*             fotonilähteitä, muuten normeeraus menee pieleen.              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessSources:"

/*****************************************************************************/

void ProcessSources()
{
  long src, surf, mat, lst, loc0, loc1, n, ptr;
  char *name, *str, tmpstr[MAX_STR];
  double tot;
  FILE *fp;

  /***************************************************************************/

  /***** Source for RIA calculation ******************************************/

  if ((long)RDB[DATA_PTR_RIA0] > VALID_PTR)
    {
      /* Allocate memory for source structure (not included in list) */

      src = ReallocMem(DATA_ARRAY, SRC_BLOCK_SIZE);

      /* Put source weight and type */

      WDB[src + SRC_WGT] = 1.0;
      WDB[src + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;

      /* Reset boundaries */
      
      WDB[src + SRC_XMIN] = -INFTY;
      WDB[src + SRC_XMAX] =  INFTY;
      WDB[src + SRC_YMIN] = -INFTY;
      WDB[src + SRC_YMAX] =  INFTY;
      WDB[src + SRC_ZMIN] = -INFTY;
      WDB[src + SRC_ZMAX] =  INFTY;
      
      /* Reset point */
      
      WDB[src + SRC_X0] = -INFTY;
      WDB[src + SRC_Y0] = -INFTY;
      WDB[src + SRC_Z0] = -INFTY;
      
      /* Reset energy */

      WDB[src + SRC_E] = -INFTY;
      
      /* Put file name */

      sprintf(tmpstr, "%s.src", GetText(DATA_PTR_INPUT_FNAME));
      WDB[src + SRC_READ_PTR_FILE] = (double)PutText(tmpstr);

      /* Allocate memory */

      ptr = ReallocMem(DATA_ARRAY, SRC_FILE_BUF_SIZE*SRC_BUF_BLOCK_SIZE);
      
      /* Put pointer and buffer size */
      
      WDB[src + SRC_READ_PTR_BUF] = (double)ptr;
      WDB[src + SRC_READ_BUF_SZ] = (double)SRC_FILE_BUF_SIZE;
      
      /* Reset file and buffer indexes */
      
      WDB[src + SRC_READ_FILE_POS] = 0.0;
      WDB[src + SRC_READ_BUF_IDX] = (double)SRC_FILE_BUF_SIZE;
      
      /* Put type */

      WDB[src + SRC_READ_FILE_TYPE] = (double)SRC_FILE_TYPE_S1_RENORM;

      /* Put pointer */

      WDB[DATA_PTR_RIA_SRC] = (double)src;
    }

  /***************************************************************************/

  /***** Source for calculation with precursors ******************************/

  /* Check if there is a source input */

  if ((loc0 = (long)RDB[DATA_PTR_PREC_DET]) > VALID_PTR)
    if ((long)RDB[loc0 + PRECDET_PTR_IN_FNAME] > VALID_PTR)
      {
	/*******************************************/
	/* Create initial precursor neutron source */
	/*******************************************/

	/* Allocate memory for source structure (not included in list) */

	src = ReallocMem(DATA_ARRAY, SRC_BLOCK_SIZE);

	/* Reset weight */

	WDB[src + SRC_WGT] = 1.0;
	      
	/* Reset boundaries */

	WDB[src + SRC_XMIN] = -INFTY;
	WDB[src + SRC_XMAX] =  INFTY;
	WDB[src + SRC_YMIN] = -INFTY;
	WDB[src + SRC_YMAX] =  INFTY;
	WDB[src + SRC_ZMIN] = -INFTY;
	WDB[src + SRC_ZMAX] =  INFTY;

	/* Reset point */

	WDB[src + SRC_X0] = -INFTY;
	WDB[src + SRC_Y0] = -INFTY;
	WDB[src + SRC_Z0] = -INFTY;

	/* Set default type to neutron */

	WDB[src + SRC_TYPE] = (double)PARTICLE_TYPE_NEUTRON;

	/* Reset energy */

	WDB[src + SRC_E] = -INFTY;

	/* Source name */
	  
	WDB[src + SRC_PTR_NAME] = (double)PutText("DynsrcPrec");

	/* Print precursor filename to tmpstr */

	sprintf(tmpstr, "%s.precpoints", GetText(loc0 + PRECDET_PTR_IN_FNAME));

	/* Set file name */
	      
	WDB[src + SRC_READ_PTR_FILE] = (double)PutText(tmpstr);
	      
	/* Allocate memory */

	ptr = ReallocMem(DATA_ARRAY, SRC_FILE_BUF_SIZE*SRC_BUF_BLOCK_SIZE);
      
	/* Put pointer and buffer size */
      
	WDB[src + SRC_READ_PTR_BUF] = (double)ptr;
	WDB[src + SRC_READ_BUF_SZ] = (double)SRC_FILE_BUF_SIZE;
      
	/* Reset file and buffer indexes */
      
	WDB[src + SRC_READ_FILE_POS] = 0.0;
	WDB[src + SRC_READ_BUF_IDX] = (double)SRC_FILE_BUF_SIZE;
      
	/* Put type */

	WDB[src + SRC_READ_FILE_TYPE] = (double)SRC_FILE_TYPE_SERPENT1;

	/* Set binary input on */

	WDB[src + SRC_READ_BINARY] = (double)YES;

	/* Put pointer */

	WDB[loc0 + PRECDET_PTR_PREC_SRC] = (double)src;
      }

  /***************************************************************************/

  /***** Process sources *****************************************************/

  /* Check that source definition exists in external source mode */
  
  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) &&
      ((long)RDB[DATA_PTR_SRC0] < VALID_PTR))
    Error(0, "External source mode without source definition");

  /* Check that source definitions exist */

  if ((long)RDB[DATA_PTR_SRC0] < VALID_PTR)
    return;
  
  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];  
  while (src > VALID_PTR)
    {
      /* Check type and reaction mode */

      if (((long)RDB[src + SRC_TYPE] == PARTICLE_TYPE_GAMMA) &&
	  ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT))
	Error(src, "Gamma source %s in criticality souce simulation",
	      GetText(src + SRC_PTR_NAME));

      /***********************************************************************/
      
      /***** Link materials to sources ***************************************/
      
      /* Check pointer */
      
      if ((long)RDB[src + SRC_PTR_MAT] > VALID_PTR)
	{
	  /* Find material */

	  mat = RDB[DATA_PTR_M0];
	  if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, 
				 GetText(src + SRC_PTR_MAT))) > VALID_PTR)
	    {
	      /* Set pointer */
		      
	      WDB[src + SRC_PTR_MAT] = (double)mat;

	      /* Set used-flag */

	      SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
	    }
	  else
	    Error(src, "Material %s in source %s not defined", 
		  GetText(src + SRC_PTR_MAT), 
		  GetText(src + SRC_PTR_NAME));
	}
      
      /***********************************************************************/

      /***** Link radioactive materials to sources ***************************/
      
      /* Check pointer */
      
      if ((long)RDB[src + SRC_PTR_RAD_SRC_MAT] > VALID_PTR)
	{
	  /* Check that decay data library is set (NOTE: this does not */
	  /* guaranee that the data is read?) */

	  if ((long)RDB[DATA_PTR_DECDATA_FNAME_LIST] < VALID_PTR)
	    Error(src, "Decay library must be defined for radiation source");
	    
	  /* Check special */

	  if (!strcmp(GetText(src + SRC_PTR_RAD_SRC_MAT), "-1"))
	    WDB[src + SRC_PTR_RAD_SRC_MAT] = -1.0;
	  else
	    {
	      /* Find material */

	      mat = RDB[DATA_PTR_M0];
	      if ((mat = SeekListStr(mat, MATERIAL_PTR_NAME, 
				     GetText(src + SRC_PTR_RAD_SRC_MAT))) 
		  > VALID_PTR)
		{
		  /* Set pointer */
		  
		  WDB[src + SRC_PTR_RAD_SRC_MAT] = (double)mat;

		  /* Set used-flag */

		  SetOption(mat + MATERIAL_OPTIONS, OPT_USED);
		}
	      else
		Error(src, "Material %s in source %s not defined", 
		      GetText(src + SRC_PTR_RAD_SRC_MAT), 
		      GetText(src + SRC_PTR_NAME));
	    }

	  /* Put global pointer for normalization */

	  if ((long)RDB[DATA_NORM_PTR_RAD_SRC_MAT] != 0)
	    Error(src, "Multiple radioactive decay sources not allowed");
	  else
	    WDB[DATA_NORM_PTR_RAD_SRC_MAT] = RDB[src + SRC_PTR_RAD_SRC_MAT];  
	}
      
      /***********************************************************************/
      
      /***** Link surfaces to sources ****************************************/

      /* Check if surface is given */

      if ((long)RDB[src + SRC_PTR_SURF] > VALID_PTR)
	{
	  /* Get name */
	  
	  str = GetText(src + SRC_PTR_SURF);
	  
	  if (*str == '-')
	    name = &str[1];
	  else
	    name = str;
	  
	  /* Find surface */

	  surf = RDB[DATA_PTR_S0];
	  if ((surf = SeekListStr(surf, SURFACE_PTR_NAME, name)) > VALID_PTR)
	    {
	      /* Set pointer */
		      
	      WDB[src + SRC_PTR_SURF] = (double)surf;

	      /* Set used-flag */

	      SetOption(surf + SURFACE_OPTIONS, OPT_USED);
	    }
	  else
	    Error(src, "Surface %s in source %s not defined", name,
		  GetText(src + SRC_PTR_NAME));

	  /* Set direction */
	  
	  if (*str == '-')
	    WDB[src + SRC_SURF_SIDE] = -1.0;
	  else
	    WDB[src + SRC_SURF_SIDE] = 1.0;
	}

      /***********************************************************************/

      /***** Source files ****************************************************/

      if ((long)RDB[src + SRC_READ_PTR_FILE] > VALID_PTR)
	{
	  /* Open file */

	  if ((fp = fopen(GetText(src + SRC_READ_PTR_FILE), "r")) == NULL)
	    Error(src, "Unable to open source file \"%s\"", 
		  GetText(src + SRC_READ_PTR_FILE));
	  else
	    fclose(fp);

	  /* Test file format */

	  if (!((long)RDB[src + SRC_READ_BINARY]))
	    TestDOSFile(GetText(src + SRC_READ_PTR_FILE));

	  /* Allocate memory */

	  ptr = ReallocMem(DATA_ARRAY, SRC_FILE_BUF_SIZE*SRC_BUF_BLOCK_SIZE);

	  /* Put pointer and buffer size */

	  WDB[src + SRC_READ_PTR_BUF] = (double)ptr;
	  WDB[src + SRC_READ_BUF_SZ] = (double)SRC_FILE_BUF_SIZE;

	  /* Reset file and buffer indexes */

	  WDB[src + SRC_READ_FILE_POS] = 0.0;
	  WDB[src + SRC_READ_BUF_IDX] = (double)SRC_FILE_BUF_SIZE;
	}

      /***********************************************************************/

      /***** Fusion plasma sources *******************************************/

      /* Check type */

      if ((long)RDB[src + SRC_READ_FILE_TYPE] == SRC_FILE_TYPE_FUSION_PLASMA)
	{
	  /* Read and process plasma source */

	  ReadPlasmaSrc(src);
	}

      /***********************************************************************/

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Renormalize directions **********************************************/

  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Calculate total */

      tot = sqrt(RDB[src + SRC_U0]*RDB[src + SRC_U0] + 
		 RDB[src + SRC_V0]*RDB[src + SRC_V0] + 
		 RDB[src + SRC_W0]*RDB[src + SRC_W0]);

      /* Normalize vectors */

      if (tot > 0.0)
	{
	  WDB[src + SRC_U0] = RDB[src + SRC_U0]/tot;
	  WDB[src + SRC_V0] = RDB[src + SRC_V0]/tot;
	  WDB[src + SRC_W0] = RDB[src + SRC_W0]/tot;
	}
      else
	{
	  /* Set x-component to infinity to indicate isotropic */

	  WDB[src + SRC_U0] = INFTY;
	}

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Renormalize weights *************************************************/

  /* Reset total */

  tot = 0.0;

  /* Loop over sources to calculate total */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Add to sum */

      tot = tot + RDB[src + SRC_WGT];

      /* Next source */

      src = NextItem(src);
    }

  /* Check sum */

  if (tot == 0.0)
    Error(0, "Sum of source weights is zero");

  /* Loop over sources to normalize */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Divide by total */

      WDB[src + SRC_WGT] = RDB[src + SRC_WGT]/tot;

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/

  /***** Check energy bin order and normalize ********************************/

  /* Loop over sources */

  src = (long)RDB[DATA_PTR_SRC0];
  while (src > VALID_PTR)
    {
      /* Pointer to energies bins */

      if ((lst = (long)RDB[src + SRC_PTR_EBINS]) > VALID_PTR)
	{
	  /* Loop over list */

	  for (n = 0; n < ListSize(lst) - 1; n++)
	    {
	      /* Get pointers */

	      loc0 = ListPtr(lst, n);
	      loc1 = ListPtr(lst, n + 1);

	      /* Compare */

	      if (RDB[loc0 + SRC_EBIN_EMAX] > RDB[loc1 + SRC_EBIN_EMAX])
		Error(src, "Energy bins are not in ascending order");
	    }
	}

      /* Next source */

      src = NextItem(src);
    }

  /***************************************************************************/
}  

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
