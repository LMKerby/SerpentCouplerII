/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printinterfaceoutput.c                         */
/*                                                                           */
/* Created:       2012/02/15 (JLe)                                           */
/* Last modified: 2015/05/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Prints output for multi-physics interface                    */
/*                                                                           */
/* Comments:    -2014/03/13 Added angular and time dependence for fuel perf. */
/*               interface                                                   */
/*              -TODO: relaxation for fast flux                              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "PrintInterfaceOutput:"

/*****************************************************************************/

void PrintInterfaceOutput()
{
  long loc0, loc1, loc2, nz, nr, n, ptr, i, j, k, uni, mfile, na, nt, l, idx;
  double zmin, zmax, zl, vol;
  char outfile[MAX_STR], object[MAX_STR], location[MAX_STR];
  FILE *fp;

  /* Check that interfaces are defined */

  if ((long)RDB[DATA_PTR_IFC0] < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* Loop over interfaces */

  loc0 = (long)RDB[DATA_PTR_IFC0];
  while (loc0 > VALID_PTR)
    {
      /* Check if printed */
      
      if ((long)RDB[loc0 + IFC_CALC_OUTPUT] == NO)
	{
	  /* Pointer to next */

	  loc0 = NextItem(loc0);

	  /* Cycle loop */

	  continue;
	}

      /***********************************************************************/

      /***** Parse output file name ******************************************/

      /* Get file name */

      sprintf(outfile, "%s", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

      /* Check for Matlab m-file format */

      mfile = NO;

      if (strlen(outfile) > 2)
	if ((outfile[strlen(outfile) - 1] == 'm') &&
	    (outfile[strlen(outfile) - 2] == '.'))
	  {
	    /* Adjust name */
	    
	    outfile[strlen(outfile) - 2] = '\0';
    
	    /* Set type flag */

	    mfile = YES;
	  }

      /* Check if burnup mode */

      if ((long)RDB[DATA_BURN_PTR_DEP] > VALID_PTR)
	sprintf(&outfile[strlen(outfile)], "%ld",
		(long)RDB[DATA_BURN_STEP]);
      
      /* Check for matlab format */
      
      if (mfile == YES)
	{
	  /* Adjust file name */
	  
	  sprintf(&outfile[strlen(outfile)], ".m");
	}
      
      /***********************************************************************/

      /***** Write data ******************************************************/
	  
      /* Open file for writing (m*/
	  
      if ((fp = fopen(outfile, "w")) == NULL)
	Error(loc0, "Unable to open output file for writing");

      /* Print variable name in Matlab mode */

      if (mfile == YES)
	fprintf(fp, "IFC%ld = [\n", (long)RDB[loc0 + IFC_IDX]);
      
      /* Check type */

      if (((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FUEP) || 
	  ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_FPIP))
	{
	  /***************************************************************/

	  /***** Interface to fuel perfomance codes **********************/

	  /* Get pointer */
	      
	  loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
              
	  /* Loop over pins */

	  while (loc1 > VALID_PTR)
	    {
	      /* Pointer to universe */
		  
	      uni = (long)RDB[loc1 + IFC_FUEP_PTR_UNI];
	      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

	      /* Get dimensions */
 
	      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];

	      nz = (long)RDB[ptr + FUEP_NZ];
	      na = (long)RDB[ptr + FUEP_NA];
	      nr = (long)RDB[ptr + FUEP_NR];	  
	      nt = (long)RDB[ptr + FUEP_NT];
		  
	      /* Loop over axial and radial zones */
	      for (l = 0; l < nt; l++)
		for (i = 0; i < nz; i++)
		  for (j = 0; j < na; j++)
		    for (k = 0; k < nr; k++)
		      {
			/* Print universe name and indexes */

			fprintf(fp, "1 %6s %4ld %4ld %4ld %4ld ",
				GetText(uni + UNIVERSE_PTR_NAME), i + 1, 
				j + 1, k + 1, l+1); 
			
			/* Print z-coordinates */

			ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_Z];
			CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
			fprintf(fp, "%12.5E %12.5E ", RDB[ptr + i], 
				RDB[ptr + i + 1]);

			/* Print angular limits */

			ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_PHI];
			CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
			fprintf(fp, "%12.5E %12.5E ", RDB[ptr + j]*360.0/2/PI, 
				RDB[ptr + j + 1]*360.0/2/PI);

			/* Print radii */

			ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_R2];
			CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
			fprintf(fp, "%12.5E %12.5E ", sqrt(RDB[ptr + k]), 
				sqrt(RDB[ptr + k + 1]));


			/* Print result */

			if(RDB[DATA_RUN_CC] == YES)
			  {
			    /* In coupled calculation print relaxed power */

			    /* Calculate bin index */

			    idx = i + j*nz + k*nz*na + l*nz*na*nr;

			    /* Get pointer to relaxed power */

			    ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER_REL];
			    CheckPointer(FUNCTION_NAME, "(Pptr1)", DATA_ARRAY, ptr);

			    /* Print relaxed power */

			    fprintf(fp, "%12.5E %7.5f\n", RDB[ptr + idx], 
				    -1.0);			
			  }
			else
			  {  
			    /* Otherwise print normal power */

			    ptr = (long)RDB[loc1 + IFC_FUEP_PTR_POWER];
			    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
			    fprintf(fp, "%12.5E %7.5f\n", Mean(ptr, i, j, k, l), 
				    RelErr(ptr, i, j, k, l));
			  }
		      }


	      /* Next pin */
		  
	      loc1 = NextItem(loc1);
	    }

	  /* Get pointer */
	      
	  loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);
              
	  /* Loop over pins */

	  while (loc1 > VALID_PTR)
	    {
	      /* Pointer to universe */
		  
	      uni = (long)RDB[loc1 + IFC_FUEP_PTR_UNI];
	      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

	      /* Fast flux */
	      /* Get dimensions */
 
	      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];

	      nz = (long)RDB[ptr + FUEP_NZ];
	      na = (long)RDB[ptr + FUEP_NA];
	      nr = (long)RDB[ptr + FUEP_NR];	  
	      nt = (long)RDB[ptr + FUEP_NT];
		  
	      /* Loop over axial and radial zones */
	      for (l = 0; l < nt; l++)
		for (i = 0; i < nz; i++)
		  for (j = 0; j < na; j++)
		    for (k = 0; k < nr; k++)
		      {
		      /* Print universe name and indexes */

		      fprintf(fp, "2 %6s %4ld %4ld %4ld %4ld ",
			      GetText(uni + UNIVERSE_PTR_NAME), i + 1, 
			      j + 1, k + 1, l + 1); 
			
		      /* Print z-coordinates */

		      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FZ];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      fprintf(fp, "%12.5E %12.5E ", RDB[ptr + i], 
			      RDB[ptr + i + 1]);

		      /* Print angular limits */

		      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FPHI];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      fprintf(fp, "%12.5E %12.5E ", RDB[ptr + j]*360.0/2/PI, 
			      RDB[ptr + j + 1]*360.0/2/PI);

		      /* Print radii */

		      ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FR2];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      fprintf(fp, "%12.5E %12.5E ", sqrt(RDB[ptr + k]), 
			      sqrt(RDB[ptr + k + 1]));

		      /* Print result */

		      ptr = (long)RDB[loc1 + IFC_FUEP_PTR_FLUX];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

		      fprintf(fp, "%12.5E %7.5f\n", Mean(ptr, i, j, k, l), 
			      RelErr(ptr, i, j, k, l));
		    }

	      /* Next pin */
		  
	      loc1 = NextItem(loc1);
	    }

	  /***************************************************************/
	}	  
      else if ((long)RDB[loc0 + IFC_TYPE] == IFC_TYPE_TET_MESH)
	{	      
	  /***** Unstructured tetrahedral mesh *******************************/

	  /*******************************************************************/

	  /* Get pointer to statistics */

	  if(RDB[DATA_RUN_CC] == YES)
	    {
	      /* In coupled calculation get pointer to the relaxed power */

	      ptr = (long)RDB[loc0 + IFC_PTR_STAT_REL];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	    }
	  else
	    {
	      ptr = (long)RDB[loc0 + IFC_PTR_STAT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    }

	  /* Get pointer to cells */
	  
	  if((loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH_PRNTS]) < VALID_PTR)
	    loc1 = (long)RDB[loc0 + IFC_PTR_TET_MSH];

	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  /* Get pointer to volumes */

	  loc2 = (long)RDB[loc0 + IFC_PTR_STAT_VOL];
	  CheckPointer(FUNCTION_NAME, "(loc2)", DATA_ARRAY, loc2);

	  /* Check Matlab format */
	  
	  if (mfile == NO)
	    {
	      /* Get output file path */

	      sprintf(location, "%s", GetText(loc0 + IFC_PTR_OUTPUT_FNAME));

	      /* Separate object and location */

	      i = strlen(location) - 1;
	      while (i > 0)
		{
		  if (location[i] == '/')
		    break;
		  
		  i--;
		}

	      /* Copy string */

	      if (i > 0)
		sprintf(object, "%s", &location[i + 1]);
	      else
		sprintf(object, "%s", &location[i]);

	      location[i] = '\0';

	      /* Print some header data */
	      
	      fprintf(fp, "FoamFile\n{\n");
	      fprintf(fp, "    version     2.0;\n");
	      fprintf(fp, "    format      ascii;\n");
	      fprintf(fp, "    class       volScalarField;\n");
	      fprintf(fp, "    location    \"%s\";\n", location);
	      fprintf(fp, "    object      %s;\n}\n", object);
	      fprintf(fp, "\ndimensions      [1 -1 -3 0 0 0 0];\n");
	      fprintf(fp, "\ninternalField   nonuniform List<scalar>\n");
	      fprintf(fp, "%ld\n(\n", (long)RDB[loc0 + IFC_STAT_NREG]);
	    }

	  /* Loop over cells */

	  while (loc1 > VALID_PTR)
	    {
	      /* Check stat index */

	      if ((i = (long)RDB[loc1 + IFC_TET_MSH_STAT_IDX]) > -1)
		{
		  /* Get volume */
		    
		  vol = RDB[loc2 + i];

		  /* Check Matlab mode */
		    
		  if (mfile == YES)
		    {
		      if(RDB[DATA_RUN_CC] == NO )
			fprintf(fp, "%6ld %6ld %12.5E %12.5E %7.5f\n",
				(long)RDB[loc1 + IFC_TET_MSH_IDX] + 1,  
				i + 1, vol, Mean(ptr, i), RelErr(ptr, i));

		      else
			fprintf(fp, "%6ld %6ld %12.5E %12.5E %7.5f\n",
				(long)RDB[loc1 + IFC_TET_MSH_IDX] + 1,  
				i + 1, vol, RDB[ptr + i], 0.0);
		    }
		  else if (vol > 0.0)
		    {
		      /* Print in W/m3 */


		      if(RDB[DATA_RUN_CC] == NO )
			fprintf(fp, "%12.5E\n", 1E+6*Mean(ptr, i)/vol);
		      else
			fprintf(fp, "%12.5E\n", 1E+6*RDB[ptr + i]/vol);

		    }
		  else if (Mean(ptr, i) > 0.0)
		    Die(FUNCTION_NAME, "Zero volume");
		  else
		    fprintf(fp, "%12.5E\n", 0.0);
		}
	      
	      /* Next cell */

	      loc1 = NextItem(loc1);
	    }
	  
	  /* Print remaining OpenFOAM stuff */

	  if (mfile == NO)
	    {
	      fprintf(fp, ")\n;\n\n");

	      /* Check pointer to batches */

	      if ((loc1 = (long)RDB[loc0 + IFC_PTR_OF_BATCHES]) > VALID_PTR)
		{
		  fprintf(fp, "boundaryField\n{\n");

		  /* Loop over batches */

		  while (loc1 > VALID_PTR)
		    {
		      fprintf(fp, "    %s\n    {\n", 
			      GetText(loc1 + IFC_OF_BATCH_PTR_NAME));
		      fprintf(fp, "        type            fixedValue;\n");
		      fprintf(fp, "        value           uniform 0;\n");
		      fprintf(fp, "    }\n");

		      /* Pointer to next */

		      loc1 = NextItem(loc1);
		    }

		  fprintf(fp, "}\n\n");
		}
	    }

	  /*******************************************************************/
	}
      else
	{
	  /*******************************************************************/

	  /***** Other types *************************************************/

	  if(RDB[DATA_RUN_CC] == YES)
	    {
	      /* In coupled calculation get pointer to the relaxed power */

	      ptr = (long)RDB[loc0 + IFC_PTR_STAT_REL];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	    }
	  else
	    {
	      /* Otherwise get direct pointer to statistics */

	      ptr = (long)RDB[loc0 + IFC_PTR_STAT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    }
	  
	  /* get parameters */
	  
	  nz = (long)RDB[loc0 + IFC_NZ];
	  zmin = RDB[loc0 + IFC_ZMIN];
	  zmax = RDB[loc0 + IFC_ZMAX];
	  nr = (long)RDB[loc0 + IFC_NR];
	  n = (long)RDB[loc0 + IFC_STAT_NREG];
	  
	  /* Calculate length */
	  
	  zl = (zmax - zmin)/((double)nz);
	  
	  /* Loop over output regions */
	  
	  loc1 = (long)RDB[loc0 + IFC_PTR_OUT];
	  while (loc1 > VALID_PTR)
	    {   
	      /* Pointer to statistics */
	      
	      if ((loc2 = (long)RDB[loc1 + IFC_OUT_PTR_SCORE]) > VALID_PTR)
		{
		  /* Get stat index */
		  
		  i = (long)RDB[loc2 + IFC_SCORE_STAT_IDX];
		  CheckValue(FUNCTION_NAME, "i", "", i, 0, n);
		  
		  /* Loop over axial and radial segments */
		  
		  for (j = 0; j < nz; j++)
		    for (k = 0; k < nr; k++)
		      {

			/* Calculate index for relaxed tally */

			idx = i + j*n + k*n*nz;		      

			/* Print bottom coordinates */
			
			fprintf(fp, "%12.5E %12.5E %12.5E ", 
				RDB[loc1 + IFC_OUT_X0],
				RDB[loc1 + IFC_OUT_Y0], zmin + n*zl);
			
			/* Print top coordinates */
			
			fprintf(fp, "%12.5E %12.5E %12.5E ", 
				RDB[loc1 + IFC_OUT_X0],
				RDB[loc1 + IFC_OUT_Y0], zmin + (n + 1)*zl);
			
			/* Print radii */
			
			fprintf(fp, "%12.5E ", sqrt((double)k/((double)nr))*
				RDB[loc1 + IFC_OUT_R]);
			fprintf(fp, "%12.5E ", 
				sqrt((double)(k + 1.0)/((double)nr))*
				RDB[loc1 + IFC_OUT_R]);
			
			/* Print results */

			if(RDB[DATA_RUN_CC] == (double)NO)
			  fprintf(fp, "%12.5E %7.5f ", Mean(ptr, i, j, k), 
				  RelErr(ptr, i, j, k));
			else
			  fprintf(fp, "%12.5E 0.0 ", RDB[ptr + idx]);

			/* Newline */
			
			fprintf(fp, "\n");
		      }		  
		}
	      
	      /* Next region */
	      
	      loc1 = NextItem(loc1);
	    }
	  
	  /*******************************************************************/
	}
    
      /* Print delimiter in Matlab mode */
      
      if (mfile == YES)
	fprintf(fp, "];\n");

      /* Close file */
      
      fclose(fp);
      
      /* Next interface */

      loc0 = NextItem(loc0);
    }
}

/*****************************************************************************/
