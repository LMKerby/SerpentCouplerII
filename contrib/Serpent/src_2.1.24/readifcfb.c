/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readifcfb.c                                    */
/*                                                                           */
/* Created:       2014/10/06 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Reads fuel behavior multi-physics interfaces                 */
/*                                                                           */
/* Comments:   -Split from readinterface.c for 2.1.22                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadIFCFB:"

/*****************************************************************************/

void ReadIFCFB(long loc0, long update)
{
  long loc1, loc2, loc3, ptr, type, np, n, m;
  long nz, i, nr, nax;
  long nt, na, uni, nnt, tme, o, l, tbi;
  long DFptr, T0ptr, T1ptr, CRptr, HRptr;
  double zmin, zmax;
  double T,  rmin, rmax, r1, r2;
  double tmin, tmax, amin, amax, emin, emax;
  double eps, maxeps, maxr, maxT, maxTp;
  char tmpstr[MAX_STR];
  FILE *fp;

  /* Avoid compiler warning */

  maxTp = 0.0;
  maxT = 0.0;
  maxr = 0.0;
  loc1 = -1;
  loc3 = -1;

  /* Open file for reading */
      
  if ((fp = fopen(GetText(loc0 + IFC_PTR_INPUT_FNAME), "r")) == NULL)
    Error(loc0, "Multi-physics interface file \"%s\" does not exist",
	  GetText(loc0 + IFC_PTR_INPUT_FNAME));

  /* Read interface type */

  if (fscanf(fp, "%ld", &type) == EOF)
    Die(FUNCTION_NAME, "fscanf error");


  /* Interface to fuel performance codes */

  WDB[loc0 + IFC_TYPE] = (double)type;
  WDB[loc0 + IFC_PTR_MAT] = NULLPTR;
  WDB[loc0 + IFC_CALC_OUTPUT] = (double)YES;

  /* Reset output flag to avoid going to next loop */

  n = NO;

  /* Read output file name */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  WDB[loc0 + IFC_PTR_OUTPUT_FNAME] = (double)PutText(tmpstr);

  /*******************************************************************/

  /***** Interface to fuel performance codes *************************/

  /* Reset maximum relative difference */

  maxeps = 0;

  /* Read number of pins */

  if (fscanf(fp, "%ld", &np) == EOF)
    Die(FUNCTION_NAME, "fscanf error");

  /* Loop */

  for (i = 0; i < np; i++)
    {
      /* Allocate memory for structure */
      if(!update)		
	loc1 = NewItem(loc0 + IFC_PTR_FUEP, IFC_FUEP_LIST_BLOCK_SIZE);
      else if(i>0)
	loc1 = NextItem(loc1);
      else
	loc1 = (long)RDB[loc0 + IFC_PTR_FUEP];

      /* Read universe */
	      
      if (fscanf(fp, "%s", tmpstr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      na = atol(tmpstr);

      if(na < 0)
	{
	  na = -na;
	  uni = 0;
	}
      else
	{
	  uni = na;
	  na = 1;
	}
		
      WDB[loc1 + IFC_FUEP_N_UNI] = (double)na;

      /* Allocate memory for pin universes */
      if(!update)	       
	{
	  ptr = ReallocMem(DATA_ARRAY, na);
	  WDB[loc1 + IFC_FUEP_PTR_UNI_LIST] = (double)ptr;
	}
      else
	ptr = (long)RDB[loc1 + IFC_FUEP_PTR_UNI_LIST];
	      
      if(uni != 0)
	{
	  WDB[ptr] = (double)PutText(tmpstr);
	}
      else
	{
	  for(n = 0; n < na; n++)
	    {
	      if(fscanf(fp, "%s", tmpstr) == EOF)
		Die(FUNCTION_NAME, "fscanf error");

	      WDB[ptr++] = (double)PutText(tmpstr);
	    }
	}
		  
      /* Read output meshing */
      if(update == 0)	       
	{
	  ptr = ReallocMem(DATA_ARRAY, FUEP_LIM_SIZE);
		  
	  WDB[loc1 + IFC_FUEP_OUT_PTR_LIM] = (double)ptr;
	}
      else
	{
	  ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];
	}

      /* Read number of axial zones */

      if (fscanf(fp, "%ld", &nz) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[ptr + FUEP_NZ] = (double)nz;

      /* Read axial limits if equally spaced zones */

      if (nz > 0)
	{
	  if(fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
	  WDB[ptr + FUEP_ZMIN] = zmin;

	  CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);
	  WDB[ptr + FUEP_ZMAX] = zmax;
		  
	}
      else if (nz < 0)
	{
	  WDB[ptr + FUEP_ZMIN] = 0.0;
	  WDB[ptr + FUEP_ZMAX] = 0.0;
	}
      else
	Die(FUNCTION_NAME, "Zero axial output zones");
      /* Read number of angular zones */

      if (fscanf(fp, "%ld", &na) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[ptr + FUEP_NA] = (double)na;

      /* Read angular limits if equally spaced zones */

      if (na > 0)
	{
	  if(fscanf(fp, "%lf %lf", &amin, &amax) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  CheckValue(FUNCTION_NAME, "amin", "", amin, -360.0, 360.0);
	  WDB[ptr + FUEP_AMIN] = amin;

	  CheckValue(FUNCTION_NAME, "amax", "", amax, amin, 360.0);
	  WDB[ptr + FUEP_AMAX] = amax;
		  
	}
      else if (na < 0)
	{
	  WDB[ptr + FUEP_AMIN] = 0.0;
	  WDB[ptr + FUEP_AMAX] = 0.0;
	}
      else
	Die(FUNCTION_NAME, "Zero angular output zones");

      /* Read number of radial zones */

      if (fscanf(fp, "%ld", &nr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[ptr + FUEP_NR] = (double)nr;

      /* Read radial limits if equal area zones */

      if (nr > 0)
	{
	  if(fscanf(fp, "%lf %lf", &rmin, &rmax) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  CheckValue(FUNCTION_NAME, "rmin", "", rmin, 0.0, INFTY);
	  WDB[ptr + FUEP_RMIN] = rmin;

	  CheckValue(FUNCTION_NAME, "amax", "", rmax, rmin, INFTY);
	  WDB[ptr + FUEP_RMAX] = rmax;
		  
	}
      else if (nr < 0)
	{
	  WDB[ptr + FUEP_RMIN] = 0.0;
	  WDB[ptr + FUEP_RMAX] = 0.0;
	}
      else
	Die(FUNCTION_NAME, "Zero radial output zones");

      /* TODO: read nnt from somewhere */

      nnt = 1;

      if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
	{
	  WDB[ptr + FUEP_TMIN] = RDB[DATA_DYN_TMIN];
	  WDB[ptr + FUEP_TMAX] = RDB[DATA_DYN_TMAX];
	  WDB[ptr + FUEP_NT] = (double)(RDB[DATA_DYN_NB]*nnt);
	}
      else
	{
	  WDB[ptr + FUEP_TMIN] = -INFTY/10.0;
	  WDB[ptr + FUEP_TMAX] =  INFTY/10.0;
	  WDB[ptr + FUEP_NT] = (double)1;
	}

      WDB[loc1 + IFC_FUEP_TYPE] = (double)type;
	      
      /* Reset pointers */
      if(!update)
	{
	  WDB[loc1 + IFC_FUEP_OUT_PTR_Z] = NULLPTR;
	  WDB[loc1 + IFC_FUEP_OUT_PTR_PHI] = NULLPTR;
	  WDB[loc1 + IFC_FUEP_OUT_PTR_R2] = NULLPTR;
	}

      loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];

      /* Read axial meshing from file */
      if( (nz = (long)RDB[loc2 + FUEP_NZ]) < 0)
	{
	  nz = -1*nz;
	  WDB[loc2 + FUEP_NZ] = (double)nz;

	  ptr = ReadIFCBins(nz, fp, -INFTY/10.0, update);

	  if(!update)
	    WDB[loc1 + IFC_FUEP_OUT_PTR_Z] = (double)ptr;
	  else
	    ptr = RDB[loc1 + IFC_FUEP_OUT_PTR_Z];
		  
	  WDB[loc2 + FUEP_ZMIN] = RDB[ptr +  0];
	  WDB[loc2 + FUEP_ZMAX] = RDB[ptr + nz];
	}		    

      /* Read angular meshing from file */
      if( (na = (long)RDB[loc2 + FUEP_NA]) < 0)
	{
	  na = -1*na;
	  WDB[loc2 + FUEP_NA] = (double)na;

	  ptr = ReadIFCBins(na, fp, -360.0, update);

	  if(!update)
	    WDB[loc1 + IFC_FUEP_OUT_PTR_PHI] = (double)ptr;
	  else
	    ptr = RDB[loc1 + IFC_FUEP_OUT_PTR_PHI];

	  WDB[loc2 + FUEP_AMIN] = RDB[ptr +  0];
	  WDB[loc2 + FUEP_AMAX] = RDB[ptr + na];
	}		    

      /* Read radial meshing from file */
      if( (nr = (long)RDB[loc2 + FUEP_NR]) < 0)
	{
	  nr = -1*nr;
	  WDB[loc2 + FUEP_NR] = (double)nr;

	  ptr = ReadIFCBins(nr, fp, 0.0, update);

	  if(!update)
	    WDB[loc1 + IFC_FUEP_OUT_PTR_R2] = (double)ptr;
	  else
	    ptr = RDB[loc1 + IFC_FUEP_OUT_PTR_R2];

	  WDB[loc2 + FUEP_RMIN] = RDB[ptr +  0];
	  WDB[loc2 + FUEP_RMAX] = RDB[ptr + nr];
	}		    

      /* Read output meshing (flux) */
      if(!update)
	{
	  ptr = ReallocMem(DATA_ARRAY, FUEP_LIM_SIZE);
	  WDB[loc1 + IFC_FUEP_OUT_PTR_FLIM] = (double)ptr;
	}
      else
	ptr = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_FLIM];
      /* TODO: read nnt from somewhere */

      nnt = 1;

      if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
	{
	  WDB[ptr + FUEP_TMIN] = RDB[DATA_DYN_TMIN];
	  WDB[ptr + FUEP_TMAX] = RDB[DATA_DYN_TMAX];
	  WDB[ptr + FUEP_NT] = (double)(RDB[DATA_DYN_NB]*nnt);
	}
      else
	{
	  WDB[ptr + FUEP_TMIN] = -INFTY/10.0;
	  WDB[ptr + FUEP_TMAX] =  INFTY/10.0;
	  WDB[ptr + FUEP_NT] = (double)1;
	}

      /* Read number of axial zones */

      if (fscanf(fp, "%ld", &nz) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[ptr + FUEP_NZ] = (double)nz;

      /* Read axial limits if equally spaced zones */

      if (nz > 0)
	{
	  if(fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);
	  WDB[ptr + FUEP_ZMIN] = zmin;

	  CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);
	  WDB[ptr + FUEP_ZMAX] = zmax;
		  
	}
      else if (nz < 0)
	{
	  WDB[ptr + FUEP_ZMIN] = 0.0;
	  WDB[ptr + FUEP_ZMAX] = 0.0;
	}
      else
	Die(FUNCTION_NAME, "Zero axial output zones");
      /* Read number of angular zones */

      if (fscanf(fp, "%ld", &na) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[ptr + FUEP_NA] = (double)na;

      /* Read angular limits if equally spaced zones */

      if (na > 0)
	{
	  if(fscanf(fp, "%lf %lf", &amin, &amax) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  CheckValue(FUNCTION_NAME, "amin", "", amin, -360.0, 360.0);
	  WDB[ptr + FUEP_AMIN] = amin;

	  CheckValue(FUNCTION_NAME, "amax", "", amax, amin, 360.0);
	  WDB[ptr + FUEP_AMAX] = amax;
		  
	}
      else if (na < 0)
	{
	  WDB[ptr + FUEP_AMIN] = 0.0;
	  WDB[ptr + FUEP_AMAX] = 0.0;
	}
      else
	Die(FUNCTION_NAME, "Zero angular output zones");

      /* Read number of radial zones */

      if (fscanf(fp, "%ld", &nr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[ptr + FUEP_NR] = (double)nr;

      /* Read radial limits if equal area zones */

      if (nr > 0)
	{
	  if(fscanf(fp, "%lf %lf", &rmin, &rmax) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  CheckValue(FUNCTION_NAME, "rmin", "", rmin, 0.0, INFTY);
	  WDB[ptr + FUEP_RMIN] = rmin;

	  CheckValue(FUNCTION_NAME, "amax", "", rmax, rmin, INFTY);
	  WDB[ptr + FUEP_RMAX] = rmax;
		  
	}
      else if (nr < 0)
	{
	  WDB[ptr + FUEP_RMIN] = 0.0;
	  WDB[ptr + FUEP_RMAX] = 0.0;
	}
      else
	Die(FUNCTION_NAME, "Zero radial output zones");

      /* Read fast flux energy limits */

      if (fscanf(fp, "%lf %lf", &emin, &emax) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      WDB[ptr + FUEP_EMIN] = emin;
      WDB[ptr + FUEP_EMAX] = emax;

	      
      /* Reset pointers */
      if(!update)
	{
	  WDB[loc1 + IFC_FUEP_OUT_PTR_FZ] = NULLPTR;
	  WDB[loc1 + IFC_FUEP_OUT_PTR_FPHI] = NULLPTR;
	  WDB[loc1 + IFC_FUEP_OUT_PTR_FR2] = NULLPTR;
	}

      /* Read flux meshings */

      loc2 = (long)RDB[loc1 + IFC_FUEP_OUT_PTR_LIM];

      /* Read axial meshing from file */
      if( (nz = (long)RDB[loc2 + FUEP_NZ]) < 0)
	{
	  nz = -1*nz;
	  WDB[loc2 + FUEP_NZ] = (double)nz;

	  ptr = ReadIFCBins(nz, fp, -INFTY/10.0, update);

	  if(!update)
	    WDB[loc1 + IFC_FUEP_OUT_PTR_FZ] = (double)ptr;
	  else
	    ptr = RDB[loc1 + IFC_FUEP_OUT_PTR_FZ];
		  
	  WDB[loc2 + FUEP_ZMIN] = RDB[ptr +  0];
	  WDB[loc2 + FUEP_ZMAX] = RDB[ptr + nz];
	}		    

      /* Read angular meshing from file */
      if( (na = (long)RDB[loc2 + FUEP_NA]) < 0)
	{
	  na = -1*na;
	  WDB[loc2 + FUEP_NA] = (double)na;

	  ptr = ReadIFCBins(na, fp, -360.0, update);

	  if(!update)
	    WDB[loc1 + IFC_FUEP_OUT_PTR_FPHI] = (double)ptr;
	  else
	    ptr = RDB[loc1 + IFC_FUEP_OUT_PTR_FPHI];

	  WDB[loc2 + FUEP_AMIN] = RDB[ptr +  0];
	  WDB[loc2 + FUEP_AMAX] = RDB[ptr + na];
	}		    

      /* Read radial meshing from file */
      if( (nr = (long)RDB[loc2 + FUEP_NR]) < 0)
	{
	  nr = -1*nr;
	  WDB[loc2 + FUEP_NR] = (double)nr;

	  ptr = ReadIFCBins(nr, fp, 0.0, update);

	  if(!update)
	    WDB[loc1 + IFC_FUEP_OUT_PTR_FR2] = (double)ptr;
	  else
	    ptr = RDB[loc1 + IFC_FUEP_OUT_PTR_FR2];

	  WDB[loc2 + FUEP_RMIN] = RDB[ptr +  0];
	  WDB[loc2 + FUEP_RMAX] = RDB[ptr + nr];
	}		    


      /* Create time intervals */
	      
      /* Put input time limits */
	      
      /* Allocate time arrays */

      nt = (long)RDB[DATA_DYN_NB];
      tmin = RDB[DATA_DYN_TMIN];
      tmax = RDB[DATA_DYN_TMAX];

      if(nt==0)
	nt = 1;
      else
	nt = nt*nnt;

      WDB[loc1 + IFC_FUEP_N_T] = (double)nt;
      WDB[loc1 + IFC_FUEP_TMIN] = tmin;
      WDB[loc1 + IFC_FUEP_TMAX] = tmax;


      if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
	{
	  WDB[loc1 + IFC_FUEP_TMIN] = RDB[DATA_DYN_TMIN];
	  WDB[loc1 + IFC_FUEP_TMAX] = RDB[DATA_DYN_TMAX];
	  WDB[loc1 + IFC_FUEP_N_T] = (double)(RDB[DATA_DYN_NB]*nnt);
	}
      else
	{
	  WDB[loc1 + IFC_FUEP_TMIN] = -INFTY;
	  WDB[loc1 + IFC_FUEP_TMAX] =  INFTY;
	  WDB[loc1 + IFC_FUEP_N_T] = (double)1;
	}

      /* Loop over time intervals */

      nt = (long)RDB[loc1 + IFC_FUEP_N_T];	     

      tme = (long)RDB[DATA_DYN_PTR_TIME_BINS];	     

      tme = (long)RDB[tme + TME_PTR_BINS];

      /* Allocate memory for time bin limit list */

      if (update == 0)
	{

	  /* Allocate memory for time zones */

	  ptr = ReallocMem(DATA_ARRAY, nt + 1);
	  WDB[loc1 + IFC_FUEP_LIM_PTR_T] = (double)ptr;

	}


      /* Reset time bin */

      tbi = 0;

      for(o = 0; o < nt; o++)
	{
		  

	  /* Allocate memory */
	  if(!update)
	    tbi = NewItem(loc1 + IFC_FUEP_PTR_T, 
			  IFC_FUEP_T_BLOCK_SIZE);
	  else if(o > 0)
	    tbi = NextItem(tbi);
	  else
	    tbi = (long)RDB[loc1 + IFC_FUEP_PTR_T];

	  if (((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC) ||
	      ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN))
	    {
	      if(o == 0)
		WDB[tbi + IFC_FUEP_T_TMIN] = 0.0;
	      else
		WDB[tbi + IFC_FUEP_T_TMIN] = RDB[tme + (long)(o/nnt)] + 
		  o%(long)nnt*(RDB[tme + (long)(o/nnt) + 1] - RDB[tme + (long)(o/nnt)])/nnt;

	      WDB[tbi + IFC_FUEP_T_TMAX] = RDB[tme + (long)((o+1)/nnt)] + 
		(o+1)%(long)nnt*(RDB[tme + (long)(o/nnt) + 1] - RDB[tme + (long)(o/nnt)])/nnt;

	      /* Put limits for time zone */

	      ptr = (long)RDB[loc1 + IFC_FUEP_LIM_PTR_T];

	      WDB[ptr + o]     = WDB[tbi + IFC_FUEP_T_TMIN];
	      WDB[ptr + o + 1] = WDB[tbi + IFC_FUEP_T_TMAX];
	      
	    }
	  else
	    {
	      WDB[tbi + IFC_FUEP_T_TMIN] = (double)-INFTY;
	      WDB[tbi + IFC_FUEP_T_TMAX] = (double)INFTY;

	      WDB[ptr]     = WDB[tbi + IFC_FUEP_T_TMIN];
	      WDB[ptr + 1] = WDB[tbi + IFC_FUEP_T_TMAX];

	    }

	  /* Read number of axial zones */

	  if (fscanf(fp, "%ld", &nax) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  /* Put value */

	  WDB[tbi + IFC_FUEP_T_N_AX] = (double)nax;

	  /* Loop over axial zones */
	      
	  for (n = 0; n < nax; n++)
	    {
	      /* Allocate memory */

	      if(!update)
		loc2 = NewItem(tbi + IFC_FUEP_T_PTR_AX, 
			       IFC_FUEP_AX_BLOCK_SIZE);
	      else if (n > 0)
		loc2 = NextItem(loc2);
	      else
		loc2 = (long)RDB[tbi + IFC_FUEP_T_PTR_AX];

	      /* Read minimum and maximum height */
		  
	      if (fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
		Die(FUNCTION_NAME, "fscanf error");

	      /* Put values */

	      WDB[loc2 + IFC_FUEP_AX_ZMIN] = zmin;
	      WDB[loc2 + IFC_FUEP_AX_ZMAX] = zmax;

	      /* Read number of angular zones */

	      if (fscanf(fp, "%ld", &na) == EOF)
		Die(FUNCTION_NAME, "fscanf error");

	      /* Put value */

	      WDB[loc2 + IFC_FUEP_AX_N_ANG] = (double)na;

	      /* Loop over angular zones */

	      for (m = 0; m < na; m++)
		{
		  /* Allocate memory */
		  
		  if(!update)
		    loc3 = NewItem(loc2 + IFC_FUEP_AX_PTR_ANG, 
				   IFC_FUEP_ANG_BLOCK_SIZE);
		  else if(m > 0)
		    loc3 = NextItem(loc3);
		  else
		    loc3 = (long)RDB[loc2 + IFC_FUEP_AX_PTR_ANG];

		  /* Read minimum and maximum angle */

		  if (fscanf(fp, "%lf %lf", &amin, &amax) == EOF)
		    Die(FUNCTION_NAME, "Was expecting angular input limits");

		  /* Put values */
		      
		  WDB[loc3 + IFC_FUEP_ANG_AMIN] = amin*2*PI/360.0;
		  WDB[loc3 + IFC_FUEP_ANG_AMAX] = amax*2*PI/360.0;

		  /* Read number of radial points */

		  if (fscanf(fp, "%ld", &nr) == EOF)
		    Die(FUNCTION_NAME, "Was expecting the number of radial points");
		      
		  /* Put number */

		  WDB[loc3 + IFC_FUEP_ANG_N_RAD] = (double)nr;

		  if (fscanf(fp, "%lf %lf %lf", &r1, &r2, &T) 
		      == EOF)
		    Die(FUNCTION_NAME, "fscanf error");

		  if(r1 != 0.0)
		    nr=nr+1;

		  /* Add point outside cladding */

		  nr = nr + 1;

		  WDB[loc3 + IFC_FUEP_ANG_N_RAD] = nr;

		  /* Allocate arrays for the points */

		  if(!update)
		    {
		      /* Density factor */

		      DFptr = ReallocMem(DATA_ARRAY, nr);
		      WDB[loc3 + IFC_FUEP_ANG_PTR_DF] = (double)DFptr;

		      /* Temperature at the beginning of time step */

		      T0ptr  = ReallocMem(DATA_ARRAY, nr);
		      WDB[loc3 + IFC_FUEP_ANG_PTR_TEMP0] = (double)T0ptr;

		      /* Temperature at the end of time step */	

		      T1ptr  = ReallocMem(DATA_ARRAY, nr);
		      WDB[loc3 + IFC_FUEP_ANG_PTR_TEMP1] = (double)T1ptr;

		      /* Cold radius */

		      CRptr = ReallocMem(DATA_ARRAY, nr);
		      WDB[loc3 + IFC_FUEP_ANG_PTR_COLD_R2] = (double)CRptr;

		      /* Hot radius */

		      HRptr = ReallocMem(DATA_ARRAY, nr);
		      WDB[loc3 + IFC_FUEP_ANG_PTR_HOT_R2] = (double)HRptr;
		    }
		  else
		    {
		      /* Density factor */
		      DFptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_DF];

		      /* Temperature at the beginning of time step */
		      T0ptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_TEMP0];

		      /* Temperature at the end of time step */	
		      T1ptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_TEMP1];

		      /* Cold radius */
		      CRptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_COLD_R2];

		      /* Hot radius */
		      HRptr = (long)RDB[loc3 + IFC_FUEP_ANG_PTR_HOT_R2];

		    }
		  /* T1 values are written in processinterface */

		  if(r1 != 0.0)
		    {
		      if(update == YES)
			{
			  /* Calculate convergence criterion */
			  eps = fabs((RDB[T0ptr]-T)/RDB[T0ptr]);

			  /* Compare to maximum */
			  if(eps > maxeps)
			    {
			      maxeps = eps;
			      maxr = r1;
			      maxT = T;
			      maxTp = RDB[T0ptr];
			    }
			}

		      WDB[T0ptr]  = T;
		      WDB[CRptr] = 0.0;
		      WDB[HRptr] = 0.0;
		      WDB[T0ptr+1]  = T;
		      WDB[CRptr+1] = r1;
		      WDB[HRptr+1] = r2;
		      l=2;
		    }
		  else
		    {

		      if(update == YES)
			{
			  /* Calculate convergence criterion */
			  eps = fabs((RDB[T0ptr]-T)/RDB[T0ptr]);

			  /* Compare to maximum */
			  if(eps > maxeps)
			    {
			      maxeps = eps;
			      maxr = r1;
			      maxT = T;
			      maxTp = RDB[T0ptr];
			    }

			}

		      WDB[T0ptr]  = T;
		      WDB[CRptr] = r1;
		      WDB[HRptr] = r2;
		      l=1;

		    }

		  for(; l < nr - 1; l++)
		    {
		      /* Read coldradius, hotradius and temperature */
			  
		      if (fscanf(fp, "%lf %lf %lf", &r1, &r2, &T) 
			  == EOF)
			Die(FUNCTION_NAME, "fscanf error");

		      if(update == YES)
			{
			  /* Calculate convergence criterion */
			  eps = fabs((RDB[T0ptr + l]-T)/RDB[T0ptr + l]);

			  /* Compare to maximum */
			  if(eps > maxeps)
			    {
			      maxeps = eps;
			      maxr = r1;
			      maxT = T;
			      maxTp = RDB[T0ptr + l];
			    }
			}
			
		      /* Put values */

		      WDB[DFptr + l] = 1.0;
		      WDB[T0ptr + l]  = T;
		      WDB[CRptr + l] = r1;
		      WDB[HRptr + l] = r2;
		    }

		  /* Add static point outside cladding       */
		  /* This allows to do coordinate transforms */
		  /* Even if cladding creeps inwards         */

		  WDB[DFptr + l] = 1.0;
		  WDB[T0ptr + l] = T;
		  WDB[CRptr + l] = r1*2.0;
		  WDB[HRptr + l] = r1*2.0;

		}
	    }

	}

      /* Close list for time bins */
      /* Direct pointer access speeds up search in ifcpoint.c */
      

      if (!update)
	{
	  /* Get pointer to time binning */

	  tbi = (long)RDB[loc1 + IFC_FUEP_PTR_T];
	  CloseList(tbi);

	}
    }

  /* Set TMS on */
  if(!update)
    WDB[DATA_TMS_MODE] = (double)TMS_MODE_CE;

  if(update)
    {
      /*
      printf("Max Eps was %E (%f -> %f) at %E\n", maxeps, maxTp, maxT, maxr);
      printf("%E %%Epsmax\n", maxeps);
      */
    }

  /*******************************************************************/

  /* Close file */

  fclose(fp);

}

/*****************************************************************************/
