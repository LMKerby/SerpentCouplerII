/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : testxs.c                                       */
/*                                                                           */
/* Created:       2010/12/26 (JLe)                                           */
/* Last modified: 2012/05/19 (TVi)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Test routine for reading and processing cross sections etc.  */
/*                                                                           */
/* Comments: - Toimii vain yhdellä xsdata directory filellä                  */
/*           - En testannut, toimiiko (TVi 2015-05-19)                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TestXS:"

#define TOFILE


void TestDistributions(long, double, long);

/*****************************************************************************/

void TestXS()
{
#ifdef TESTING
  long nmat, nnuc, nace, n, m, i, ace, ptr, mat, iso, count, N, maxt, nt, opti, loc;
  long seed0, mix, nth, mode, tl, th, tg, max_opti, nucmax, id;
  char tmpstr[MAX_STR], xsdata[MAX_STR], **nuc;
  double E, tw0, tc0, x, y, z, u, v, w, wgt;
  FILE *inpfile;

  fprintf(out, "\nRunning cross section tester routine...\n\n");

  /* Mode 1: cross section and reaction testing without burnup */
  /* Mode 2: like 1, but with OpenMP */
  /* Mode 3: Test of external burnup routines */
  /* Mode 4: Test of internal burnup routines */
  /* Mode 5: Test of processing routines only */

  mode = 1; 

  /* Number of test points */

  if (mode == 5)
    N = 10000;
  else
    N = 100000;

  /* Get maximum number of threads */

#ifdef OPEN_MP

  if (mode == 2)
    maxt = OMP_MAX_THREADS;
  else
    maxt = 1;

#else
  
  maxt = 1;

#endif

  if (mode == 5)
    max_opti = 0;
  else if (mode > 2)
    max_opti = 1;
  else
    max_opti = 3;
  
  /* Avoid warning msg */

  tw0 = 0.0;
  tc0 = 0.0;

  /***************************************************************************/

  /***** Read available nuclide names ****************************************/

  /* Read ace directory file */

  ReadDirectoryFile();

  /* Loop over ace files and count */

  nace = 0;
  ace = (long)RDB[DATA_PTR_ACE0];

  while (ace > VALID_PTR)
    {
      /* Add counter */

      nace++;

      /* Next */

      ace = (long)ACE[ace + ACE_PTR_NEXT];
    }

  /* Allocate memory */

  uc = (char **)Mem(MEM_ALLOC, nace, sizeof(char *));

  for(n = 0; n < nace; n++)
    nuc[n] = (char *)Mem(MEM_ALLOC, 200, sizeof(char));

  /* Read data */

  n = 0;
  ace = (long)RDB[DATA_PTR_ACE0];

  while (ace > VALID_PTR)
    {
      /* Copy name */

      if((ACE[ace + ACE_TYPE] == NUCLIDE_TYPE_SAB) || 
	 (ACE[ace + ACE_TEMP] == 300.0))
	{
	  RDB[DATA_DUMMY] = ACE[ace + ACE_PTR_NAME];
	  sprintf(nuc[n++], "%s", GetText(DATA_DUMMY));
	}

      /* Next */

      ace = (long)ACE[ace + ACE_PTR_NEXT];
    }

  nace = n;

  /* Remember directory file name */

  ptr = (long)RDB[DATA_PTR_ACEDATA_FNAME_LIST];
  sprintf(xsdata, "%s", GetText(ptr));

  /* Intial seed */

  seed0 = parent_seed;

  /***************************************************************************/

  /***** Loop ****************************************************************/
  
  /* Loop over repeats and counts */

  for (count = 0; count < 100000; count++)
    {
      for (opti = 0; opti < max_opti + 1; opti++)
	for (nt = 1; nt < maxt + 1; nt++)
	  {
	    /***** Reset everything *****************************************/
	  
	    /* Set number of threads */
	  
#ifdef OPEN_MP

	    omp_set_num_threads(nt);
	    
#endif

	    /* Free memory */
	    
	    FreeMem();
	    
	    /* Init data array */
	    
	    InitData();

	    /* Reset thermal library counter */

	    nth = 0;
	    
	    /* Init random number generator */

	    parent_seed = seed0 + 100000*count;
	    srand48(parent_seed);
	    
	    /* Allocate memory for directory file list */
	    
	    ptr = ReallocMem(DATA_ARRAY, 2);
	    
	    /* Put pointer */
	    
	    RDB[DATA_PTR_ACEDATA_FNAME_LIST] = (double)ptr;
	    
	    /* Put name */
	    
	    RDB[ptr++] = (double)PutText(xsdata);
	    
	    /* Put null */
	    
	    RDB[ptr] = NULLPTR;
	    
	    /* Put input file name */
	    
	    RDB[DATA_PTR_INPUT_FNAME] = (double)PutText("tester");
	    
	    if (mode > 2)
	      {
		/* Decay library */
		
		n = (long)(5.0*drand48());
		
		if (n == 0)
		  RDB[DATA_PTR_DECAY_FNAME] = 
		    (double)PutText("sss_jef22.dec");
		else if (n == 1)
		  RDB[DATA_PTR_DECAY_FNAME] = 
		    (double)PutText("sss_jeff31.dec");
		else if (n == 2)
		  RDB[DATA_PTR_DECAY_FNAME] = 
		    (double)PutText("sss_jeff311.dec");
		else if (n == 3)
		  RDB[DATA_PTR_DECAY_FNAME] = 
		    (double)PutText("sss_endfb68.dec");
		else if (n == 4)
		  RDB[DATA_PTR_DECAY_FNAME] = 
		    (double)PutText("sss_endfb7.dec");

		/* Fission yield library */
		
		n = (long)(5.0*drand48());
		
		if (n == 0)
		  RDB[DATA_PTR_NFY_FNAME] = 
		    (double)PutText("sss_jef22.nfy");
		else if (n == 1)
		  RDB[DATA_PTR_NFY_FNAME] = 
		    (double)PutText("sss_jeff31.nfy");
		else if (n == 2)
		  RDB[DATA_PTR_NFY_FNAME] = 
		    (double)PutText("sss_jeff311.nfy");
		else if (n == 3)
		  RDB[DATA_PTR_NFY_FNAME] = 
		    (double)PutText("sss_endfb68.nfy");
		else if (n == 4)
		  RDB[DATA_PTR_NFY_FNAME] = 
		    (double)PutText("sss_endfb7.nfy");
	      }

	    /* Grid thinning */
	    
	    if (drand48() < 10.1)
	      RDB[DATA_ERG_TOL] = -1.0;      
	    else
	      RDB[DATA_ERG_TOL] = 1E-4*drand48();

	    /* Delta-tracking (construct majorant) */

	    if (mode == 5)
	      RDB[DATA_OPT_USE_DT] = (double)NO;
	    else
	      RDB[DATA_OPT_USE_DT] = (double)YES;

	    /* Grid minimum and maximum */
	    
	    RDB[DATA_NEUTRON_EMIN] = 
	      -exp(drand48()*(log(1E-6) - log(1E-12)) + log(1E-12));
	    
	    RDB[DATA_NEUTRON_EMAX] = -drand48()*35.0;
	    
	    /* Ures data (TODO: Add modes) */
	    
	    if (drand48() < 0.5)
	      RDB[DATA_USE_URES] = (double)YES;
	    else
	      RDB[DATA_USE_URES] = (double)NO;
	    
	    /* Set optimization mode */
	    
	    RDB[DATA_OPTIMIZATION_MODE] = (double)opti;
	  
	    /* Set running mode */

	    if (mode == 3)
	      RDB[DATA_RUNNING_MODE] = RUNNING_MODE_EXT_BURN;
	    else if (mode == 4)
	      RDB[DATA_RUNNING_MODE] = RUNNING_MODE_INT_BURN;
	    else 
	      RDB[DATA_RUNNING_MODE] = RUNNING_MODE_TRANSPORT;

	    /* Open files for writing */
	  
#ifdef TOFILE
	    
	    out = fopen("tester.log", "w");
	    err = out;
	    wrn = out;
	    
#endif
	    
	    /* Cross section data plotter */
	    /*
	      RDB[DATA_XSPLOT_NE] = drand48()*100 + 10;
	    */
	    /*****************************************************************/
	  
	    /***** Create materials ******************************************/
	  
	    /* Number of materials and nuclides per material */
	    
	    if (mode == 5)
	      {
		nmat = 1;
		nucmax = 200;
	      }
	    else if (mode > 2)
	      {
		nmat = (long)(1000.0*drand48()) + 1;
		nucmax = 5;
	      }
	    else
	      {
		nmat = 5;
		nucmax = 10;
	      }

	    /* Loop over materials */
	    
	    for (n = 0; n < nmat; n++)
	      {
		/* Create material */
		
		mat = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);
		
		/* Set used-flag */

		SetOption(mat + MATERIAL_OPTIONS, OPT_USED);

		/* Put name */
		
		sprintf(tmpstr, "mat%ld", n + 1);
		RDB[mat + MATERIAL_PTR_NAME] = (double)PutText(tmpstr);
		
		/* Density */
		
		RDB[mat + MATERIAL_ADENS] = drand48();
		
		if ((n > 0) && (mode > 2) && (mode != 5))
		  SetOption(mat + MATERIAL_OPTIONS, OPT_BURN_MAT);
		
		/* Reset temperature */
		
		RDB[mat + MATERIAL_TEMP] = -1.0;

		/* Reset thermal flags */
		
		tl = 0;
		th = 0;
		tg = 0;
				
		/* Number of nuclides */

		if ((mode > 2) & (mode != 5))
		  nnuc = (long)(drand48()*nucmax) + 1;
		else
		  nnuc = nucmax;

		/* Loop over composition */
		
		for (m = 0; m < nnuc; m++)
		  {
		    /* Create new nuclide*/
		    
		    iso = NewItem(mat + MATERIAL_PTR_COMP, 
				  COMPOSITION_BLOCK_SIZE);
		    

		    /* Sample nuclide */
		    
		    do
		      {
			i = (long)(drand48()*nace);
			
			/* Test thermal scattering */
			
			if ((nuc[i][0] == 'l') && (tl == 0))
			  {
			    /* Put flag */
			    
			    tl++;
			    
			    /* Put nuclide name */
			    
			    RDB[iso + COMPOSITION_PTR_NUCLIDE] = 
			      (double)PutText("1001.03c");
			    
			    /* Sample fraction */
			    
			    RDB[iso + COMPOSITION_ADENS] = drand48();
			    
			    /* S(a,b) name */
			    
			    sprintf(tmpstr, "lwtr%ld", nth++);
			    
			    /* Add S(a,b) data */
			    
			    ptr = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);
			    RDB[ptr + THERM_PTR_ALIAS] = 
			      (double)PutText(tmpstr);
			    
			    loc = NewItem(ptr + THERM_PTR_SAB, SAB_BLOCK_SIZE);
			    WDB[loc + SAB_PTR_NAME] = (double)PutText(nuc[i]);

			    RDB[ptr + THERM_T] = -1.0;
			    
			    /* Add material entry */
			    
			    ptr = NewItem(mat + MATERIAL_PTR_SAB, 
					  THERM_BLOCK_SIZE);
			    RDB[ptr + THERM_PTR_ALIAS] = 
			      (double)PutText(tmpstr);
			    RDB[ptr + THERM_ZA] = 1001.0;
			  }
			else if ((nuc[i][0] == 'h') && (th == 0))
			  {
			    /* Put flag */
			    
			    th++;
			    
			    /* Put nuclide name */
			    
			    RDB[iso + COMPOSITION_PTR_NUCLIDE] = 
			      (double)PutText("1002.03c");
			    
			    /* Sample fraction */
			    
			    RDB[iso + COMPOSITION_ADENS] = drand48();
			    
			    /* S(a,b) name */
			    
			    sprintf(tmpstr, "hwtr%ld", nth++);
			    
			    /* Add S(a,b) data */
			    
			    ptr = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);
			    RDB[ptr + THERM_PTR_ALIAS] = 
			      (double)PutText(tmpstr);

			    loc = NewItem(ptr + THERM_PTR_SAB, SAB_BLOCK_SIZE);
			    WDB[loc + SAB_PTR_NAME] = (double)PutText(nuc[i]);
			   
			    RDB[ptr + THERM_T] = -1.0;		    
			    
			    /* Add material entry */
			    
			    ptr = NewItem(mat + MATERIAL_PTR_SAB, 
					  THERM_BLOCK_SIZE);
			    RDB[ptr + THERM_PTR_ALIAS] = 
			      (double)PutText(tmpstr);
			    RDB[ptr + THERM_ZA] = 1002.0;
			  }
			else if ((nuc[i][0] == 'g') && (tg == 0))
			  {
			    /* Put flag */
			    
			    tg++;
			    
			    /* Put nuclide name */
			    
			    RDB[iso + COMPOSITION_PTR_NUCLIDE] = 
			      (double)PutText("6000.03c");
			    
			    /* Sample fraction */
			    
			    RDB[iso + COMPOSITION_ADENS] = drand48();
			    
			    /* S(a,b) name */
			    
			    sprintf(tmpstr, "grph%ld", nth++);
			    
			    /* Add S(a,b) data */
			    
			    ptr = NewItem(DATA_PTR_T0, THERM_BLOCK_SIZE);
			    RDB[ptr + THERM_PTR_ALIAS] = 
			      (double)PutText(tmpstr);

			    loc = NewItem(ptr + THERM_PTR_SAB, SAB_BLOCK_SIZE);
			    WDB[loc + SAB_PTR_NAME] = (double)PutText(nuc[i]);

			    RDB[ptr + THERM_T] = -1.0;		    
			
			    /* Add material entry */

			    ptr = NewItem(mat + MATERIAL_PTR_SAB, 
					  THERM_BLOCK_SIZE);
			    RDB[ptr + THERM_PTR_ALIAS] = 
			      (double)PutText(tmpstr);
			    RDB[ptr + THERM_ZA] = 6000.0;
			  }
			else if (nuc[i][strlen(nuc[i]) - 1] == 'c')
			  {
			    /* Put nuclide name */
			    
			    RDB[iso + COMPOSITION_PTR_NUCLIDE] = 
			      (double)PutText(nuc[i]);
			    
			    /* Sample fraction */
			    
			    RDB[iso + COMPOSITION_ADENS] = drand48();
			  }
			else
			  i = -1;
		      }
		    while (i == -1);
		    
		  }
	      }

	    if (mode < 3)
	      {
		/* Create mixture */

		mat = NewItem(DATA_PTR_M0, MATERIAL_BLOCK_SIZE);
		
		/* Set used-flag */

		SetOption(mat + MATERIAL_OPTIONS, OPT_USED);

		/* Put name */
	    
		sprintf(tmpstr, "mixer");
		RDB[mat + MATERIAL_PTR_NAME] = (double)PutText(tmpstr);
	    
		/* Put first material */

		ptr = (long)RDB[DATA_PTR_M0];

		mix = NewItem(mat + MATERIAL_PTR_MIX, MIXTURE_BLOCK_SIZE);
		RDB[mix + MIXTURE_PTR_MAT] = RDB[ptr + MATERIAL_PTR_NAME];
		RDB[mix + MIXTURE_VFRAC] = drand48();

		/* Loop over others */
	    
		while (ptr < mat)
		  {
		    if ((drand48() < 0.5) &&
			!((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT))
		      {
			mix = NewItem(mat + MATERIAL_PTR_MIX, 
				      MIXTURE_BLOCK_SIZE);
			RDB[mix + MIXTURE_PTR_MAT] = 
			  RDB[ptr + MATERIAL_PTR_NAME];
			RDB[mix + MIXTURE_VFRAC] = drand48();
		      }
		    
		    ptr = NextItem(ptr);
		  }
	      }
	    
	    /*****************************************************************/

	    /***** Generate input file ***************************************/

	    inpfile = fopen("tester", "w");
	    
	    /* Loop over materials */
	    
	    mat = (long)RDB[DATA_PTR_M0];
	    while (mat > VALID_PTR)
	      {
		/* Check if material is material or mixture */
		
		if ((iso = (long)RDB[mat + MATERIAL_PTR_COMP]) > VALID_PTR)
		  {
		    fprintf(inpfile, "\nmat %s %E ", 
			    GetText(mat + MATERIAL_PTR_NAME),
			    RDB[mat + MATERIAL_ADENS]);
		    
		    
		    ptr = (long)RDB[mat + MATERIAL_PTR_SAB];
		    while (ptr > VALID_PTR)
		      {
			fprintf(inpfile, "moder %s %ld ",
				GetText(ptr + THERM_PTR_ALIAS), 
				(long)RDB[ptr + THERM_ZA]);
			
			ptr = NextItem(ptr);
		      }
		    
		    if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
		      fprintf(inpfile, "burn 1\n");
		    else
		      fprintf(inpfile, "\n");
		    
		    /* Loop over composition */
		    
		    while(iso > VALID_PTR)
		      {
			if (RDB[iso + COMPOSITION_ADENS] > 100.0)
			  
			  fprintf(inpfile, "%s %ld\n", 
				  GetText(iso + COMPOSITION_PTR_NUCLIDE),
				  (long)RDB[iso + COMPOSITION_ADENS]);
			else
			  fprintf(inpfile, "%s %f\n", 
				  GetText(iso + COMPOSITION_PTR_NUCLIDE),
				  RDB[iso + COMPOSITION_ADENS]);
			
			/* Next isotope */
			
			iso = NextItem(iso);
		      }
		  }
		else if ((mix = (long)RDB[mat + MATERIAL_PTR_MIX]) > VALID_PTR)
		  {
		    fprintf(inpfile, "\nmix %s\n", 
			    GetText(mat + MATERIAL_PTR_NAME));
		    
		    /* Loop over composition */
		    
		    while(mix > VALID_PTR)
		      {
			fprintf(inpfile, "%s %E\n", 
				GetText(mix + MIXTURE_PTR_MAT),
				RDB[mix + MIXTURE_VFRAC]);
			
			/* Next isotope */
			
			mix = NextItem(mix);
		      }
		  }

		
		/* Next material */
		
		mat = NextItem(mat);
	      }

	    fprintf(inpfile, "\n");

	    ptr = (long)RDB[DATA_PTR_T0];
	    while (ptr > VALID_PTR)
	      {
		loc = (long)RDB[ptr + THERM_PTR_SAB];

		fprintf(inpfile, "therm %s %s\n",
			GetText(ptr + THERM_PTR_ALIAS), 
			GetText(loc + SAB_PTR_NAME));

		ptr = NextItem(ptr);
	      }
	    
	    fprintf(inpfile, "\nset acelib \"%s\"\n", xsdata);
	    
	    if ((long)RDB[DATA_PTR_DECAY_FNAME] > VALID_PTR)
	      fprintf(inpfile, "set declib \"%s\"\n",
		      GetText(DATA_PTR_DECAY_FNAME));

	    if ((long)RDB[DATA_PTR_NFY_FNAME] > VALID_PTR)
	      fprintf(inpfile, "set nfylib \"%s\"\n",
		      GetText(DATA_PTR_NFY_FNAME));

	    if ((long)RDB[DATA_PTR_SFY_FNAME] > VALID_PTR)
	      fprintf(inpfile, "set sfylib \"%s\"\n",
		      GetText(DATA_PTR_SFY_FNAME));


	    if (RDB[DATA_ERG_TOL] > 0.0)
	      fprintf(inpfile, "\nset egrid %E %E %E\n", RDB[DATA_ERG_TOL],
		      -RDB[DATA_NEUTRON_EMIN], -RDB[DATA_NEUTRON_EMAX]);
	    
	    fprintf(inpfile, "\nset opti %ld\n", 
		    (long)RDB[DATA_OPTIMIZATION_MODE]);
	    fprintf(inpfile, "set ures %ld\n", 
		    (long)RDB[DATA_USE_URES]);
	    
	    fprintf(inpfile, "set seed %lu\n", parent_seed);
	    
	    fprintf(inpfile, "\nset xsplot 1000 %E %E\n", -RDB[DATA_NEUTRON_EMIN],
		    -RDB[DATA_NEUTRON_EMAX]);
	    
	    if ((long)RDB[DATA_RUNNING_MODE] == RUNNING_MODE_INT_BURN)
	      fprintf(inpfile, "\ndep\n");

	    fclose(inpfile);
	    
	    /*****************************************************************/
	    
	    /***** Process data **********************************************/

	    /* Re-init random number generator */
	    
	    srand48(parent_seed);
	    
	    /* Process input */
	    
	    ProcessInput();
	    
	    fprintf(out, "Testing with %ld samples...\n", N);
	    
	    if (nt == 1)
	      fprintf(out, "\nLoop %ld, optimization mode %ld, mem %1.2f Mb...\n\n", 
		     count + 1, (long)RDB[DATA_OPTIMIZATION_MODE],
		     RDB[DATA_TOTAL_BYTES]/MEGA);
	    
	    /* Reset statistics */

	    /* Start timing */
	    
	    ResetTimer(0);
	    StartTimer(0);

	    mat = (long)RDB[DATA_PTR_M0];
	    while (mat > VALID_PTR)
	      {
		MakeBurnMatrix(mat);
		mat = NextItem(mat);
	      }

#ifdef OPEN_MP 
	    
#pragma omp parallel private(n, m, i, mat, E, x, y, z, u, v, w, wgt, id)
	    {	
	  
#pragma omp for schedule(dynamic)
	    
#else
	      {
#endif

		/* Get Open MP thread id */

		id = OMP_THREAD_NUM;

		for (n = 0; n < N; n++)
		  {
		    /* Sample target material */
		    
		    m = (long)((double)nmat*drand48());
		    
		    mat = (long)RDB[DATA_PTR_M0];
		    for (i = 0; i < m; i++)
		      mat = NextItem(mat);
		    
		    /* Sample energy */

		    E = exp(drand48()*(log(RDB[DATA_NEUTRON_EMAX]) - 
				     log(RDB[DATA_NEUTRON_EMIN])) 
			    + log(RDB[DATA_NEUTRON_EMIN]));
		    
		    /*
		    E = exp(((double)n/((double)N))*(log(RDB[DATA_NEUTRON_EMAX]) - 
						     log(RDB[DATA_NEUTRON_EMIN])) 
			    + log(RDB[DATA_NEUTRON_EMIN]));
		    */

		    /* Reaction list sums */

		    
		    CheckReaListSum(mat, E, NO);
		    

		    /* Sample reaction */
		    /*
		    ptr = SampleReaction(mat, E, 1.0);  		    
		    */

		    if (mode == 5)
		      {
			/* Loop over composition */

			iso = (long)RDB[mat + MATERIAL_PTR_COMP];
			while (iso > VALID_PTR)
			  {
			    /* Pointer to nuclide */
			    
			    ptr = (long)RDB[iso + COMPOSITION_PTR_NUCLIDE];

			    /* Test distributions */

			    TestDistributions(ptr, E);
			    
			    /* Next nuclide */

			    iso = NextItem(iso);
			  }
		      }
		    else
		      {
			/* Set dummy variables */
			/*
			x = 0.0;
			y = 0.0;
			z = 0.0;

			IsotropicDirection(&u, &v, &w);

			wgt = 1.0;

			Collision(mat, &x, &y, &z, &u, &v, &w, &E, &wgt);
			*/
		      }
		    
		    
		  }     
	      }
	      
	  fprintf(out, "OK.\n\n");
	      
	  /* Stop timing */

	  StopTimer(0);

	  /* Get single-CPU times */
	  
	  if (nt == 1)
	    {
	      tw0 = TimerVal(0);
	      tc0 = TimerCPUVal(0);
	    }
	  
#ifdef TOFILE
	  
	  /* Close file */
	  
	  fclose(out);
	  
#endif
	  
	  /* Print timer data */
	  
	  fprintf(out, "Threads: %ld ", nt);      
	  fprintf(out, "CPU time: %6.2f ", TimerCPUVal(0));
	  fprintf(out, "(%1.2f x) ", tc0/TimerCPUVal(0));
	  fprintf(out, "Wall-clock time: %6.2f ", TimerVal(0));
	  fprintf(out, "(%1.2f x)\n", tw0/TimerVal(0));

	  /* Collect results */

	  ptr = (long)RDB[DATA_PTR_SCORE0];
	  while (ptr > VALID_PTR)
	    {
	      /* Collect data */
	      
	      Die(FUNCTION_NAME, 
		  "Tässä oli collectbuf-funktio, mutta ei oo enää, Hä hää!");
	      
	      /* Next */
	      
	      ptr = NextItem(ptr);
	    }
	  
	  /* Reaction rate output */
	  
	  RROutput((long)RDB[DATA_OPTIMIZATION_MODE], nt);
	  
	  /*******************************************************************/
	    }
	  }
 
   /*************************************************************************/
    
  /* Free memory */

  FreeMem();

  /* Free nuclide array */
      
  for(n = 0; n < nace; n++)
    Mem(MEM_FREE, nuc[n]);
  
  Mem(MEM_FREE, nuc);

#endif
    }


/*****************************************************************************/

  void TestDistributions(long nuc, double E, long id)
{
  long rea;
  double Eout, mu;

  /* Loop over reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Sample random value for mu to prevent warnings */

      mu = drand48();

      /* Sample energy distribution */

      SampleENDFLaw(rea, -1, E, &Eout, &mu, id);
      
      /* Sample angular distribution */
      /*
      SampleMu(rea, E);
      */
      /* Next reaction */

      rea = NextItem(rea);
    }
}

/*****************************************************************************/

