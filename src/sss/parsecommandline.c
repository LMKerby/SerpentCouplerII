/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : parsecommandline.c                             */
/*                                                                           */
/* Created:       2010/11/21 (JLe)                                           */
/* Last modified: 2016/03/04 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Handles command line input                                   */
/*                                                                           */
/* Comments: - From Serpent 1.1.14                                           */
/*           - Volume checker rutiinin käsittely puuttuu                     */
/*           - Täältä kutsutaan exit(-1)                                     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ParseCommandLine:"

/*****************************************************************************/

long ParseCommandLine(int argc, char **argv)
{
  long n, mode, np, ptr;
  char str[MAX_STR], fname[MAX_STR];
  FILE *fp;

#ifdef MPI

  char tmpstr[MAX_STR];
  tmpstr[0] = '\0';

#endif

  /* Avoid warning messages */

  np = 0;

  /***** Check number of arguments *******************************************/

  if ((argc < 2))
    {
      fprintf(out, "\nUsage: %s <inputfile> [options]\n\n", argv[0]);
      fprintf(out, "       Where <inputfile> is the file path for the main input ");
      fprintf(out, "file and\n       the available options are:\n\n");

      fprintf(out, "       -version           :  print version information and ");
      fprintf(out, "exit\n");
      fprintf(out, "       -replay            :  run simulation using random ");
      fprintf(out, "number seed from a\n                             previous run\n");
      fprintf(out, "       -his               :  run only burnup history in ");
      fprintf(out, "coefficient calculation\n");
      fprintf(out, "       -coe               :  run only restarts in ");
      fprintf(out, "coefficient calculation\n");
      fprintf(out, "       -plot              :  stop after geometry plot\n");
      fprintf(out, "       -checkvolumes <N>  :  calculate Monte Carlo estimates ");
      fprintf(out, "for material\n                             volumes\n");
      fprintf(out, "       -checkstl <N> <M>  :  check for holes and errors in ");
      fprintf(out, "STL geometries\n                             by sampling ");
      fprintf(out, "<M> directions in <N> points\n");
      fprintf(out, "       -mpi <N>           :  run simulation in MPI mode using ");
      fprintf(out, "<N> parallel\n                             tasks\n");
      fprintf(out, "       -omp <M>           :  run simulation in OpenMP mode ");
      fprintf(out, "using <M> parallel\n                             threads\n");
      fprintf(out, "       -disperse          :  generate random particle or ");
      fprintf(out, "pebble distribution\n                             files for ");
      fprintf(out, "HTGR calculations\n");
      fprintf(out, "       -rdep              :  read binary depletion file from ");
      fprintf(out, "previous\n                             calculation and print ");
      fprintf(out, "new output according to\n");
      fprintf(out, "                             inventory list\n");
      fprintf(out, "       -tracks <N>        :  draw particle tracks in ");
      fprintf(out, "the geometry plots\n");
      fprintf(out, "       -comp <mat>        :  print standard material ");
      fprintf(out, "composition that can be\n");
      fprintf(out, "                             copy-pasted into the input");
      fprintf(out, "file\n");
      fprintf(out, "       -elem <mat> <dens> :  decomposes elemental ");
      fprintf(out, "composition into isotopic\n");
      fprintf(out, "       -qp                :  quick plot mode (ignore overlaps)\n");
      fprintf(out, "       -nofatal           :  ignore fatal errors\n");
      fprintf(out, "\n");
      fprintf(out, "       The input file can be omitted for version and disperse ");
      fprintf(out, "options.\n\n");

      exit(-1);
    }

  /***************************************************************************/

  /***** Check version request ***********************************************/

  if (!strcasecmp(argv[1], "-version"))
    {
      /* Print title and exit */

      PrintTitle();

      exit(-1);
    }

  /***************************************************************************/

  /***** Print standard composition ******************************************/

  if (!strcasecmp(argv[1], "-comp"))
    {
      /* Print composition and exit */

      if (argc > 3)
	StdComp(argv[2], argv[3]);
      else if (argc > 2)
	StdComp(argv[2], NULL);
      else
	fprintf(out, "\nMissing composition name\n\n");

      exit(-1);
    }

  /***************************************************************************/

  /***** Decompose elements **************************************************/

  if (!strcasecmp(argv[1], "-elem"))
    {
      /* Print composition and exit */

      if (argc > 4)
	Element(argv[2], argv[3], argv[4]);
      else if (argc > 3)
	Element(argv[2], argv[3], NULL);
      else if (argc > 2)
	fprintf(out, "\nMissing density\n\n");
      else
	fprintf(out, "\nMissing element name\n\n");

      exit(-1);
    }

  /***************************************************************************/

  /***** Check xs test mode **************************************************/

  if (!strcasecmp(argv[1], "-testxs"))
    {
      /* Set file name */

      if (argc > 2)
	{
	  /* Allocate memory for list */

	  ptr = ReallocMem(DATA_ARRAY, 2);

	  /* Put pointer */

	  WDB[DATA_PTR_ACEDATA_FNAME_LIST] = (double)ptr;

	  /* Put name */

	  WDB[ptr++] = (double)PutText(argv[2]);

	  /* Put null */

	  WDB[ptr] = NULLPTR;

	  /* Read file name */

	  if (argc > 3)
	    WDB[DATA_PTR_XSTEST_FNAME] = (double)PutText(argv[3]);
	}

      /* Test cross sections */

      TestXS();

      /* Exit subroutine */

      return -1;
    }

  /***************************************************************************/

  /***** Treat rest of the arguments *****************************************/

  /* Reset running mode and file pointer */

  mode = 0;
  fp = NULL;

  /* Loop over arguments */

  for (n = 1; n < argc; n++)
    {
      /* Test options */

      if (!strcmp(argv[n], "-replay"))
	mode |= MODE_REPLAY;
      else if (!strcmp(argv[n], "-his"))
	WDB[DATA_COEF_CALC_SPECIAL_MODE] = (double)SPECIAL_COEF_MODE_HIS_ONLY;
      else if (!strcmp(argv[n], "-nofatal"))
	WDB[DATA_TERMINATE_ON_DIE] = (double)NO;
      else if (!strcmp(argv[n], "-coe"))
	WDB[DATA_COEF_CALC_SPECIAL_MODE] = (double)SPECIAL_COEF_MODE_COE_ONLY;
      else if (!strcmp(argv[n], "-plot"))
	WDB[DATA_STOP_AFTER_PLOT] = (double)STOP_AFTER_PLOT_GEOM;
      else if (!strcmp(argv[n], "-tracks"))
	{
	  /* Track plotter */

	  if (argc < n + 2)
	    {
	      fprintf(err, "\nNumber of histories not given.\n\n");
	      exit(-1);
	    }
	  else if ((np = atol(argv[n + 1])) < 1)
	    {
	      fprintf(err, "\nInvalid number of histories for track plotter.\n\n");
	      exit(-1);
	    }
	  else
	    {
	      WDB[DATA_STOP_AT_BOUNDARY] = (double)YES;
	      WDB[DATA_STOP_AFTER_PLOT] = (double)STOP_AFTER_PLOT_TRACKS;
	      WDB[DATA_TRACK_PLOTTER_HIS] = (double)np;
	    }

	  n++;
	}
      else if (!strcmp(argv[n], "-testgeom"))
	{
	  fprintf(err, "Geometry tester is not available in Serpent 2.\n");
	  exit(-1);
	}
      else if (!strcasecmp(argv[n], "-rdep"))
	WDB[DATA_PARTICLE_REDEPLETE_MODE] = YES;
      else if (!strcmp(argv[n], "-disperse"))
	WDB[DATA_PARTICLE_DISPERSER_MODE] = YES;
      else if (!strcmp(argv[n], "-qp"))
	WDB[DATA_QUICK_PLOT_MODE] = YES;
      else if ((!strcmp(argv[n], "-checkvolume")) ||
	       (!strcmp(argv[n], "-checkvolumes")))
	{
	  /* Material volume test. */

	  if (argc < n + 2)
	    {
	      fprintf(err, "Number of test points not given.\n");
	      exit(-1);
	    }
	  else
	    {
	      /* Get number of random points */

	      WDB[DATA_VOLUME_MC_NMAX] = atof(argv[++n]);

	      /* Set running mode */

	      WDB[DATA_VOLUME_CALCULATION_MODE] = YES;
	    }
	}
      else if (!strcmp(argv[n], "-checkstl"))
	{
	  /* STL geometry test */

	  if (argc < n + 3)
	    {
	      fprintf(err, "Number of test points not given.\n");
	      exit(-1);
	    }
	  else
	    {
	      /* Get number of random points and directions */

	      WDB[DATA_STL_TEST_N_PTS] = atof(argv[++n]);
	      WDB[DATA_STL_TEST_N_DIR] = atof(argv[++n]);
	    }
	}

#ifdef MPI

      else if (!strcmp(argv[n], "-mpi"))
	{
	  mode |= MODE_MPI;

	  /* Get number of tasks */

	  if (argc < n + 2)
	    {
	      fprintf(err, "Number of MPI tasks not given.\n");
	      exit(-1);
	    }
	  else if ((np = atoi(argv[n + 1])) < 2)
	    {
	      fprintf(err, "Invalid MPI task number \"%s\".\n",argv[n + 1]);
	      exit(-1);
	    }
	  n++;
	}

#else

      else if (!strcmp(argv[n], "-mpi"))
	{
	  fprintf(err, "Parallel calculation not available in compilation.\n");
	  exit(-1);
	}

#endif

#ifdef OPEN_MP

      else if (!strcmp(argv[n], "-omp"))
	{
	  /* Get number of threads */

	  if (argc < n + 2)
	    {
	      fprintf(err, "Number of OpenMP threads not given.\n");
	      exit(-1);
	    }
	  else if (!strcmp(argv[n + 1], "max"))
	    {
	      /* Use maximum available number of threads */

	      WDB[DATA_OMP_MAX_THREADS] = -1.0;
	    }
	  else if (atoi(argv[n + 1]) < 1)
	    {
	      fprintf(err, "Invalid OpenMP thread number \"%s\".\n",
		      argv[n + 1]);
	      exit(-1);
	    }
	  else if (atoi(argv[n + 1]) > MAX_OMP_THREADS)
	    {
	      fprintf(err, "Maximum number of OpenMP threads is %d",
		      MAX_OMP_THREADS);
	      exit(-1);
	    }
	  else
	    {
	      /* Set number of threads */

	      WDB[DATA_OMP_MAX_THREADS] = (double)atoi(argv[n + 1]);
	    }

	  n++;
	}

#else

      else if (!strcmp(argv[n], "-omp"))
	{
	  fprintf(err, "OpenMP mode not available in compilation.\n");
	  exit(-1);
	}

#endif
/* LMK added to couple to MOOSE 7/2016 */
/* Called by external program */
      else if (!strcmp(argv[n], "-ext"))
  {
     WDB[DATA_EXT_MODE] = (double)YES;
     n++;
  }
/* end LMK */
      else if (fp == NULL)
	{

	  /* Argument is a file name, exit if disperse mode */

	  if ((long)RDB[DATA_PARTICLE_DISPERSER_MODE] == YES)
	    {
	      WDB[DATA_PTR_INPUT_FNAME] = PutText(argv[n]);
	      return 0;
	    }

	  /* Check that file exists */

	  else if ((fp = fopen(argv[n], "r")) == NULL)
	    {
	      fprintf(err, "Input file \"%s\" does not exist.\n", argv[n]);
	      exit(-1);
	    }

	  /* Remember file name */

	  strcpy(fname, argv[n]);
	}
      else
	{
	  /* One file name already read */

	  fprintf(err, "Invalid parameter \"%s\".\n", argv[n]);
	  exit(-1);
	}
    }

  /***************************************************************************/

  /***** Check that input file is given **************************************/

  /* Exit if disperse mode */

  if ((long)RDB[DATA_PARTICLE_DISPERSER_MODE] == YES)
    return 0;

  if (fp == NULL)
    {
      fprintf(err, "No input file given.\n");
      exit(-1);
    }
  else
    fclose(fp);

  /***************************************************************************/

  /***** Check for MPI-mode **************************************************/

#ifdef MPI

  if (mode & MODE_MPI)
    {
      /* Parse command string */

      sprintf(str, "%s -np %ld %s %s ", MPIRUN_PATH, np, argv[0], fname);

      /* Add some special options */

      if (mode & MODE_REPLAY)
	strcat(str, "-replay ");

      if ((long)RDB[DATA_OMP_MAX_THREADS] == -1)
	strcat(str, "-omp max ");
      else if ((long)RDB[DATA_OMP_MAX_THREADS] > 1)
	{
	  sprintf(tmpstr, "-omp %ld ", (long)RDB[DATA_OMP_MAX_THREADS]);
	  strcat(str, tmpstr);
	}

      /* Call recursively */

      if (system(str) == -1)
	Die(FUNCTION_NAME, "Recursive call failed");

      /* Free memory */

      FreeMem();

      /* Finalize MPI */

      FinalizeMPI();

      /* Exit */

      return -1;
    }

#endif

  /***************************************************************************/

  /***** Check that seed file exists if in replay mode ***********************/

  if (mode & MODE_REPLAY)
    {
      sprintf(str, "%s.seed", fname);

      if ((fp = fopen(str, "r")) != NULL)
	{
	  if (fscanf(fp, "%lu", &parent_seed) == EOF)
	    Die(FUNCTION_NAME, "fscanf error");

	  fclose(fp);
	}
      else
	{
	  fprintf(err, "Seed file \"%s\" does not exist.\n", str);
	  exit(-1);
	}

      /* Set replay option */

      WDB[DATA_OPTI_REPLAY] = (double)YES;
    }
  else
    {
      /* Write seed file */

      sprintf(str, "%s.seed", fname);

      if ((fp = fopen(str, "w")) == NULL)
	{
	  fprintf(err, "%s Unable to open seed file for writing.\n\n",
		  FUNCTION_NAME);
	  exit(-1);
	}

      fprintf(fp, "%lu\n", parent_seed);

      fclose(fp);
    }

  /***************************************************************************/

  /* Set file name */

  WDB[DATA_PTR_INPUT_FNAME] = PutText(fname);

  /* Init random number generator */

  srand48(parent_seed);

  /* Exit subroutine */

  return 0;
}

/*****************************************************************************/
