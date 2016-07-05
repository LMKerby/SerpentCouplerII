/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : warn.c                                         */
/*                                                                           */
/* Created:       2010/09/14 (JLe)                                           */
/* Last modified: 2014/01/23 (JLe)                                           */
/* Version:       2.1.17                                                     */
/*                                                                           */
/* Description: Prints warning message                                       */
/*                                                                           */
/* Comments: - Tähän yhdeksi argumentiksi pointteri laskuriin?               */
/*           - Noi messaget vois kerätä myös erilliseen fileen               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Warn:"

/*****************************************************************************/

void Warn(char *func, ...)
{
  /* char fname[MAX_STR]; */
  FILE *fp;
  va_list argp;
  va_start (argp, func);

  /* Open file for writing (tota input filea ei välttämättä oo) */

  /*
    sprintf(fname, "%s.wrn", GetText(DATA_PTR_INPUT_FNAME));

  if ((fp = fopen(fname, "a")) == NULL)
    Die(FUNCTION_NAME, "Unable to open file for writing");
  */

  fp = stdout;

  /* Print warning message */

  if (mpitasks > 1)
    fprintf(fp, "***** %s (seed = %lu, MPI task = %d)\n", TimeStamp(), 
	    parent_seed, mpiid);
  else
    fprintf(fp, "***** %s (seed = %lu)\n", TimeStamp(), parent_seed);
  
  fprintf(fp, "Warning message from function %s\n\n", func);
  vfprintf(fp, va_arg(argp, char *), argp);
  fprintf(fp, "\n\n");

  /* Close file */
  /*
  fclose(fp);
  */
}

/*****************************************************************************/
