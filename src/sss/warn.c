#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : warn.c                                         */
/*                                                                           */
/* Created:       2010/09/14 (JLe)                                           */
/* Last modified: 2015/10/30 (JLe)                                           */
/* Version:       2.1.25                                                     */
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
  va_list argp;
  va_start (argp, func);

  /* Print warning message */

  if (mpitasks > 1)
    fprintf(err, "\n***** %s (seed = %lu, MPI task = %d)\n", TimeStamp(), 
	    parent_seed, mpiid);
  else
    fprintf(err, "\n***** %s (seed = %lu)\n", TimeStamp(), parent_seed);
  
  fprintf(err, "Warning message from function %s\n\n", func);
  vfprintf(err, va_arg(argp, char *), argp);
  fprintf(err, "\n\n");
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
