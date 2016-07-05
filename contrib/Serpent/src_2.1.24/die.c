/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : die.c                                          */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2012/05/29 (JLe)                                           */
/* Version:       2.1.6                                                      */
/*                                                                           */
/* Description: Terminates run in fatal error                                */
/*                                                                           */
/* Comments: - From Serpent 1.1.8                                            */
/*           - Use this for errors in code, user errors terminate the run    */
/*             with Error()                                                  */
/*           - This and Error() are the only functions that should make      */
/*             a call to exit().                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "Die:"

/*****************************************************************************/

int Die(char *func, ...)
{
  va_list argp;
  va_start (argp, func);

  /* Print error message */
  
  fprintf(err, "\n***** %s (seed = %lu, MPI task = %d, OMP thread = %d)\n\n", 
	  TimeStamp(), parent_seed, mpiid, OMP_THREAD_NUM);

  fprintf(err, "Fatal error in function %s\n\n", func);
  vfprintf(err, va_arg(argp, char *), argp);
  fprintf(err, "\n\n");

  /*
  if ((int)RDB[DATA_RUNNING_MODE] == RUNNING_MODE_XSTEST)
    return 0;
  */

  fprintf(err, "Simulation aborted.\n\n");

  /* Exit with value -1 to terminate all MPI tasks */

  exit(-1);
}

/*****************************************************************************/
