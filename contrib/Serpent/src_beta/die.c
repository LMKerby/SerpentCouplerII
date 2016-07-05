/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : die.c                                          */
/*                                                                           */
/* Created:       2010/11/19 (JLe)                                           */
/* Last modified: 2016/03/04 (JLe)                                           */
/* Version:       2.1.26                                                     */
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
  long ptr, n;
  va_list argp;
  va_start (argp, func);

  /* Print error message */
  
  fprintf(err, "\n***** %s:\n\n", TimeStamp());
  
  fprintf(err, " - MPI task         = %d\n", mpiid);
  fprintf(err, " - OpenMP thread    = %d\n", OMP_THREAD_NUM);
  fprintf(err, " - RNG parent seed  = %lu\n", parent_seed);
  if ((SEED0 != NULL) && (SEED0[OMP_THREAD_NUM*RNG_SZ] > 0))
    fprintf(err, " - RNG history seed = %lu\n", 
	    SEED0[OMP_THREAD_NUM*RNG_SZ]);

  /* Check if simulation is running */

  if ((long)RDB[DATA_SIMULATION_COMPLETED] == NO)
    {
      /* Get history index */
      
      n = -1;
      ptr = (long)RDB[DATA_PTR_PRIVA_HIS_IDX];
      if (ptr > VALID_PTR)
	if ((long)RDB[DATA_PRIVA_MEM_READY] == YES)
	  n = GetPrivateData(ptr, OMP_THREAD_NUM);
      
      if (n > -1)
	fprintf(err, " - RNG history idx  = %ld\n\n", n);
      else
	fprintf(err, "\n");
    }
  else
    fprintf(err, "\n");
  
  fprintf(err, "Fatal error in function %s\n\n", func);
  vfprintf(err, va_arg(argp, char *), argp);
  fprintf(err, "\n\n");

  /*
  if ((int)RDB[DATA_RUNNING_MODE] == RUNNING_MODE_XSTEST)
    return 0;
  */

  if ((long)RDB[DATA_TERMINATE_ON_DIE] == YES)
    {
      fprintf(err, "Simulation aborted.\n\n");
      
      /* Exit with value -1 to terminate all MPI tasks */
      
      exit(-1);
    }
  else
    return 0;
}

/*****************************************************************************/
