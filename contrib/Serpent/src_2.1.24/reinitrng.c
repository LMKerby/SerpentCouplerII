/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : reinitrng.c                                    */
/*                                                                           */
/* Created:       2011/03/03 (TVi)                                           */
/* Last modified: 2015/04/04 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: This routine is used to skip n*2^STRIDE steps forward in the */
/*              random number sequence beginning from parent seed.           */ 
/*                                                                           */
/* Comments: Algorithm taken from MCNP5. Basically calculates                */
/*           g^(n*2^STRIDE)*parentseed+(g^(n*2^STRIDE)-1)/(g-1)*12345mod2^64 */
/*           efficiently.                                                    */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReInitRNG:"

/*****************************************************************************/

unsigned long ReInitRNG(long n)
{
  unsigned long gen, g, inc, c, gp;
  
  /* Adjust index in MPI mode */
  
  if ((mpiid > 0) && ((long)RDB[DATA_OPTI_MPI_REPRODUCIBILITY] == NO))
    n = n + (long)(mpiid*RDB[DATA_MPI_TOT_PARTICLES]);

  /* Re-initialize RNG */

  n = n<<STRIDE;

  gen = 1;
  g = 2862933555777941757;
  inc = 0;
  c = 12345;

  while(n > 0)
    {
      if((n%2) == 1)
	{
	  gen *= g;
	  inc = inc*g + c;
	}
      
      gp = g + 1;
      g *= g;
      c *= gp;
      
      n = n>>1;
    }  

  return parent_seed*gen + inc;
}

/*****************************************************************************/
