#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : freemem.c                                      */
/*                                                                           */
/* Created:       2010/09/15 (JLe)                                           */
/* Last modified: 2015/04/04 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Frees memory allocated to data blocks                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "FreeMem:"

/*****************************************************************************/

void FreeMem()
{
  /* Free FINIX data */

#ifdef FINIX

  FreeFinix();

#endif

  /* Free data arrays */

  if (ACE != NULL)
    Mem(MEM_FREE, ACE);

  if (WDB != NULL)
    Mem(MEM_FREE, WDB);

  if (RES1 != NULL)
    Mem(MEM_FREE, RES1);

  if (RES2 != NULL)
    Mem(MEM_FREE, RES2);

  if (BUF != NULL)
    Mem(MEM_FREE, BUF);

  if (PRIVA != NULL)
    Mem(MEM_FREE, PRIVA);

  if (ASCII != NULL)
    Mem(MEM_FREE, ASCII);

  if (SEED != NULL)
    Mem(MEM_FREE, SEED);

  if (mpiid > 0)
    fclose(out);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
