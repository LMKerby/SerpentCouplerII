/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ccsmatrixprint.c                               */
/*                                                                           */
/* Created:       2011/05/02 (MPu)                                           */
/* Last modified: 2012/01/06 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Printtaa matriisin                                           */
/*                                                                           */
/* Comments: - Käyttö debuggaukseen                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "ccsMatrixPrint:"

/*****************************************************************************/

void ccsMatrixPrint(struct ccsMatrix *cmat)
{
  long k, j;

  fprintf(out, "ccs matrix\n\n");

  for (k = 0; k < cmat->m; k++)
    {
      fprintf(out, "column = %ld\n", k + 1); 
    
      for (j = cmat->colptr[k]; j < cmat->colptr[k + 1]; j++)
	fprintf(out, "(%+-12.5e, %+-12.5e), %ld\n", cmat->values[j].re, 
	       cmat->values[j].im, cmat->rowind[j] + 1);
    }
}

/*****************************************************************************/
