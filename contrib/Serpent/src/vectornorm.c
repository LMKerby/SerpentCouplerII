/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : vectorNorm.c                                   */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2014/06/07 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/*                                                                           */
/* Description: Compute 2-norm of complex-valued vector                      */
/*                                                                           */
/* Comments: -                                                               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "vectorNorm:"

/*****************************************************************************/
/*--------------------------------------------------------------------*/
double vectorNorm(long n, complex *v){

  long i; 
  double nrm;  

  nrm = 0.0; 

  for (i=0; i<n; i++){
   
    CheckValue(FUNCTION_NAME, "v[i].re","", v[i].re, -INFTY, INFTY);
    CheckValue(FUNCTION_NAME, "v[i].re","", v[i].im, -INFTY, INFTY);

    nrm = nrm + c_norm( v[i] ) * c_norm( v[i] ); 

  }

  nrm = sqrt(nrm); 
  CheckValue(FUNCTION_NAME, "nrm","", nrm, 0.0, INFTY);

  return nrm; 
}
/*--------------------------------------------------------------------*/
