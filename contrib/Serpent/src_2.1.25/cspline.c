/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : cspline.c                                      */
/*                                                                           */
/* Created:       2014/11/11 (TKa)                                           */
/* Last modified: 2014/11/25 (TKa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Functions for cubic spline interpolation                     */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "CSpline:"

/* TODO: function descriptions */

/* Local function definitions */
static double CSplineF(double **coeff, double x, long idx);



/*****************************************************************************/


/*****************************************************************************/
void CSplineConstruct(const double *x0, const double *f0, long N0, double d2f0,
                     double d2fend, double **coeff, long F0flag, double *F0) {

  /* Construct cubic splines (Uses Thomas algorithm for solving tridiagonal
   *  system) */
  /* TODO: better description,  A d2f = D */
  long i;
  double *d2f, *h, *d, *H, *D, *u, *v;
  double Hhu, Fa, Fb, Ftot;

  /* Check the array size */
  if (N0 < 3)
    Die(FUNCTION_NAME, "Number of elements too small: %ld", N0);

  /* Check the array is sorted */
  for (i = 1; i < N0; i++)
    if (x0[i-1] > x0[i])
      Die(FUNCTION_NAME, "Input array not sorted");

  /* Allocate memory */
  d2f = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
  h = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
  d = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
  H = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
  D = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
  u = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
  v = (double *)Mem(MEM_ALLOC, N0, sizeof(double));

  /* Calculate the elements of the tridiagonal matrix (h and H) and vector d */
  for (i = 0; i < N0-1; i++) {
    h[i] = x0[i+1] - x0[i];
    d[i] = (f0[i+1] - f0[i])/h[i];
  }

  for (i = 1; i < N0-1; i++)
    H[i] = 2.0*(h[i-1] + h[i]);

  /* Set the second derivative of the first and last points */
  d2f[0] = d2f0;
  d2f[N0-1] = d2fend;

  /* More stuff */
  D[1] = 6.0*(d[1] - d[0]) - h[0]*d2f[0];
  D[N0-2] = 6.0*(d[N0-2] - d[N0-3]) - h[N0-2]*d2f[N0-1];

  for (i = 2; i < N0-2; i++)
    D[i] = 6.0*(d[i] - d[i-1]);

  /* Forward substitution */
  u[1] = h[1]/H[1];
  v[1] = D[1]/H[1];
  for (i = 2; i < N0-1; i++) {
    if ((Hhu = H[i] - h[i-1]*u[i-1]) == 0) {
      /* This should not happen! */
      Die(FUNCTION_NAME, "H[i] - h[i-1]*u[i-1] is zero, can not construct spline");
    }
    u[i] = h[i]/Hhu;
    v[i] = (D[i] - h[i-1]*v[i-1])/Hhu;
  }

  /* Back substitution */
  d2f[N0-2] = v[N0-2];
  for (i = N0-3; i > 0; i--) {
    d2f[i] = v[i] - u[i]*d2f[i+1];
  }

  /* Calculate spline coefficients */
  for (i = 0; i < N0-1; i++) {
    coeff[i][0] = (d2f[i]*x0[i+1]*x0[i+1]*x0[i+1] - d2f[i+1]*x0[i]*x0[i]*x0[i]
                  + 6.0*(f0[i]*x0[i+1] - f0[i+1]*x0[i]))/(6.0*h[i])
                  + h[i]/6.0*(d2f[i+1]*x0[i] - d2f[i]*x0[i+1]);
    coeff[i][1] = (d2f[i+1]*x0[i]*x0[i] - d2f[i]*x0[i+1]*x0[i+1] + 2.0*(f0[i+1] - f0[i]))
                  /(2.0*h[i]) + h[i]/6.0*(d2f[i] - d2f[i+1]);
    coeff[i][2] = (d2f[i]*x0[i+1] - d2f[i+1]*x0[i])/(2.0*h[i]);
    coeff[i][3] = (d2f[i+1] - d2f[i])/(6.0*h[i]);
  }


  /* Construct cumulative spline integral if flag is set */
  if (F0flag) {
    Ftot = 0.0;

    for (i = 0; i < N0-1; i++) {
      Fa = CSplineF(coeff, x0[i], i);
      Fb = CSplineF(coeff, x0[i+1], i);

      Ftot += Fb - Fa;
      F0[i] = Ftot;
    }
  }


  /* Free memory */
  Mem(MEM_FREE, d2f);
  Mem(MEM_FREE, h);
  Mem(MEM_FREE, d);
  Mem(MEM_FREE, H);
  Mem(MEM_FREE, D);
  Mem(MEM_FREE, u);
  Mem(MEM_FREE, v);

}
/*****************************************************************************/


/*****************************************************************************/
double CSplineF(double **coeff, double x, long idx) {
  /* Helper function */
  return x*(coeff[idx][0] + x*(coeff[idx][1]/2.0 + x*(coeff[idx][2]/3.0
         + x*coeff[idx][3]/4.0)));
}
/*****************************************************************************/


/*****************************************************************************/
void CSplineInterpolate(double **coeff, double *x0, long N0, double *x1,
                        double *f1, long N1) {
  /* Interpolation using cubic splines */

  long i, lo;

  /* Check the array size */
  if (N0 < 3)
    Die(FUNCTION_NAME, "Number of elements too small: %ld", N0);

  for (i = 0; i < N1; i++) {

    /* Find interval for xb */
    if ((lo = SearchArray(x0, x1[i], N0)) == -1) {
      if (x1[i] == x0[N0-1])
        lo = N0 - 2;
      else
        Die(FUNCTION_NAME, "CSplineInterpolate: x[%ld] %.5E not found\n", i, x1[i]);
    }

    f1[i] = coeff[lo][0] + x1[i]*(coeff[lo][1] + x1[i]*(coeff[lo][2]
            + x1[i]*coeff[lo][3]));

  }
}
/*****************************************************************************/


/*****************************************************************************/
void CSplineInterpolate0(double *x0, double *f0, long N0, double d2f0,
                         double d2fend, double *x1, double *f1, long N1,
                         long mode) {
  /* Wrapper for interpolation when spline coefficients etc. are not needed */

  long i;
  double tmp;
  double *lx0, *lf0, *lx1;
  double **coeff;

  /* Check the array size */
  if (N0 < 3)
    Die(FUNCTION_NAME, "Number of elements too small: %ld", N0);

  /* Allocate memory for the coefficients */
  coeff = (double **)Mem(MEM_ALLOC, N0-1, sizeof(double*));
  for (i = 0; i < N0-1; i++)
    coeff[i] = (double *)Mem(MEM_ALLOC, 4, sizeof(double));

  if (mode == 1) {

    /* Linear */

    /* Construct splines */
    CSplineConstruct(x0, f0, N0, d2f0, d2fend, coeff, 0, &tmp);

    /* Interpolate */
    CSplineInterpolate(coeff, x0, N0, x1, f1, N1);
  }
  else if (mode == 2) {

    /* Log-log */

    lx0 = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
    lf0 = (double *)Mem(MEM_ALLOC, N0, sizeof(double));
    for (i = 0; i < N0; i++) {
      lx0[i] = log(x0[i]);
      lf0[i] = log(f0[i]);
    }

    lx1 = (double *)Mem(MEM_ALLOC, N1, sizeof(double));
    for (i = 0; i < N1; i++)
      lx1[i] = log(x1[i]);

    /* Construct splines */
    CSplineConstruct(lx0, lf0, N0, d2f0, d2fend, coeff, 0, &tmp);

    /* Interpolate */
    CSplineInterpolate(coeff, lx0, N0, lx1, f1, N1);

    for (i = 0; i < N1; i++)
      f1[i] = exp(f1[i]);

    Mem(MEM_FREE, lx0);
    Mem(MEM_FREE, lf0);
    Mem(MEM_FREE, lx1);

  }
  else
    Die(FUNCTION_NAME, "Wrong mode");

  /* Free memory */
  for (i = 0; i < N0-1; i++)
    Mem(MEM_FREE, coeff[i]);
  Mem(MEM_FREE, coeff);

}
/*****************************************************************************/


/*****************************************************************************/
double CSplineIntegrate(double **coeff, double *x0, double *F0, long N0,
                     double xa, double xb) {
  /* Integration using cubic splines */

  long ia, ib;
  double Fa, Fb, Fab;

  /* Check the array size */
  if (N0 < 3)
    Die(FUNCTION_NAME, "Number of elements too small: %ld", N0);

  /* Find interval for xa */
  if ((ia = SearchArray(x0, xa, N0)) == -1) {
    if (xa == x0[N0-1])
      ia = N0 - 2;
    else
      Die(FUNCTION_NAME, "CSplineIntegrate: xa %.5E not found\n", xa);
  }

  /* Find interval for xb */
  if ((ib = SearchArray(x0, xb, N0)) == -1) {
    if (xb == x0[N0-1])
      ib = N0 - 2;
    else
      Die(FUNCTION_NAME, "CSplineIntegrate: xb %.5E not found\n", xb);
  }

  /* Calculate integral */
  Fa = CSplineF(coeff, xa, ia);
  Fb = CSplineF(coeff, xb, ib);

  if (ia == ib) {
    Fab = Fb - Fa;
  }
  else {
    Fa = CSplineF(coeff, x0[ia+1], ia) - Fa;
    Fb = Fb - CSplineF(coeff, x0[ib], ib);
    Fab = F0[ib-1] - F0[ia] + Fa + Fb;
  }

  return Fab;
}
/*****************************************************************************/


/*****************************************************************************/
double CSplineIntegrate0(double *x0, double *f0, long N0, double d2f0,
                         double d2fend, double xa, double xb) {
  /* Wrapper for integration when spline coefficients etc. are not needed */

  long i;
  double Fab;
  double *F0;
  double **coeff;

  /* Check the array size */
  if (N0 < 3)
    Die(FUNCTION_NAME, "Number of elements too small: %ld", N0);

  /* Allocate memory for the coefficients */
  coeff = (double **)Mem(MEM_ALLOC, N0-1, sizeof(double*));
  for (i = 0; i < N0-1; i++)
    coeff[i] = (double *)Mem(MEM_ALLOC, 4, sizeof(double));

  /* Allocate memory for integral array */
  F0 = (double *)Mem(MEM_ALLOC, N0, sizeof(double));

  /* Construct splines */
  CSplineConstruct(x0, f0, N0, d2f0, d2fend, coeff, 1, F0);

  /* Integrate */
  Fab = CSplineIntegrate(coeff, x0, F0, N0, xa, xb);

  /* Free memory */
  for (i = 0; i < N0-1; i++)
    Mem(MEM_FREE, coeff[i]);
  Mem(MEM_FREE, coeff);
  Mem(MEM_FREE, F0);

  return Fab;
}
/*****************************************************************************/
