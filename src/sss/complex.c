#ifdef __cplusplus 
extern "C" { 
#endif 
#include "header.h"

/* complex.c: modified 7.6.2014: (MPu) */
/* c_sqrt and c_exp added              */

/*****************************************************************************/

/* kompleksiluvut */

/*--------------------------------------------------------------------*/
complex c_add(complex x, complex y){
  complex z;
  z.re = x.re + y.re;
  z.im = x.im + y.im; 
  return z;
}
/*--------------------------------------------------------------------*/
complex c_sub(complex x, complex y){
  complex z;
  z.re = x.re - y.re;
  z.im = x.im - y.im;
  return z;
}
/*--------------------------------------------------------------------*/
complex c_mul(complex x, complex y){
  complex z;
  z.re = x.re * y.re - x.im * y.im; 
  z.im = x.re * y.im + x.im * y.re; 
  return z;
}
/*--------------------------------------------------------------------*/
double c_norm(complex x){
  double norm; 
  norm = x.re * x.re + x.im * x.im;
  norm = sqrt(norm); 
  return norm;
}
/*--------------------------------------------------------------------*/
complex c_div(complex x, complex y){
  complex z;
  double n = c_norm(y); 
  n = n*n;
  z.re = (x.re * y.re + x.im * y.im) / n;
  z.im = (x.im * y.re - x.re * y.im) / n; 
  return z;
}
/*--------------------------------------------------------------------*/
complex c_con(complex x){
  complex z;
  z.re = x.re;
  z.im = -x.im; 
  return z;
}
/*--------------------------------------------------------------------*/
complex c_sqrt(complex z){

  double x,y,r,theta; 
  double eps; 
  complex w;

  eps = 1E-30; 

  x = z.re; 
  y = z.im; 

  r = sqrt(x*x + y*y); 


  /* Case x = 0 :*/
  if (y > eps){
    theta =  0.5 * M_PI; /* Pi/2 */
  }
  if (y < -eps){
    theta = -0.5 * M_PI; /* -Pi/2 */
  }
  if (y > -eps && y < eps){ /* 0.0*/
    theta = 0.0; 
  }


  /* Case x != 0 */


  /* -Pi/2 < atan < Pi/2 */

  if (x > eps){
    theta = atan(y/x); /* Right half-plane of C */
  }
  if (x < -eps){ 
    if (y > eps){ /* second quarter, Pi/2 < theta < Pi  */
      theta = atan(abs(y/x)) + 0.5 * M_PI; 
    }
    if (y < -eps){ /* third quarter, -Pi < theta < -Pi/2 */
      theta = atan(abs(y/x)) - M_PI; 
    }
    if (-eps < y && y < eps){ /* y = 0, x < 0 => theta = Pi */
      theta = M_PI; 
    }
  }
  
  w.re = sqrt(r) * cos( 0.5 * theta ); 
  w.im = sqrt(r) * sin( 0.5 * theta ); 

  return w;
}
/*--------------------------------------------------------------------*/


complex c_exp(complex z){

  double x,y; 
  complex w;

  x = z.re; 
  y = z.im; 
  w.re = exp(x) * cos(y); 
  w.im = exp(x) * sin(y); 


  return w;
}
/*--------------------------------------------------------------------*/


/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
