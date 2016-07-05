/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : dfsolver.c                                     */
/*                                                                           */
/* Created:       2014/06/07 (MPu)                                           */
/* Last modified: 2014/06/07 (MPu)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/*                                                                           */
/* Description: solve homogeneous diffusion equation                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "dfsolver:"

void bsfunN(long, double *, double *, complex **, complex **, complex **);

/*****************************************************************************/
void dfSolver(long ifCorn, long nG, long nJ, double a, double dc, 
	      double **surfs_n, double **surfs_r, 
              complex **T, complex **U, complex *b, complex *c)
{
  long i,j,k,ii,jj,ind,ind1,ind2,Nr,nS;
  double r[3], tn[3], tn1[3], tn2[3]; 
  double *x,*y,*t, *x1, *x2, *y1, *y2; 
  double del, nrm, nrm1, nrm2; 
  complex z; 
  complex **E, **Mr, **M; 


  /* Integration parameter */

  Nr = 100; 

  if (nJ == 8 || nJ == 12)
    nS = nJ/2; 
  else
    nS = nJ; 

  
  /* Allocate temporary variables */
  x  = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  x1 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  x2 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  y  = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  y1 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  y2 = (double *)Mem(MEM_ALLOC, Nr, sizeof(double)); 
  t  = (double *)Mem(MEM_ALLOC, Nr, sizeof(double));
  
  E  = (complex **)Mem(MEM_ALLOC, nG   , sizeof(complex));
  Mr = (complex **)Mem(MEM_ALLOC, nG*nG, sizeof(complex)); 
  M  = (complex **)Mem(MEM_ALLOC, nG*nJ, sizeof(complex)); 
  for (i=0; i<nG; i++)
    E[i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex)); 
  for (i=0; i<nG*nG; i++)
    Mr[i] = (complex *)Mem(MEM_ALLOC, Nr, sizeof(complex)); 
  for (i=0; i<nG*nJ; i++)
    M[i] = (complex *)Mem(MEM_ALLOC, nG*nJ, sizeof(complex)); 


  /*---------------------------------------------------------------------*/

  /* Integration variables */

  del = a/(Nr-1); 
  for(i=1; i<Nr; i++){
    t[i] = t[i-1] + del; 
  }

  /* --- Form boundary conditions --- */

  for (i=0; i<nS; i++){ /* surfaces */ 

      /* Tangent vector of the surface  */

      tn[0] = -surfs_n[i][1]; /* x component */
      tn[1] =  surfs_n[i][0]; /* y component */
      tn[2] =  0.0; 

      /* Take Nr points from the boundary */

/*       printf("surfs_n = (%f %f %f)\n", surfs_n[i][0], surfs_n[i][1], surfs_n[i][2]);  */
/*       printf("surfs_r = (%f %f %f)\n", surfs_r[i][0], surfs_r[i][1], surfs_r[i][2]);  */

      for (k=0; k<Nr; k++){

	/* x(t) = x0 + tn_x * t */
	

	x[k] = surfs_r[i][0] + t[k] * tn[0]; 
	y[k] = surfs_r[i][1] + t[k] * tn[1]; 

      }

      for (j=0; j<nJ; j++){ /* functions */

      /* surfs_n[j] = basis function vector */
      /* surfs_n[i] = surface normal vector */

      /* Compute inner product of vectors */

      nrm = -surfs_n[j][0] * surfs_n[i][0] - surfs_n[j][1] * surfs_n[i][1] ; 

      /* Divide by surface area (this is for integration later) */
      
      nrm = nrm / a;


      /* Integrate over surface i:  */

      if ( nrm > 1E-20 || nrm < -1E-20){ /* nrm = 0.0 => zero matrix */

	for (k=0; k<Nr; k++){
	  
	  r[0] = x[k]; 
	  r[1] = y[k]; 
	  r[2] = 0.0; 

	  /* Normal derivative at r */

	  for (ii=0; ii<nG; ii++){
	    memset(E[ii], 0, nG*sizeof(complex)); 
	  }
	  bsfunN(nG, surfs_n[j], r, T, U, E); 

	  /* Store row by row */

	  ind = 0; 
	  for (ii=0; ii<nG; ii++){
	    for (jj=0; jj<nG; jj++){
	      Mr[ind][k] = E[ii][jj]; 
	      ind = ind + 1; 
	    }
	  }
	}  /* End of k loop */

	/* Integrate  => E */

	ind = 0; 
	for (ii=0; ii<nG; ii++){
	  for (jj=0; jj<nG; jj++){

	    z = trapz(Nr, t, Mr[ind]); 

	    E[ii][jj].re = z.re * nrm; 
	    E[ii][jj].im = z.im * nrm; 
	  
	    ind = ind + 1; 
	  }
	}

	/* Store to matrix M:  */

	/* -> functions */

	/* |           */
	/* V  surfaces */

	for (ii=0; ii<nG; ii++){
	  for (jj=0; jj<nG; jj++){

	    M[i*nG + ii][j*nG + jj] = E[ii][jj]; 

	  }
	}
   
      } /* if nrm != 0 ...*/

    } /* functions */

  } /* surfaces */

  /*-----------------------------------------------------------------*/


  if (ifCorn == 1){

    /* --- Corner integration --- */

    del = (dc/2.0)/(Nr-1);  /* half corner */
    for(i=1; i<Nr; i++){
      t[i] = t[i-1] + del; 
    }

    for (i=0; i<nS; i++){ /* corners */

      /* First  part of corner on surface i     */
      /* Second part of corner on surface i + 1 */

      ind1 = i; 
      if (i == nS-1){  /* Last corner */
	ind2 = 0; 
      }
      else{
	ind2 = i+1; 
      }
      

      /* Tangent vectors  */

      /* Surface 1: */

      tn1[0] = -surfs_n[ind1][1]; /* x component */
      tn1[1] =  surfs_n[ind1][0]; /* y component */

      /* Surface 2 */

      tn2[0] = -surfs_n[ind2][1]; /* x component */
      tn2[1] =  surfs_n[ind2][0]; /* y component */

      /* Take Nr points */

      for (k=0; k<Nr; k++){

	x1[k] = surfs_r[ind1][0] + (0.95 * a + t[k]) * tn1[0]; 
	y1[k] = surfs_r[ind1][1] + (0.95 * a + t[k]) * tn1[1]; 

	x2[k] = surfs_r[ind2][0] + t[k] * tn2[0]; 
	y2[k] = surfs_r[ind2][1] + t[k] * tn2[1];
      }

      for (j=0; j<nJ; j++){ /* functions */     

	/* surfs_n[j] = basis function vector */
	/* surfs_n[i] = surface normal vector */


	/* Compute inner product of vectors */

	nrm1 = -surfs_n[j][0] * surfs_n[ind1][0] - surfs_n[j][1] * surfs_n[ind1][1] ; 
	nrm2 = -surfs_n[j][0] * surfs_n[ind2][0] - surfs_n[j][1] * surfs_n[ind2][1] ; 

	/* Divide by surface area (this is for integration later) */
      
	nrm1 = nrm1 / dc;
	nrm2 = nrm2 / dc;


	/* Integrate over boundary area */


	/* -- Firt part of the corner -- */

	if ( nrm1 > 1E-20 || nrm1 < -1E-20){ /* nrm = 0.0 => zero matrix */

	  /* Compute bsfun value at Nr points */

	  for (k=0; k<Nr; k++){
	  
	    r[0] = x1[k]; 
	    r[1] = y1[k]; 

	    /* Normal derivative at r */

	    for (ii=0; ii<nG; ii++){
	      memset(E[ii], 0, nG*sizeof(complex)); 
	    }
	    bsfunN(nG,surfs_n[j], r, T, U, E);


	    /* Store row by row */

	    ind = 0; 
	    for (ii=0; ii<nG; ii++){
	      for (jj=0; jj<nG; jj++){
		Mr[ind][k] = E[ii][jj]; 
		ind = ind + 1; 
	      }
	    }
	  }  /* End of k loop */

	  /* Integrate  => E*/

	  ind = 0; 
	  for (ii=0; ii<nG; ii++){
	    for (jj=0; jj<nG; jj++){

	      z = trapz(Nr, t, Mr[ind]); 

	      E[ii][jj].re = z.re * nrm1; 
	      E[ii][jj].im = z.im * nrm1; 
	  
	      ind = ind + 1; 
	    }
	  }

	  /* Store to matrix M:  */
     
	  for (ii=0; ii<nG; ii++){
	    for (jj=0; jj<nG; jj++){
	      M[nS*nG + i*nG + ii][j*nG + jj] = E[ii][jj]; 
	    }
	  }         
	} /* if nrm1 != 0 ...*/


	/* -- Second part of the corner -- */

	if ( nrm2 > 1E-20 || nrm2 < -1E-20){ /* nrm = 0.0 => zero matrix */

	  /* Compute bsfun value at Nr points */
	
	  for (k=0; k<Nr; k++){
	  
	    r[0] = x2[k]; 
	    r[1] = y2[k]; 

	    /* Normal derivative at r */

	    for (ii=0; ii<nG; ii++){
	      memset(E[ii], 0, nG*sizeof(complex)); 
	    }
	    bsfunN(nG,surfs_n[j], r, T, U, E);


	    /* Store row by row */

	    ind = 0; 
	    for (ii=0; ii<nG; ii++){
	      for (jj=0; jj<nG; jj++){
		Mr[ind][k] = E[ii][jj];        
		ind = ind + 1; 
	      }
	    }
	  }  

	  /* Integrate  => E */

	  ind = 0; 
	  for (ii=0; ii<nG; ii++){
	    for (jj=0; jj<nG; jj++){

	      z = trapz(Nr, t, Mr[ind]); 

	      E[ii][jj].re = z.re * nrm2; 
	      E[ii][jj].im = z.im * nrm2; 
	  
	      ind = ind + 1; 
	    }
	  }


	  /* Store to matrix M:  */

	  for (ii=0; ii<nG; ii++){
	    for (jj=0; jj<nG; jj++){

	      M[nS*nG + i*nG + ii][j*nG + jj] = c_add( M[nS*nG + i*nG + ii][j*nG + jj], E[ii][jj] ); 

	      CheckValue(FUNCTION_NAME, "M[ii][jj].re","", M[nS*nG + i*nG + ii][j*nG + jj].re, -INFTY, INFTY);
	      CheckValue(FUNCTION_NAME, "M[ii][jj].im","", M[nS*nG + i*nG + ii][j*nG + jj].im, -INFTY, INFTY);	      

	    }
	  }      
   
	} /* if nrm2 != 0 ...*/
      } /* functions */
    } /* corners */
  } /* End of "if(ifCorn == 1)" */


  /*
  for (jj=0; jj<nJ*nG; jj++){
    for (ii=0; ii<nJ*nG; ii++){
      printf("M(%ld,%ld) = %e %e\n", ii+1, jj+1, M[ii][jj].re, M[ii][jj].im);
    }
  }

  */
  /* Solve c from  Mc = b */

  
  LUdecomposition(nG*nJ, M, b, c);
  /*

  for (i=0; i<nG*nJ; i++)
    printf("c(%ld) = (%e,%e)\n", i, c[i].re, c[i].im); 
  */

  /* OK, everything done here! */

  /*-----------------------------------------------------------------*/

  /* Free temporary variables */

  Mem(MEM_FREE, x); 
  Mem(MEM_FREE, x1); 
  Mem(MEM_FREE, x2); 
  Mem(MEM_FREE, y); 
  Mem(MEM_FREE, y1); 
  Mem(MEM_FREE, y2); 
  Mem(MEM_FREE, t); 

  for (i=0; i<nG; i++)
    Mem(MEM_FREE, E[i]);
  for (i=0; i<nG*nG; i++)
    Mem(MEM_FREE, Mr[i]);
  for (i=0; i<nG*nJ; i++)
    Mem(MEM_FREE, M[i]);

  Mem(MEM_FREE, E);
  Mem(MEM_FREE, Mr);  
  Mem(MEM_FREE, M);  

 /*-----------------------------------------------------------------*/
  

  /* Return */

  return; 

}

/**************************************************************************************/
void bsfunN(long nG, double *v, double *r, complex **T, complex **U, complex **E){

  /* Compute normal derivative of DE basis function at r */
  /* function: expm( v'*r sqrtm(A) )                     */
  /* A = U*T*U', T tirangular, U unitary                 */

  long i, j, dim; 
  double prod; 
  complex z; 

  complex *d, **As,**B; 


  /* Allocate auxiliary variables */

  d  = (complex * )Mem(MEM_ALLOC, nG, sizeof(complex  )); 
  As = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *)); 
  B  = (complex **)Mem(MEM_ALLOC, nG, sizeof(complex *)); 
  for (i=0; i<nG; i++){
    As[i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex)); 
    B [i] = (complex *)Mem(MEM_ALLOC, nG, sizeof(complex)); 
  }

  /*-----------------------------------------------------------------------*/


  /* 3-D space */

  dim = 3; 

  /* Compute inner product v' * r:  */ 

  prod = 0.0; 
  for (i=0; i<dim; i++){
    prod = prod + v[i] * r[i]; 
  }

  /*-----------------------------------------------------------------------*/

  /* Use Parlett method to compute e^( prod * sqrt( A ) ) */


  /* Compute B = expm (prod * sqrtm (A)) */

  if (prod > 1E-20 || prod <-1E-20){ /* Compute normal derivative on the boundary */

    for (i=0; i<nG; i++){
      z    = T[i][i]; 
      z    = c_sqrt(z); 
      z.re = z.re * prod; 
      z.im = z.im * prod; 

      CheckValue(FUNCTION_NAME, "z.re","", z.re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "z.im","", z.im, -INFTY, INFTY);
      
      d[i] = c_exp(z);     
    }

  /* Parlett algorithm */

    parlett(nG, d, T, U, B); 

  }
  else{ /* expm of zero  matrix => identity matrix */
  
    for (i=0; i<nG; i++){
      B[i][i].re = 1.0; 
    }    
  }

  for (i=0; i<nG; i++){
    for (j=0; j<nG; j++){
      CheckValue(FUNCTION_NAME, "B[i][j].re","", B[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "B[i][j].im","", B[i][j].im, -INFTY, INFTY);
    }
  }


  /* Compute As = sqrtm( A ) */

  for (i=0; i<nG; i++){
    z    = T[i][i]; 
    d[i] = c_sqrt(z); 
  }

   
  parlett(nG, d, T, U, As); 

  for (i=0; i<nG; i++){
    for (j=0; j<nG; j++){
      CheckValue(FUNCTION_NAME, "As[i][j].re","", As[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "As[i][j].im","", As[i][j].im, -INFTY, INFTY);
    }
  }

  matProduct(nG,nG,nG,As,B,E); 

  for (i=0; i<nG; i++){
    for (j=0; j<nG; j++){
      CheckValue(FUNCTION_NAME, "E[i][j].re","", E[i][j].re, -INFTY, INFTY);
      CheckValue(FUNCTION_NAME, "E[i][j].im","", E[i][j].im, -INFTY, INFTY);
    }
  }

  /* Free auxiliary variables */  

  Mem(MEM_FREE, d); 
  for (i=0; i<nG; i++){
    Mem(MEM_FREE,As[i]);
    Mem(MEM_FREE,B[i]);      
  }

  Mem(MEM_FREE,As); 
  Mem(MEM_FREE,B); 

  /* Return */

  return; 

}

/**************************************************************************************/
