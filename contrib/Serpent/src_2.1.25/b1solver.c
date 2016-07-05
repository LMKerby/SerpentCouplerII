/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : b1solver.c                                     */
/*                                                                           */
/* Created:       2011/06/11 (JLe)                                           */
/* Last modified: 2013/04/04 (JLe)                                           */
/* Version:       2.1.20                                                     */
/*                                                                           */
/* Description: Forms and solves B1 equations                                */
/*                                                                           */
/* Comments: - Chi should match infinite spectrum calculation?               */
/*                                                                           */
/*           - Micro-group cross sections are calculated in                  */
/*             calcmicrogroupxs.c and cleared in clearmicrogroupxs.c         */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "B1Solver:"

/* Setting maximum iterations to zero should reproduce infinite */
/* spectrum results. */

#define MAX_ITER 1000
#define B2_LIM 1E-8

void PrintB1Data(long, double *, complex *, double *, double *, double *, 
		 double *, double *, complex *, double **, double **, 
		 double **);

/*****************************************************************************/

void B1Solver()
{
  return;

  Die(FUNCTION_NAME, "Tää pitää korjata");

#ifdef mmmmmmmmmmmmmmmmmm

  long ptr, nmg, nfg, gcu, uni, n, m, i, j, iter, nnz;
  double *flx0, *tot, *fiss, *abs, *nsf, *mubar, *x, *alpha, J, *Db1, **p0;
  double *Xep, *Ip, *Smp, *Pmp, *Ia, *Xea, *Pma, *Sma;
  const double *imap;
  double **np0, **np1, norm, flx, sum1, sum2, val, kinf, L2, B2, B, err, kb1;
  complex theta0, *ei, *di, **D, *chi, *flx1,*LUval, *chival;
  struct ccsMatrix *Dinv, *Dinv2, *A; 

  /* Check active cycle */
  
  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

  /* Check if group constants are generated */

  if((long)RDB[DATA_OPTI_GC_CALC] == NO)
    return;

  /* Check B1 calculation */

  if ((long)RDB[DATA_B1_CALC] == NO)
    return;

  /* Check micro-group mode */

  if (((long)RDB[DATA_SIMULATION_COMPLETED] == NO) && 
      ((long)RDB[DATA_MICRO_CALC_MODE] == MICRO_CALC_MODE_END))
    return;
  else if (((long)RDB[DATA_SIMULATION_COMPLETED] == YES) && 
	   ((long)RDB[DATA_MICRO_CALC_MODE] == MICRO_CALC_MODE_CYCLE))
    return;

  /* Get pointer to micro-group structure */
  
  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of groups */

  nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
  nfg = (long)RDB[DATA_ERG_FG_NG];

  /* Index map */

  ptr = (long)RDB[DATA_MICRO_PTR_IDX_MAP];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  imap = &RDB[ptr];

  /* Reduce results buffer */

  WDB[DATA_RES2_REDUCED] = (double)NO;
  ReducePrivateRes();

  /* Allocate memory for matrix np0 */

  np0 = (double **)Mem(MEM_ALLOC, nmg, sizeof(double *));

  for(n = 0; n < nmg; n++)
    np0[n] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));

  /* Allocate memory for matrix np1 */

  np1 = (double **)Mem(MEM_ALLOC, nmg, sizeof(double *));

  for(n = 0; n < nmg; n++)
    np1[n] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));

  /* Allocate memory for matrix p0 */

  p0 = (double **)Mem(MEM_ALLOC, nmg, sizeof(double *));

  for(n = 0; n < nmg; n++)
    p0[n] = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));

  /* Allocate memory for temporary variables */

  x = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));
  alpha = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));
  
  /* Ei tartte kopioida koko LU-matriisia vaan ainoastaan values joka kerta*/

  LUval   = (complex *)Mem(MEM_ALLOC, nmg*nmg, sizeof(complex));

  /* Allocate memory for matrix etc. */

  Dinv = ccsMatrixNew(nmg, nmg, nmg*nmg);
  Dinv2 = ccsMatrixNew(nmg, nmg, nmg*nmg);
  A = ccsMatrixNew(nmg, nmg, nmg*nmg);
  ei = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));   
  di = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));   
  chi = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));  
  chival = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));  
  flx1 = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));   
  Db1 = (double *)Mem(MEM_ALLOC, nmg, sizeof(double));   

  /* Matriisi D kannattaa tallentaa suoraan taulukkona */

  D = (complex **)Mem(MEM_ALLOC, nmg, sizeof(complex *)); 

  for(n = 0; n < nmg; n++)
    D[n] = (complex *)Mem(MEM_ALLOC, nmg, sizeof(complex));

  /* Get infinite flux normalization factor (tätä pitää kutsua täällä */
  /* koska NormCoef():sta kutsutaan ReduceBuffer():ia, joka asettaa   */
  /* seuraavassa resetoitavan flägin.) */

  if ((long)RDB[DATA_SIMULATION_COMPLETED] == NO)
    norm = NormCoef(PARTICLE_TYPE_NEUTRON);
  else
    {
      ptr = (long)RDB[RES_NORM_COEF];
      norm = Mean(ptr, 0)/RDB[DATA_CRIT_CYCLES];
    }

  /* Mark scoring buffer unreduced (tän pitäisi toimia koska  */
  /* käsitellään sellasia muuttujia joita ei lisätä bufferiin */
  /* muualla). */
  
  WDB[DATA_BUF_REDUCED] = (double)NO;

  /***************************************************************************/

  /***** Perform B1 calculation **********************************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Pointer to universe */

      uni = (long)RDB[gcu + GCU_PTR_UNIV];
      CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);

      if ((long)RDB[DATA_SIMULATION_COMPLETED] == YES)
	fprintf(out, "Performing B1 iteration in universe %s...\n\n",
		GetText(uni + UNIVERSE_PTR_NAME));

      /***********************************************************************/

      /***** Calculate micro-group cross sections ****************************/
      
      /* Get pointers cross section data */

      ptr = (long)RDB[gcu + GCU_MICRO_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      flx0 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_TOT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      tot = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_FISS];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      fiss = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      abs = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_NSF];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      nsf = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_MUBAR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      mubar = &RES2[ptr];

      /* Check poison calculation flag */

      ptr = (long)RDB[gcu + GCU_MICRO_I135_PROD];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Ip = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_XE135_PROD];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Xep = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_PM149_PROD];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Pmp = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_SM149_PROD];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Smp = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_I135_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Ia = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_XE135_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Xea = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_PM149_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Pma = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SM149_ABS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      Sma = &RES2[ptr];
    
      /* Read data to scattering matrixes */

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT_PROD_MTX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

      for (n = 0; n < nmg; n++)
	for (m = 0; m < nmg; m++)
	  np0[n][m] = RES2[ptr + m*nmg + n];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT_MTX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);

      for (n = 0; n < nmg; n++)
	for (m = 0; m < nmg; m++)
	  p0[n][m] = RES2[ptr + m*nmg + n];

      /* Read fission spectrum */

      ptr = (long)RDB[gcu + GCU_MICRO_CHI];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      
      for (n = 0; n < nmg; n++)
	{
	  chi[n].re = RES2[ptr + n];
	  chi[n].im = 0.0;
	}

      /* Calculate np1 scattering matrix */
      
      for (n = 0; n < nmg; n++)
	for (m = 0; m < nmg; m++)
	  {
	    if (flx0[m] > 0.0)
	      np1[n][m] = mubar[m]*np0[n][m];
	    else
	      np1[n][m] = 0.0;
	  }

      /***********************************************************************/

      /***** Calculate initial guess for kinf and B2 *************************/

      /* Calculate source term and loss terms*/

      sum1 = 0.0;
      sum2 = 0.0;

      for (n = 0; n < nmg; n++)
	{
	  sum1 = sum1 + flx0[n]*nsf[n];
	  sum2 = sum2 + flx0[n]*abs[n];
	}

      /* Check values */

      CheckValue(FUNCTION_NAME, "sum1", "", sum1, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "sum2", "", sum2, ZERO, INFTY);

      /* Calculate k-inf */

      kinf = sum1/sum2;

      /* Calculate absorption xs */

      sum1 = 0.0;
      sum2 = 0.0;

      for (n = 0; n < nmg; n++)
	{
	  sum1 = sum1 + flx0[n]*abs[n];
	  sum2 = sum2 + flx0[n];
	}

      /* Check values */

      CheckValue(FUNCTION_NAME, "sum1", "", sum1, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "sum2", "", sum2, ZERO, INFTY);
      
      val = sum1/sum2;

      /* Calculate diffusion coefficient (approximate removal by absorption) */

      sum1 = 0.0;
      sum2 = 0.0;

      for (n = 0; n < nmg; n++)
	{
	  sum1 = sum1 + flx0[n];
	  sum2 = sum2 + 3.0*flx0[n]*(tot[n] - mubar[n]*(tot[n] - abs[n]));
	}

      /* Check values */

      CheckValue(FUNCTION_NAME, "sum1", "", sum1, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "sum2", "", sum2, ZERO, INFTY);

      /* Migration area */

      L2 = sum1/sum2/val;

      /* Buckling */

      B2 = (kinf - 1.0)/L2;

      /* Toi ylläoleva on materiaalin kupevuus, jonka pitäisi olla vähän */
      /* lähempänä geometrista kupevuutta johon iteraatio konvergoituu */
      /* kuin kiinnitetyn alkuarvauksen, joka saattaa olla monta */
      /* kertaluokkaa pielessä. Toimii ainakin Serpent 1:ssä. Jos tämä */
      /* rutiini toteutetaan niin että sitä kutsutaan jokaisen sukupolven */
      /* jälkeen, niin alkuarvauksena kannattaa käyttää edellisen sukupolven */
      /* konvergoitunutta arvoa (JLe). */
      
      /* Asetetaan alkuarvaus s.e. 0 <= x <= 1*/
      /*
      if (kinf >= 1.0)
	{
	B2 = 1e-5;
	val = B2; 
	for (n=0; n < nmg; n++)
	  if (tot[n]*tot[n] < val)
	    val = tot[n]*tot[n]; 
	B2 = val; 
	}
      else
	{
	B2 = -1e-5;
	val = -B2; 
	for (n=0; n < nmg; n++)
	  if (tot[n]*tot[n] < val)
	    val = tot[n]*tot[n]; 
	B2 = -val; 	
	}
      */

      /***********************************************************************/

      /***** Main loop *******************************************************/

      /* Avoid compiler warning */

      kb1 = kinf;

      /* Reset temporary variables */

      memset(x, 0.0, nmg*sizeof(double));
      memset(alpha, 0.0, nmg*sizeof(double));

      /* Reset error */

      err = 1.0;

      /* Iteration loop */

      /* Kaikki varsinainen laskenta tapahtuu tän silmukan sisällä.    */
      /* MAX_ITER on nyt asetetettu nollaksi, eli silmukkaa ei käydä   */
      /* ollenkaan läpi. Tuloksena saadut ryhmävakiot ovat tällöin     */
      /* NSF:ää lukuunottamatta numeroarvoltaan samat kuin infinite-   */
      /* spectrum -laskussa, mitä voi hyvin käyttää myöhemmin tehtävän */
      /* uudelleenhomogenisoinnin tarkastamiseen. Toi silmukan sisällä */
      /* tehtävä laskenta on vielä vähän mitä sattuu. Pääsin alkuun,   */
      /* mutta en saanut kaikista osista tolkkua. */

      for (iter = 0; iter < MAX_ITER; iter++)
	{
	  /*******************************************************************/

	  /***** Set temporary variables *************************************/

	  /* Check if too close to criticality */
	  
	  if (fabs(B2) < B2_LIM)
	    break; 
	
	  /* Square root of buckling */

	  if (B2 >= 0.0)
	    B = sqrt(B2);
	  else
	    B = sqrt(-B2);

	  /* Calculate temporary variable x (toi voi mennä INF:ksi!) */

	  for (n = 0; n < nmg; n++)
	    {	      
	      x[n] = B/tot[n];

	      /* NOTE: Tää voi mennä huonon statistiikan takai yli 1, */
	      /* mikä aiheuttaa NAN:in tuolla myöhemmin jos B2 < 0.   */
	      /* Katkaistaan tällöin luuppi (JLE 2.1.4 / 3.4.2012).   */

	      if (x[n] > 1.0)
		break;

	      CheckValue(FUNCTION_NAME, "x[n]", "", x[n], 0.0, 1.0);
	    }       

	  /* Tarkistetaan jos luuppi katkaistu, asetetaan iteraatioiden */
	  /* määrä maksimiin --> spektriä ei päivitetä, ja katkaistaan  */
	  /* iteraatio (JLE 2.1.4 / 3.4.2012) */

	  if (n < nmg)
	    {
	      iter = MAX_ITER;
	      break;
	    }
	
	  /* Calculate alpha  */
	  
	  if (B2 > 0.0)
	    { 	   
	      for (n = 0; n < nmg; n++)
		{
		  alpha[n] = x[n]*x[n];
		  alpha[n] = alpha[n]*atan(x[n])/(x[n] - atan(x[n])); 
		  CheckValue(FUNCTION_NAME, "alpha[n] (1)", "", alpha[n], 
			     -INFTY, INFTY);
		}
	    }
	  else if(B2 < 0.0)
	    {
	      for (n = 0; n < nmg; n++)
		{
		  alpha[n] = x[n]*x[n];
		  alpha[n] = alpha[n]*log((1.0 + x[n])/(1.0 - x[n]))/ 
		    (log((1.0 + x[n])/(1.0 - x[n])) - 2.0*x[n]);
		  
		  CheckValue(FUNCTION_NAME, "alpha[n] (2)", "", alpha[n], 
			     -INFTY, INFTY);
		}
	    }
	  else
	    {
	      for (n = 0; n < nmg; n++)
		alpha[n] = 3.0; 
	    }

	  /*******************************************************************/

	  /***** Form matrix Dinv ********************************************/

	  /* Reset size and pointer to first column */

	  nnz = 0;
	  Dinv->colptr[0] = 0 ;

	  /* Loop over columns */

	  for (m = 0; m < nmg; m++)
	    {
	      /* Loop over rows */

	      for (n = 0; n < nmg; n++)
		{
		  /* Check that element is non-zero */
		  
		  CheckValue(FUNCTION_NAME, "np1[n][m]", "", np1[n][m], -INFTY,
			     INFTY);

		  /* Set value */

		  Dinv->values[nnz].re = -3.0*np1[n][m];

		  /* Add diagonal */
		  
		  if (n == m)
		    {
		      CheckValue(FUNCTION_NAME, "alpha[n]", "", alpha[n], 
				 -INFTY, INFTY);
		      CheckValue(FUNCTION_NAME, "tot[m]", "", tot[m], 
				 -INFTY, INFTY);
		      Dinv->values[nnz].re = Dinv->values[nnz].re +
			alpha[n]*tot[m];
		    }
		  
		  /* Reset imagynary part */
		  
		  Dinv->values[nnz].im = 0.0;
		  
		  CheckValue(FUNCTION_NAME, "Dinv->values[nnz].re", "", 
			     Dinv->values[nnz].re, -INFTY, INFTY);
		  CheckValue(FUNCTION_NAME, "Dinv->values[nnz].im", "", 
			     Dinv->values[nnz].im, 0.0, 0.0);
		  
		  /* Put row index */
		  
		  Dinv->rowind[nnz] = n;
		  
		  /* Update size */
		  
		  nnz++;
		}	
	      
	      /* Put column index */
	      
	      Dinv->colptr[m + 1] = nnz;
	    }	  
	  
	  /*******************************************************************/
	  
	  /***** Form matrix D ***********************************************/
	
#ifdef DEBUG

	  for (i=0; i < Dinv->nnz; i++)
	    {
	      CheckValue(FUNCTION_NAME, "Dinv->values[i].re", "", 
			 Dinv->values[i].re, -INFTY, INFTY);
	      CheckValue(FUNCTION_NAME, "Dinv->values[i].im", "", 
			 Dinv->values[i].im, 0.0, 0.0);
	    }

#endif

	  /* Ota values talteen */
	
	  memcpy(LUval, Dinv->values, (Dinv->nnz) * sizeof(complex));

#ifdef DEBUG

	  for (i=0; i < Dinv->nnz; i++)
	    {
	      CheckValue(FUNCTION_NAME, "LUval.re", "", LUval[i].re, 
			 -INFTY, INFTY);
	      CheckValue(FUNCTION_NAME, "LUval.im", "", LUval[i].im, 
			 0.0, 0.0);
	    }

#endif

	  /* Set theta to zero */
	  
	  theta0.re = 0.0; 
	  theta0.im = 0.0; 

	  FindRowIndexes(Dinv);

	  /* Loop over energy groups */
	  
	  for (n = 0; n < nmg; n++)
	    {
	      /* Form ei */
	      
	      memset(ei, 0.0, nmg*sizeof(complex));
	      ei[n].re = 1.0; 
	      
	      memcpy(Dinv->values, LUval, (Dinv->nnz) * sizeof(complex));
	      
	      /* D:n sarakkeille varattu tila valmiiksi, Gauss ei */
	      /* varaa mitään ! */
	      
	      NumericGauss(Dinv, ei, theta0, D[n]); 

	    }
	  
	  /*******************************************************************/
	  
	  /***** Form matrix A ***********************************************/
	  
	  /* Reset size and pointer to first column */
	  
	  nnz = 0;
	  A->colptr[0] = 0 ;
	  
	  /* Loop over columns */
	  
	  for (m = 0; m < nmg; m++)
	    {
	      /* Loop over rows */
	      
	      for (n = 0; n < nmg; n++)
		{
		  /* Check that element is non-zero */
		  
		  /* if ((n == m) || (np0[n][m] != 0.0) || D[m][n].re != 0.0)*/
		  
		  /* Set value */
		  
		  A->values[nnz].re = -np0[n][m] + B2*D[m][n].re; /* D^T*/
		  
		  /* Add diagonal */
		  
		  if (n == m)
		    A->values[nnz].re = A->values[nnz].re + tot[n];
		  
		  /* Reset imaginary part */
		  
		  A->values[nnz].im = 0.0;
		  
		  /* Put row index */
		  
		  A->rowind[nnz] = n;
		  
		  /* Update size */
		  
		  nnz++;
		}

	      /* Put column index */
	      
	      A->colptr[m + 1] = nnz;
	    }
	  
	  /* Set size */
	  
	  A->nnz = nnz;
	  
	  /*******************************************************************/
	  
	  /***** Solve equation **********************************************/
	  
	  /* Solve flux */
	  
	  FindRowIndexes(A);
	  memcpy(chival, chi, nmg*sizeof(complex));

	  /* HUOM! Myös chival gaussataan... */ 

	  NumericGauss(A, chival, theta0, flx1); 
	  
	  /*B1 keff */
	  
	  kb1 = 0.0;
	  for (n = 0; n < nmg; n++)
	    kb1 = kb1 + nsf[n]*flx1[n].re; 
	  
	  /* Interpolointi */
	  
	  val = B2/(1.0/kb1 - 1.0/kinf); 
	  B2 = B2 + val - val/kb1; 

	  if ((long)RDB[DATA_SIMULATION_COMPLETED] == YES)
	    fprintf(out, "B2 = %11.5E, k = %1.5f \n", B2, kb1);  
	  
	  /* Calculate error */

	  err = kb1 - 1.0; 

	  /* Check convergence */
	  
	  if (fabs(err) < RDB[DATA_FUM_ERR_LIMIT])
	    break; 
	  
	  /*******************************************************************/
	}
    
      /***********************************************************************/

      /***** Normalize flux spectrum *****************************************/

      /* Check convergence */
      
      if (iter < MAX_ITER)
	{
	  /* Normalize infinite spectrum */

	  for (n = 0; n < nmg; n++)
	    flx0[n] = flx0[n]*norm;

	  /* Reset sums */
	  
	  sum1 = 0.0;
	  sum2 = 0.0;
	  
	  /* Check some specific normalizations */      
	  
	  if ((ptr = (long)RDB[DATA_PTR_NORM]) > VALID_PTR)
	    {
	      /* Constant flux normalization */
	      
	      if (RDB[ptr + NORM_FLUX] >= 0.0)
		{
		  /* Calculate sum over critical and infinite spectrum */
		  
		  for (n = 0; n < nmg; n++)
		    {
		      sum1 = sum1 + flx0[n];
		      sum2 = sum2 + flx1[n].re;
		    }
		}
	      
	      /* Normalization to generation rate */
	      
	      else if (RDB[ptr + NORM_GENRATE] >= 0.0)
		{
		  /* Calculate integral fission neutron production rate in */
		  /* critical and infinite spectrum */
		  
		  for (n = 0; n < nmg; n++)
		    {
		      sum1 = sum1 + nsf[n]*flx0[n];
		      sum2 = sum2 + nsf[n]*flx1[n].re;
		    }
		}
	      
	      /* Normalizatio to absorption rate */
	      
	      else if (RDB[ptr + NORM_ABSRATE] >= 0.0)
		{
		  /* Calculate integral absorption rate in critical and */
		  /* infinite spectrum */
		  
		  for (n = 0; n < nmg; n++)
		    {
		      sum1 = sum1 + abs[n]*flx0[n];
		      sum2 = sum2 + abs[n]*flx1[n].re;
		    }
		}
	    }
	  
	  /* Check sum (otherwise it is assumed that normalization is */
	  /* carried out by preserving some quantity proportional to  */
	  /* fission rate, such as power or power density.) */

	  if (sum2 == 0.0)
	    {
	      /* Calculate integral fission rate in critical and infinite */
	      /* spectrum */
	      
	      for (n = 0; n < nmg; n++)
		{
		  sum1 = sum1 + fiss[n]*flx0[n];
		  sum2 = sum2 + fiss[n]*flx1[n].re;
		}
	    }
	  
	  /* Check sums */

	  CheckValue(FUNCTION_NAME, "sum1", "", sum1, ZERO, INFTY);
	  CheckValue(FUNCTION_NAME, "sum2", "", sum2, ZERO, INFTY);
	  
	  /* Normalize critical flux*/
	  
	  for (n = 0; n < nmg; n++)
	    flx1[n].re = flx1[n].re*sum1/sum2;
	}

      /***********************************************************************/
    
      /***** Check convergence ***********************************************/

      /* Print for debugging */
      /*
      PrintB1Data(nmg, flx0, flx1, tot, abs, fiss, nsf, mubar, chi, p0, np0, 
		  np1);
      */
      /* Check convergence */

      if (iter == MAX_ITER)
	{
	  /* No convergence */
	  
	  if ((long)RDB[DATA_SIMULATION_COMPLETED] == YES)
	    fprintf(out, "No convergence, final error = %1.5E.\n\n", 
		    fabs(err));

	  /* Reset flux correction factors */

	  ptr = (long)RDB[gcu + GCU_FUM_PTR_SPEC_CORR];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  for (n = 0; n < nmg; n++)
	    WDB[ptr + n] = 1.0;
	}
      else
	{
	  /* Check status */

	  if ((long)RDB[DATA_SIMULATION_COMPLETED] == YES)
	    {
	      if (fabs(B2) < B2_LIM)
		fprintf(out, 
			"\nToo close to criticality for B1 iteration.\n\n");
	      else
		fprintf(out, "\nConvergence reached, final error = %1.5E.\n\n", 
			fabs(err));
	    }

	  /* Set flux correction factors */

	  ptr = (long)RDB[gcu + GCU_FUM_PTR_SPEC_CORR];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  for (n = 0; n < nmg; n++)
	    {
	      if (flx0[n] > 0.0)
		WDB[ptr + n] = flx1[n].re/flx0[n];
	      else
		WDB[ptr + n] = 1.0;
	    }

	  /* Copy multiplication factor */

	  kinf = kb1;
	  
	  /* Copy B1 flux */
	  
	  for (n = 0; n < nmg; n++)
	    flx0[n] = flx1[n].re;
	}

      /***********************************************************************/
      
      /***** Re-homogenize cross sections ************************************/

      /* Calculate multi-group diffusion coefficient */

      for (n = 0; n < nmg; n++)
	{
	  /* Current */

	  J = 0.0;
	  for (m = 0; m < nmg; m++)
	    J = J + D[m][n].re*flx0[m];  
      
	  /* Diffusion coefficient */
      
	  if (flx0[n] > 0.0)
	    Db1[n] = J/flx0[n];
	  else
	    Db1[n] = 0.0;
	}

      /* Calculate few-group reaction rates */

      for (n = 0; n < nmg; n++)
	{
	  /* Get few-group index */
	  
	  i = (long)imap[n];
	  
  	  /* Add reaction rates to buffers */

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_FLUX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n], 1.0, ptr, 0, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_TOTXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*tot[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*tot[n], 1.0, ptr, 0, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_ABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*abs[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*abs[n], 1.0, ptr, 0, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_FISSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*fiss[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*fiss[n], 1.0, ptr, 0, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_NSF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*nsf[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*nsf[n], 1.0, ptr, 0, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_DIFFCOEF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Db1[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Db1[n], 1.0, ptr, 0, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_CHI];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(chi[n].re, 1.0, ptr, 0, i);

	  /* Check poison calculation flag */

	  ptr = (long)RDB[gcu + GCU_FUM_FG_I135_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Ip[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Ip[n], 1.0, ptr, 0, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_XE135_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Xep[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Xep[n], 1.0, ptr, 0, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_PM149_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Pmp[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Pmp[n], 1.0, ptr, 0, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_SM149_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Smp[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Smp[n], 1.0, ptr, 0, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_I135_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Ia[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Ia[n], 1.0, ptr, 0, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_XE135_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Xea[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Xea[n], 1.0, ptr, 0, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_PM149_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Pma[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Pma[n], 1.0, ptr, 0, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_SM149_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(flx0[n]*Sma[n], 1.0, ptr, 0, i + 1);
	  AddBuf1D(flx0[n]*Sma[n], 1.0, ptr, 0, 0);
	
	  /* Scattering matrix */

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_SCATTXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  for (m = 0; m < nmg; m++)
	    {
	      j = (long)imap[m];
	      AddBuf1D(flx0[n]*p0[m][n], 1.0, ptr, 0, -1, i, j);

	      /* Tää oli ennen: (14.3.2012) */
	      /*
	      AddBuf1D(flx0[n]*np0[n][m], 1.0, ptr, 0, -1, j, i);
	      */
   	    }

	  /* Scattering production matrix */

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_SCATTPRODXS]; 
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
 
	  for (m = 0; m < nmg; m++)
	    {
	      j = (long)imap[m];
	      AddBuf1D(flx0[n]*np0[m][n], 1.0, ptr, 0, -1, i, j);
   	    }

	  /* Calculate reduced absorption rate */

	  val = flx0[n]*tot[n];
	  for (m = 0; m < nmg; m++)
	    val = val - flx0[n]*np0[m][n];

	  /* Add to statistics */

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_RABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  AddBuf1D(val, 1.0, ptr, 0, i + 1);
	  AddBuf1D(val, 1.0, ptr, 0, 0);
	}
    
      /* k-inf and B2 */

      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_KINF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(kinf, 1.0, ptr, 0, 0);

      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_BUCKLING];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      AddBuf1D(B2, 1.0, ptr, 0, 0);

      /***********************************************************************/

      /* Next universe */

      gcu = NextItem(gcu);

      /***********************************************************************/
    }
  
  /***************************************************************************/

  /***** Add to statistics ***************************************************/

  /* Reduce buffer */

  ReduceBuffer();

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /***********************************************************************/
      
      /***** Data from buffer to  statistics *********************************/

      /* Reset sum for total removal xs */

      sum1 = 0.0;

      /* Calculate few-group cross sections */

      for (i = 0; i < nfg;  i++)
	{
	  /* Reset sum for group removal xs */

	  sum2 = 0.0;

	  /* Few-group flux */

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_FLUX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  flx = BufVal(ptr, i + 1);
	  AddStat(flx, ptr, i + 1);

	  /* Divide reaction rates by flux */

	  if (flx > 0.0)
	    {
	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_TOTXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);

	      sum1 = sum1 + val;
	      sum2 = sum2 + val;
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_ABSXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);

	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_FISSXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);

	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_NSF];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);

	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_DIFFCOEF];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);

	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_RABSXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);

	      ptr = (long)RDB[gcu + GCU_FUM_FG_I135_PROD_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_XE135_PROD_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_PM149_PROD_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_SM149_PROD_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_I135_ABS_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_XE135_ABS_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_SM149_ABS_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	      
	      ptr = (long)RDB[gcu + GCU_FUM_FG_PM149_ABS_XS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i + 1);
	      AddStat(val/flx, ptr, i + 1);
	    	      
	      /* Scattering matrix */

	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_SCATTXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      for (j = 0; j < nfg; j++)
		{
		  val = BufVal(ptr, i, j);
		  AddStat(val/flx, ptr, j, i);
		  
		  sum1 = sum1 - val;

		  if (j == i)
		    sum2 = sum2 - val;
		}

	      /* Scattering production matrix */

	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_SCATTPRODXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      for (j = 0; j < nfg; j++)
		{
		  val = BufVal(ptr, i, j);
		  AddStat(val/flx, ptr, j, i);
		}

	      /* Removal cross section */

	      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_REMXS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddStat(sum2/flx, ptr, i + 1);
	    }
	  
	  /* Chi */
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_CHI];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, i);
	  AddStat(val, ptr, i);
	}

      /* Total flux */

      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      flx = BufVal(ptr, 0);
      AddStat(flx, ptr, 0);

      /* Divide total reaction rates by flux */
      
      if (flx > 0.0)
	{
	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_TOTXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_ABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_FISSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_NSF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_DIFFCOEF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_RABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_B1_REMXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddStat(sum1/flx, ptr, 0);

	  ptr = (long)RDB[gcu + GCU_FUM_FG_I135_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_XE135_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_PM149_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_SM149_PROD_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_I135_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_XE135_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_PM149_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	  
	  ptr = (long)RDB[gcu + GCU_FUM_FG_SM149_ABS_XS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  val = BufVal(ptr, 0);
	  AddStat(val/flx, ptr, 0);
	}

      /* Store k-inf and B2 */

      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_KINF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      kinf = BufVal(ptr, 0);
      AddStat(kinf, ptr, 0);

      ptr = (long)RDB[gcu + GCU_FUM_FG_B1_BUCKLING];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
      B2 = BufVal(ptr, 0);
      AddStat(B2, ptr, 0);
    
      /* Next universe */

      gcu = NextItem(gcu);

      /***********************************************************************/
    }

  /***************************************************************************/

  /* Free allocated memory */

  for(n = 0; n < nmg; n++)
    Mem(MEM_FREE, np0[n]);
  
  Mem(MEM_FREE, np0);
  
  for(n = 0; n < nmg; n++)
    Mem(MEM_FREE, np1[n]);

  Mem(MEM_FREE, np1);

  for(n = 0; n < nmg; n++)
    Mem(MEM_FREE, p0[n]);
  
  Mem(MEM_FREE, p0);

  ccsMatrixFree(Dinv);
  ccsMatrixFree(Dinv2);
  ccsMatrixFree(A);

  Mem(MEM_FREE, ei); 
  Mem(MEM_FREE, di); 
  Mem(MEM_FREE, D); 
  Mem(MEM_FREE, flx1); 
  Mem(MEM_FREE, Db1); 
  Mem(MEM_FREE, chi);
  Mem(MEM_FREE, x);
  Mem(MEM_FREE, alpha);

#endif
}

/*****************************************************************************/

void PrintB1Data(long nmg, double *flx0, complex *flx1, 
		 double *tot, double *abs, double *fiss, 
		 double *nsf, double *mubar, complex *chi, double **p0, 
		 double **np0, double **np1)
{
  long n, m, ptr, i;
  FILE *fp;

  /* Open file */

  fp = fopen("b1_super.m", "w");
  if (fp == NULL)
    Die(FUNCTION_NAME, "Could not open m-file\n"); 

  /* Get pointer to microgroup energy grid */
  
  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	

  ptr = (long)RDB[ptr + ENERGY_GRID_PTR_DATA];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	

  /* Print data */

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "E(%ld) = %e;\n", i++, RDB[ptr + nmg - n]);
      fprintf(fp, "E(%ld) = %e;\n", i++, RDB[ptr + nmg - n - 1]);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "flx0(%ld) = %e;\n", i++, flx0[n]);
      fprintf(fp, "flx0(%ld) = %e;\n", i++, flx0[n]);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "flx1(%ld) = %e;\n", i++, flx1[n].re);
      fprintf(fp, "flx1(%ld) = %e;\n", i++, flx1[n].re);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "tot(%ld) = %e;\n", i++, tot[n]);
      fprintf(fp, "tot(%ld) = %e;\n", i++, tot[n]);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "abs(%ld) = %e;\n", i++, abs[n]);
      fprintf(fp, "abs(%ld) = %e;\n", i++, abs[n]);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "fiss(%ld) = %e;\n", i++, fiss[n]);
      fprintf(fp, "fiss(%ld) = %e;\n", i++, fiss[n]);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "nsf(%ld) = %e;\n", i++, nsf[n]);
      fprintf(fp, "nsf(%ld) = %e;\n", i++, nsf[n]);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "mubar(%ld) = %e;\n", i++, mubar[n]);
      fprintf(fp, "mubar(%ld) = %e;\n", i++, mubar[n]);
    }

  i = 1;
  for (n=0; n < nmg; n++)
    {
      fprintf(fp, "chi(%ld) = %e;\n", i++, chi[n].re);
      fprintf(fp, "chi(%ld) = %e;\n", i++, chi[n].re);
    }

  for (n=0; n < nmg; n++)
    for (m=0; m < nmg; m++)
      fprintf(fp, "p0(%ld, %ld) = %e;\n", n+1, m+1, p0[n][m]);

  for (n=0; n < nmg; n++)
    for (m=0; m < nmg; m++)
      fprintf(fp, "np0(%ld, %ld) = %e;\n", n+1, m+1, np0[n][m]);

  for (n=0; n < nmg; n++)
    for (m=0; m < nmg; m++)
      fprintf(fp, "np1(%ld, %ld) = %e;\n", n+1, m+1, np1[n][m]);
  
  /* Close file */
    
  fclose(fp);
}

/*****************************************************************************/
