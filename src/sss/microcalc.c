#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : microcalc.c                                    */
/*                                                                           */
/* Created:       2013/03/04 (JLe)                                           */
/* Last modified: 2016/03/05 (JLe)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Group constant calculation using homogenized micro-group     */
/*              data                                                         */
/*                                                                           */
/* Comments: - Micro-group cross sections are calculated in                  */
/*             calcmicrogroupxs.c and cleared in clearmicrogroupxs.c         */
/*                                                                           */
/*           - Noita REP_TIME ja PROMPT_LIFE -muuttujia ei tulosteta         */
/*             matlaboutput.c:ssä --> ne voisi ehkä poistaa jos niitä ei     */
/*             oikeasti tarvii mihinkään (voi olla että cax-outputissa       */
/*             tarvitaan).                                                   */
/*                                                                           */
/*           - Calculation of microscopic poison cross sections was changed  */
/*             5.3.2016 / 2.1.26 (JLe)                                       */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "MicroCalc:"

/*****************************************************************************/

void MicroCalc()
{
  long ptr, nmg, nfg, gcu, adf, ppw, alb, n, m, i, j, np, ns, m1, m2, i1, i2;
  double flx, val, div, norm, mubar;
  const double *flx0, *imap, *dfs, *dfc, *JSin, *JSout, *JCin, *JCout;
  const double *JMin, *JMout;
  double *flx1, *D1, src, loss, kinf;
  const double *ds, *dm, *dc;
  const double *tot, *fiss, *abs, *nsf, *fissE, *invv, *chit, *chip, *chid;
  const double *s0, *sp0, *s1, *sp1, *s2, *sp2, *s3, *sp3, *s4, *sp4, *s5, *sp5;
  const double *s6, *sp6, *s7, *sp7, *fflx;
  const double *Xep, *Ip, *Smp, *Pm7p, *Pm8p, *Pm8mp, *Pm9p,*Ia, *Xea;
  const double *Pm7a, *Pm8a, *Pm8ma, *Pm9a, *Sma, *Xeam, *Smam, *ppc, *ppp;

  /* Check if access to RES2 array is allowed */

  if ((long)RDB[DATA_PRIVA_MEM_READY] == NO)
    Die(FUNCTION_NAME, "RES2 array not ready for access");

  /* Get normalization coefficient (JLe: tätä muutettiin 1.7.2015 / 2.1.25) */

  norm = RDB[DATA_MICRO_CALC_NORM]/RDB[DATA_MICRO_CALC_BATCH_SIZE]/
    RDB[DATA_MICRO_CALC_BATCH_SIZE];
  CheckValue(FUNCTION_NAME, "norm", "", norm, 0.0, INFTY);

  /* Mark scoring buffer unreduced (tän pitäisi toimia koska  */
  /* käsitellään sellasia muuttujia joita ei lisätä bufferiin */
  /* muualla). */
  
  WDB[DATA_BUF_REDUCED] = (double)NO;

  /* Get pointer to micro-group structure */
  
  ptr = (long)RDB[DATA_MICRO_PTR_EGRID];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

  /* Number of micro- and macro-groups */

  nmg = (long)RDB[ptr + ENERGY_GRID_NE] - 1;
  nfg = (long)RDB[DATA_ERG_FG_NG];

  /* Index map */

  ptr = (long)RDB[DATA_MICRO_PTR_IDX_MAP];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  imap = &RDB[ptr];

  /* Reduce results buffer */

  WDB[DATA_RES2_REDUCED] = (double)NO;
  ReducePrivateRes();

  /* Mark scoring buffer unreduced (tän pitäisi toimia koska  */
  /* käsitellään sellasia muuttujia joita ei lisätä bufferiin */
  /* muualla). */
  
  WDB[DATA_BUF_REDUCED] = (double)NO;

  /***************************************************************************/

  /***** Infinite and critical spectrum calculations *************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /***********************************************************************/

      /***** Get micro-group data ********************************************/

      /* Reaction cross sections */

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

      ptr = (long)RDB[gcu + GCU_MICRO_FISSE];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      fissE = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_INV_V];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      invv = &RES2[ptr];

      /* Fission spectra */

      ptr = (long)RDB[gcu + GCU_MICRO_CHIT];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      chit = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_CHIP];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      chip = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_CHID];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      chid = &RES2[ptr];
    
      /* Scattering matrixes */

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT0];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s0 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP0];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp0 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT1];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s1 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP1];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp1 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT2];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s2 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP2];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp2 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT3];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s3 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP3];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp3 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT4];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s4 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP4];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp4 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT5];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s5 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP5];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp5 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT6];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s6 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP6];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp6 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATT7];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      s7 = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_SCATTP7];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      sp7 = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_FISS_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      fflx = &RES2[ptr];
      
      /* Check poison calculation */

      if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
	{
	  ptr = (long)RDB[gcu + GCU_MICRO_I135_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Ip = &RES2[ptr];

	  ptr = (long)RDB[gcu + GCU_MICRO_XE135_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Xep = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM147_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm7p = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM148_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm8p = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM148M_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm8mp = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM149_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm9p = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_SM149_YIELD];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Smp = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_I135_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Ia = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_XE135_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Xea = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM147_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm7a = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM148_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm8a = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM148M_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm8ma = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_PM149_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Pm9a = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_SM149_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Sma = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_XE135_MACRO_ABS];
	  CheckPointer(FUNCTION_NAME, "(pt1)", DATA_ARRAY, ptr);
	  Xeam = &RES2[ptr];
	  
	  ptr = (long)RDB[gcu + GCU_MICRO_SM149_MACRO_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  Smam = &RES2[ptr];
	}
      else
	{
	  /* Avoid compiler warning */

	  Ip = NULL;
	  Xep = NULL;
	  Pm7p = NULL; 
	  Pm8p = NULL; 
	  Pm8mp = NULL; 
	  Pm9p = NULL; 
	  Smp = NULL; 
	  Ia = NULL; 
	  Xea = NULL;
	  Pm7a = NULL;
	  Pm8a = NULL;
	  Pm8ma = NULL;
	  Pm9a = NULL;
	  Sma = NULL;
	  Xeam = NULL;
	  Smam = NULL;
	}

      /***********************************************************************/

      /***** Infinite-spectrum calculations **********************************/

      /* Infinite spectrum */

      ptr = (long)RDB[gcu + GCU_MICRO_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      flx0 = &RES2[ptr];

      /* Loop over micro-group structure */

      for (n = 0; n < nmg; n++)
	{
	  /* Get macro-group index */

	  i = (long)imap[n];
	  CheckValue(FUNCTION_NAME, "i", "", i, 0, nfg - 1);

	  /* Flux */

	  ptr = (long)RDB[gcu + GCU_INF_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n], 1.0, ptr, 0, i);
	  
	  /* Flux in fissile zones */

	  ptr = (long)RDB[gcu + GCU_INF_FISS_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*fflx[n], 1.0, ptr, 0, i);

	  /* Total cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_TOT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*tot[n], 1.0, ptr, 0, i);
	  
	  /* Total rate for removal */

	  ptr = (long)RDB[gcu + GCU_INF_REMXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*tot[n], 1.0, ptr, 0, i);

	  /* Fission */

	  ptr = (long)RDB[gcu + GCU_INF_FISS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*fiss[n], 1.0, ptr, 0, i);

	  /* Capture */

	  ptr = (long)RDB[gcu + GCU_INF_CAPT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*(abs[n] - fiss[n]), 1.0, ptr, 0, i);

	  /* Absorption */

	  ptr = (long)RDB[gcu + GCU_INF_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*abs[n], 1.0, ptr, 0, i);

	  /* Absorption rate for reduced absorption */

	  ptr = (long)RDB[gcu + GCU_INF_RABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*abs[n], 1.0, ptr, 0, i);

	  /* Fission neutron production */

	  ptr = (long)RDB[gcu + GCU_INF_NSF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*nsf[n], 1.0, ptr, 0, i);

	  /* Fission energy production */

	  ptr = (long)RDB[gcu + GCU_INF_KAPPA];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*fissE[n], 1.0, ptr, 0, i);

	  /* 1/v */

	  ptr = (long)RDB[gcu + GCU_INF_INVV];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx0[n]*invv[n], 1.0, ptr, 0, i);

	  /* Check poison calculation */

	  if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
	    {
	      /* I-135 production */

	      ptr = (long)RDB[gcu + GCU_INF_I135_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*fiss[n]*Ip[n], 1.0, ptr, 0, i);

	      /* Xe-135 production */

	      ptr = (long)RDB[gcu + GCU_INF_XE135_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*fiss[n]*Xep[n], 1.0, ptr, 0, i);

	      /* Pm production */

	      ptr = (long)RDB[gcu + GCU_INF_PM147_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*fiss[n]*Pm7p[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_INF_PM148_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*fiss[n]*Pm8p[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_INF_PM148M_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*fiss[n]*Pm8mp[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_INF_PM149_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*fiss[n]*Pm9p[n], 1.0, ptr, 0, i);

	      /* Sm-149 production */

	      ptr = (long)RDB[gcu + GCU_INF_SM149_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*fiss[n]*Smp[n], 1.0, ptr, 0, i);

	      /* I-135 absorption */

	      ptr = (long)RDB[gcu + GCU_INF_I135_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx0[n]*fflx[n]*Ia[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx0[n]*Ia[n], 1.0, ptr, 0, i);

	      /* Xe-135 absorption */

	      ptr = (long)RDB[gcu + GCU_INF_XE135_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx0[n]*fflx[n]*Xea[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx0[n]*Xea[n], 1.0, ptr, 0, i);

	      /* Pm absorption */

	      ptr = (long)RDB[gcu + GCU_INF_PM147_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx0[n]*fflx[n]*Pm7a[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx0[n]*Pm7a[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_INF_PM148_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx0[n]*fflx[n]*Pm8a[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx0[n]*Pm8a[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_INF_PM148M_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx0[n]*fflx[n]*Pm8ma[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx0[n]*Pm8ma[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_INF_PM149_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx0[n]*fflx[n]*Pm9a[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx0[n]*Pm9a[n], 1.0, ptr, 0, i);

	      /* Sm-149 absorption */

	      ptr = (long)RDB[gcu + GCU_INF_SM149_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx0[n]*fflx[n]*Sma[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx0[n]*Sma[n], 1.0, ptr, 0, i);

	      /* Macroscopic  Xe-135 absorption */

	      ptr = (long)RDB[gcu + GCU_INF_XE135_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*Xeam[n], 1.0, ptr, 0, i);

	      /* Macroscopic  Sm-149 absorption */

	      ptr = (long)RDB[gcu + GCU_INF_SM149_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx0[n]*Smam[n], 1.0, ptr, 0, i);
	    }
	  
	  /* Calculate total P0 and P1 scattering rates */

	  val = 0.0;
	  for (m = 0; m < nmg; m++)
	    val = val + flx0[n]*s1[n*nmg + m];

	  div = 0.0;
	  for (m = 0; m < nmg; m++)
	    div = div + flx0[n]*s0[n*nmg + m];

	  /* Check */
	
	  if (div > 0.0)
	    {
	      /* Calculate mubar */

	      mubar = val/div;

	      /* Transport cross section */

	      if ((val = tot[n] - mubar*(tot[n] - abs[n])) > 0.0)
		{
		  ptr = (long)RDB[gcu + GCU_INF_TRANSPXS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  AddBuf1D(flx0[n]/val, 1.0, ptr, 0, i);
		}
	    }

	  /* Total fission spectrum */

	  ptr = (long)RDB[gcu + GCU_INF_CHIT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(chit[n], 1.0, ptr, 0, i);

	  /* Pompt fission spectrum */

	  ptr = (long)RDB[gcu + GCU_INF_CHIP];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(chip[n], 1.0, ptr, 0, i);

	  /* Delayed fission spectrum */

	  ptr = (long)RDB[gcu + GCU_INF_CHID];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(chid[n], 1.0, ptr, 0, i);
	}

      /* Loop over micro-group structures */

      for (n = 0; n < nmg; n++)
	for (m = 0; m < nmg; m++)
	  {
	    /* Get macro-group indexes */
	    
	    i = (long)imap[n];
	    j = (long)imap[m];

	    CheckValue(FUNCTION_NAME, "i", "", i, 0, nfg - 1);
	    CheckValue(FUNCTION_NAME, "j", "", j, 0, nfg - 1);

	    /* P0 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s0[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s0[n*nmg + m], 1.0, ptr, 0, i);

	    /* P0 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp0[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATTP0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp0[n*nmg + m], 1.0, ptr, 0, i);

	    /* P1 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s1[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s1[n*nmg + m], 1.0, ptr, 0, i);

	    /* P1 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp1[n*nmg + m], 1.0, ptr, 0, -1, i, j);
	    
	    ptr = (long)RDB[gcu + GCU_INF_SCATTP1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp1[n*nmg + m], 1.0, ptr, 0, i);

	    /* P2 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s2[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s2[n*nmg + m], 1.0, ptr, 0, i);

	    /* P2 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp2[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATTP2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp2[n*nmg + m], 1.0, ptr, 0, i);

	    /* P3 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s3[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s3[n*nmg + m], 1.0, ptr, 0, i);

	    /* P3 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp3[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATTP3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp3[n*nmg + m], 1.0, ptr, 0, i);

	    /* P4 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s4[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s4[n*nmg + m], 1.0, ptr, 0, i);

	    /* P4 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp4[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATTP4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp4[n*nmg + m], 1.0, ptr, 0, i);

	    /* P5 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s5[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s5[n*nmg + m], 1.0, ptr, 0, i);

	    /* P5 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp5[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATTP5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp5[n*nmg + m], 1.0, ptr, 0, i);

	    /* P6 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s6[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s6[n*nmg + m], 1.0, ptr, 0, i);

	    /* P6 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp6[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATTP6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp6[n*nmg + m], 1.0, ptr, 0, i);

	    /* P7 scattering rate */

	    ptr = (long)RDB[gcu + GCU_INF_S7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*s7[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATT7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*s7[n*nmg + m], 1.0, ptr, 0, i);

	    /* P7 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_INF_SP7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx0[n]*sp7[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_INF_SCATTP7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx0[n]*sp7[n*nmg + m], 1.0, ptr, 0, i);

	    /* Reduced absorption cross section */

	    ptr = (long)RDB[gcu + GCU_INF_RABSXS];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	    AddBuf1D(-flx0[n]*(sp0[n*nmg + m] - s0[n*nmg + m]), 1.0, ptr, 0, i);

	    /* Removal cross section */

	    if (j == i)
	      {
		ptr = (long)RDB[gcu + GCU_INF_REMXS];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
		AddBuf1D(-flx0[n]*s0[n*nmg + m], 1.0, ptr, 0, i);
	      }
	  }

      /************************************************************************/

      /***** Additional integral parameters ***********************************/
      
      /* Calculate total source term */

      src = 0.0;
      for (n = 0; n < nmg; n++)
	src = src + flx0[n]*nsf[n];

      /* Calculate total loss term */

      loss = 0.0;
      for (n = 0; n < nmg; n++)
	{
	  /* Add total reaction rate */

	  loss = loss + flx0[n]*tot[n];
      
	  /* Subtract scattering production */
	  
	  for (m = 0; m < nmg; m++)
	    loss = loss - flx0[n]*sp0[n*nmg + m];
	}

      /* Calculate k-inf */
      
      if (src == 0.0)
	kinf = 0.0;
      else if (loss > 0.0)
	kinf = src/loss;
      else
	{
	  /* Zero or negative loss rate indicates very high (n,xn) */
	  /* production, which is a bit strange -- print warning */
	  
	  Warn(FUNCTION_NAME, "Zero or negative loss rate (%E %E)", loss, src);
	  
	  /* Set k-inf to zero */

	  kinf = 0.0;
	}
 
      /* Store values */

      ptr = (long)RDB[gcu + GCU_INF_KINF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
      AddStat(kinf, ptr, 0);
      
      /* Store spectrum */

      ptr = (long)RDB[gcu + GCU_INF_MICRO_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      for (n = 0; n < nmg; n++)
	AddStat(flx0[n], ptr, n);

      /***********************************************************************/
      
      /***** Critical spectrum calculations **********************************/

      /* Check if B1 calculation is invoked */

      if ((long)RDB[DATA_B1_CALC] == NO)
	{
	  /* Not run, Pointer to next universe */
	  
	  gcu = NextItem(gcu);
	  
	  /* Cycle loop */
	  
	  continue;
	}

      /* Array for B1 spectrum */

      ptr = (long)RDB[gcu + GCU_MICRO_B1_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      flx1 = &RES2[ptr];

      /* Array for B1 diffusion coefficients */

      ptr = (long)RDB[gcu + GCU_MICRO_B1_DIFFCOEF];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      D1 = &RES2[ptr];
      
      /* Solve B1 flux */
      
      if (B1Flux(gcu, nmg, flx0, tot, abs, nsf, sp0, sp1, chit, flx1, D1) == NO)
	{
	  /* B1 calculation failed, copy infinite spectrum to B1 spectrum */

	  memcpy(flx1, flx0, nmg*sizeof(double));
   
	  /* Calculate micro-group wise diffusion coefficient */

	  for (n = 0; n < nmg; n++)
	    {
	      /* Calculate total P0 and P1 scattering rates */

	      val = 0.0;
	      for (m = 0; m < nmg; m++)
		val = val + flx1[n]*s1[n*nmg + m];
	      
	      div = 0.0;
	      for (m = 0; m < nmg; m++)
		div = div + flx1[n]*s0[n*nmg + m];
	      
	      /* Check */
	      
	      if (div > 0.0)
		{
		  /* Calculate mubar */
		  
		  mubar = val/div;
		  
		  /* Transport cross section */

		  if ((val = tot[n] - mubar*(tot[n] - abs[n])) > 0.0)
		    {
		      /* Diffusion coefficient */
		      
		      D1[n] = 1.0/(3.0*val);
		    }
		}
	    }
	}

      /* Loop over micro-group structure */

      for (n = 0; n < nmg; n++)
	{
	  /* Get macro-group index */
	  
	  i = (long)imap[n];
	  CheckValue(FUNCTION_NAME, "i", "", i, 0, nfg - 1);

	  /* Flux */

	  ptr = (long)RDB[gcu + GCU_B1_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n], 1.0, ptr, 0, i);

	  /* Flux in fissile zones */

	  ptr = (long)RDB[gcu + GCU_B1_FISS_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*fflx[n], 1.0, ptr, 0, i);

	  /* Total cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_TOT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*tot[n], 1.0, ptr, 0, i);
	  
	  /* Total rate for removal */

	  ptr = (long)RDB[gcu + GCU_B1_REMXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*tot[n], 1.0, ptr, 0, i);

	  /* Fission */

	  ptr = (long)RDB[gcu + GCU_B1_FISS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*fiss[n], 1.0, ptr, 0, i);

	  /* Capture */

	  ptr = (long)RDB[gcu + GCU_B1_CAPT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*(abs[n] - fiss[n]), 1.0, ptr, 0, i);

	  /* Absorption */

	  ptr = (long)RDB[gcu + GCU_B1_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*abs[n], 1.0, ptr, 0, i);

	  /* Absorption rate for reduced absorption */

	  ptr = (long)RDB[gcu + GCU_B1_RABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*abs[n], 1.0, ptr, 0, i);

	  /* Fission neutron production */

	  ptr = (long)RDB[gcu + GCU_B1_NSF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*nsf[n], 1.0, ptr, 0, i);

	  /* Fission energy production */

	  ptr = (long)RDB[gcu + GCU_B1_KAPPA];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*fissE[n], 1.0, ptr, 0, i);

	  /* 1/v */

	  ptr = (long)RDB[gcu + GCU_B1_INVV];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*invv[n], 1.0, ptr, 0, i);

	  /* Diffusion coefficient */

	  ptr = (long)RDB[gcu + GCU_B1_DIFFCOEF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddBuf1D(flx1[n]*D1[n], 1.0, ptr, 0, i);

	  /* Check poison calculation */

	  if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
	    {
	      /* I-135 production */

	      ptr = (long)RDB[gcu + GCU_B1_I135_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*fiss[n]*Ip[n], 1.0, ptr, 0, i);

	      /* Xe-135 production */

	      ptr = (long)RDB[gcu + GCU_B1_XE135_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*fiss[n]*Xep[n], 1.0, ptr, 0, i);

	      /* Pm production */

	      ptr = (long)RDB[gcu + GCU_B1_PM147_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*fiss[n]*Pm7p[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_B1_PM148M_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*fiss[n]*Pm8p[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_B1_PM148_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*fiss[n]*Pm8mp[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_B1_PM149_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*fiss[n]*Pm9p[n], 1.0, ptr, 0, i);

	      /* Sm-149 production */

	      ptr = (long)RDB[gcu + GCU_B1_SM149_YIELD];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*fiss[n]*Smp[n], 1.0, ptr, 0, i);

	      /* I-135 absorption */

	      ptr = (long)RDB[gcu + GCU_B1_I135_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx1[n]*fflx[n]*Ia[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx1[n]*fflx[n]*Ia[n], 1.0, ptr, 0, i);

	      /* Xe-135 absorption */

	      ptr = (long)RDB[gcu + GCU_B1_XE135_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx1[n]*fflx[n]*Xea[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx1[n]*Xea[n], 1.0, ptr, 0, i);

	      /* Pm absorption */

	      ptr = (long)RDB[gcu + GCU_B1_PM147_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx1[n]*fflx[n]*Pm7a[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx1[n]*Pm7a[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_B1_PM148_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx1[n]*fflx[n]*Pm8a[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx1[n]*Pm8a[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_B1_PM148M_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx1[n]*fflx[n]*Pm8ma[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx1[n]*Pm8ma[n], 1.0, ptr, 0, i);

	      ptr = (long)RDB[gcu + GCU_B1_PM149_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx1[n]*fflx[n]*Pm9a[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx1[n]*Pm9a[n], 1.0, ptr, 0, i);

	      /* Sm-149 absorption */

	      ptr = (long)RDB[gcu + GCU_B1_SM149_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      /*
	      AddBuf1D(flx1[n]*fflx[n]*Sma[n], 1.0, ptr, 0, i);
	      */
	      AddBuf1D(flx1[n]*Sma[n], 1.0, ptr, 0, i);

	      /* Macroscopic  Xe-135 absorption */

	      ptr = (long)RDB[gcu + GCU_B1_XE135_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*Xeam[n], 1.0, ptr, 0, i);

	      /* Macroscopic  Sm-149 absorption */

	      ptr = (long)RDB[gcu + GCU_B1_SM149_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      AddBuf1D(flx1[n]*Smam[n], 1.0, ptr, 0, i);
	    }

	  /* Calculate total P0 and P1 scattering rates (B1-transport- */
	  /* vaikutusalan laskemista muutettiin 21.8.2014 / 2.1.22.) */ 
	  
	  val = 0.0;
	  for (m = 0; m < nmg; m++)
	    val = val + flx0[n]*s1[n*nmg + m];

	  div = 0.0;
	  for (m = 0; m < nmg; m++)
	    div = div + flx0[n]*s0[n*nmg + m];

	  /* Check */
	
	  if (div > 0.0)
	    {
	      /* Calculate mubar */

	      mubar = val/div;

	      /* Transport cross section */

	      if ((val = tot[n] - mubar*(tot[n] - abs[n])) > 0.0)
		{
		  ptr = (long)RDB[gcu + GCU_B1_TRANSPXS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  AddBuf1D(flx1[n]/val, 1.0, ptr, 0, i);
		}
	    }

	  /* Total fission spectrum */

	  ptr = (long)RDB[gcu + GCU_B1_CHIT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(chit[n], 1.0, ptr, 0, i);

	  /* Pompt fission spectrum */

	  ptr = (long)RDB[gcu + GCU_B1_CHIP];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(chip[n], 1.0, ptr, 0, i);

	  /* Delayed fission spectrum */

	  ptr = (long)RDB[gcu + GCU_B1_CHID];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  AddBuf1D(chid[n], 1.0, ptr, 0, i);
	}

      /* Loop over micro-group structures */

      for (n = 0; n < nmg; n++)
	for (m = 0; m < nmg; m++)
	  {
	    /* Get macro-group indexes */
	    
	    i = (long)imap[n];
	    j = (long)imap[m];

	    CheckValue(FUNCTION_NAME, "i", "", i, 0, nfg - 1);
	    CheckValue(FUNCTION_NAME, "j", "", j, 0, nfg - 1);

	    /* P0 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s0[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s0[n*nmg + m], 1.0, ptr, 0, i);

	    /* P0 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp0[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATTP0];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp0[n*nmg + m], 1.0, ptr, 0, i);

	    /* P1 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s1[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s1[n*nmg + m], 1.0, ptr, 0, i);

	    /* P1 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp1[n*nmg + m], 1.0, ptr, 0, -1, i, j);
	    
	    ptr = (long)RDB[gcu + GCU_B1_SCATTP1];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp1[n*nmg + m], 1.0, ptr, 0, i);

	    /* P2 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s2[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s2[n*nmg + m], 1.0, ptr, 0, i);

	    /* P2 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp2[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATTP2];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp2[n*nmg + m], 1.0, ptr, 0, i);

	    /* P3 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s3[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s3[n*nmg + m], 1.0, ptr, 0, i);

	    /* P3 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp3[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATTP3];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp3[n*nmg + m], 1.0, ptr, 0, i);

	    /* P4 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s4[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s4[n*nmg + m], 1.0, ptr, 0, i);

	    /* P4 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp4[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATTP4];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp4[n*nmg + m], 1.0, ptr, 0, i);

	    /* P5 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s5[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s5[n*nmg + m], 1.0, ptr, 0, i);

	    /* P5 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp5[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATTP5];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp5[n*nmg + m], 1.0, ptr, 0, i);

	    /* P6 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s6[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s6[n*nmg + m], 1.0, ptr, 0, i);

	    /* P6 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp6[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATTP6];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp6[n*nmg + m], 1.0, ptr, 0, i);

	    /* P7 scattering rate */

	    ptr = (long)RDB[gcu + GCU_B1_S7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*s7[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATT7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*s7[n*nmg + m], 1.0, ptr, 0, i);

	    /* P7 scattering production rate */

	    ptr = (long)RDB[gcu + GCU_B1_SP7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf(flx1[n]*sp7[n*nmg + m], 1.0, ptr, 0, -1, i, j);

	    ptr = (long)RDB[gcu + GCU_B1_SCATTP7];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	    AddBuf1D(flx1[n]*sp7[n*nmg + m], 1.0, ptr, 0, i);

	    /* Reduced absorption cross section */

	    ptr = (long)RDB[gcu + GCU_B1_RABSXS];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	    AddBuf1D(-flx1[n]*(sp0[n*nmg + m] - s0[n*nmg + m]), 1.0, ptr, 0, i);

	    /* Removal cross section */

	    if (j == i)
	      {
		ptr = (long)RDB[gcu + GCU_B1_REMXS];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
		AddBuf1D(-flx1[n]*s0[n*nmg + m], 1.0, ptr, 0, i);
	      }
	  }
    
      /************************************************************************/

      /***** Additional integral parameters ***********************************/

      /* Calculate total source term */

      src = 0.0;
      for (n = 0; n < nmg; n++)
	src = src + flx1[n]*nsf[n];

      /* Calculate total loss term */

      loss = 0.0;
      for (n = 0; n < nmg; n++)
	{
	  /* Add total reaction rate */

	  loss = loss + flx1[n]*tot[n];
      
	  /* Subtract scattering production */
	  
	  for (m = 0; m < nmg; m++)
	    loss = loss - flx1[n]*sp0[n*nmg + m];
	}

      /* Calculate k-inf */
      
      if (src == 0.0)
	kinf = 0.0;
      else if (loss > 0.0)
	kinf = src/loss;
      else
	{	  
	  /* Set k-inf to zero (warning message was printed above) */

	  kinf = 0.0;
	}

      /* Store value */

      ptr = (long)RDB[gcu + GCU_B1_KINF];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
      AddStat(kinf, ptr, 0);
    
      /* Store spectrum */

      ptr = (long)RDB[gcu + GCU_B1_MICRO_FLX];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

      for (n = 0; n < nmg; n++)
	AddStat(flx1[n], ptr, n);

      /***********************************************************************/

      /* Next universe */

      gcu = NextItem(gcu);

      /************************************************************************/
    }
  
  /****************************************************************************/

  /***** Discontinuity factors ************************************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check pointer */

      if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) < VALID_PTR)
	{
	  /* Next universe */

	  gcu = NextItem(gcu);

	  /* Cycle loop */

	  continue;
	}

      /***********************************************************************/
      
      /***** Get micro-group data ********************************************/
      
      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      dfs = &RES2[ptr];
      
      if ((ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_FLUX]) > VALID_PTR)
	dfc = &RES2[ptr];
      else
	dfc = NULL;
      
      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      JSin = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_ADF_SURF_OUT_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      JSout = &RES2[ptr];

      ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      JMin = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_ADF_MID_OUT_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      JMout = &RES2[ptr];

      if ((ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_IN_CURR]) > VALID_PTR)
	JCin = &RES2[ptr];
      else
	JCin = NULL;

      if ((ptr = (long)RDB[gcu + GCU_MICRO_ADF_CORN_OUT_CURR]) > VALID_PTR)
	JCout = &RES2[ptr];
      else
	JCout = NULL;

      ptr = (long)RDB[gcu + GCU_MICRO_ADF_CELL_FLUX];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      flx0 = &RES2[ptr];

      /***********************************************************************/

      /***** Fluxes and currents for discontinuity factors *******************/

      /* Loop over micro-group structure */
      
      for (n = 0; n < nmg; n++)
	{
	  /* Get few-group index */
	  
	  i = (long)imap[n];
	  CheckValue(FUNCTION_NAME, "i", "", i, 0, nfg - 1);

	  /* Loop over vertices */

	  for (m = 0; m < (long)RDB[adf + ADF_NSURF]; m++)
	    {
	      /* Surface flux */

	      ptr = RDB[gcu + GCU_RES_FG_DF_HET_SURF_FLUX];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(dfs[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Inward surface averaged current */

	      ptr = RDB[gcu + GCU_RES_FG_DF_SURF_IN_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JSin[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Outward surface averaged current */

	      ptr = RDB[gcu + GCU_RES_FG_DF_SURF_OUT_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JSout[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Net surface averaged current */
	      
	      ptr = RDB[gcu + GCU_RES_FG_DF_SURF_NET_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JSin[n + m*nmg] - JSout[n + m*nmg], 1.0, ptr, 0, 
		     -1, m, i);

	      /* Inward surface mid-point current */

	      ptr = RDB[gcu + GCU_RES_FG_DF_MID_IN_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JMin[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Outward surface mid-point current */

	      ptr = RDB[gcu + GCU_RES_FG_DF_MID_OUT_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JMout[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Net surface mid-point current */
	      
	      ptr = RDB[gcu + GCU_RES_FG_DF_MID_NET_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JMin[n + m*nmg] - JMout[n + m*nmg], 1.0, ptr, 0, 
		     -1, m, i);
	    }

	  /* Loop over corners */

	  for (m = 0; m < (long)RDB[adf + ADF_NCORN]; m++)
	    {
	      /* Corner flux */

	      ptr = RDB[gcu + GCU_RES_FG_DF_HET_CORN_FLUX];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	      
	      AddBuf(dfc[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Inward corner current */

	      ptr = RDB[gcu + GCU_RES_FG_DF_CORN_IN_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JCin[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Outward corner current */

	      ptr = RDB[gcu + GCU_RES_FG_DF_CORN_OUT_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JCout[n + m*nmg], 1.0, ptr, 0, -1, m, i);

	      /* Net corner current */
	      
	      ptr = RDB[gcu + GCU_RES_FG_DF_CORN_NET_CURR];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddBuf(JCin[n + m*nmg] - JCout[n + m*nmg], 1.0, ptr, 0, 
		     -1, m, i);
	    }

	  /* Cell flux */
	  
	  ptr = RDB[gcu + GCU_RES_FG_DF_HET_VOL_FLUX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	      
	  AddBuf1D(flx0[n]/RDB[adf + ADF_VOL], 1.0, ptr, 0, i);
	}

      /***********************************************************************/
    
      /* Next universe */

      gcu = NextItem(gcu);
    }

  /* For EDo (18.1.2014) */

  DiffCoefED(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  /****************************************************************************/

  /***** Coordinates for pin-wise homogeneous flux ****************************/

 /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check pointer */

      if ((ppw = (long)RDB[gcu + GCU_PTR_PPW]) > VALID_PTR)
	{
	  /* Power array */

	  ptr = (long)RDB[gcu + GCU_MICRO_PPW_POW];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  ppp = &RES2[ptr];

	  /* Coordinate array */

	  ptr = (long)RDB[gcu + GCU_MICRO_PPW_XYZ];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  ppc = &RES2[ptr];

	  /* Get number of pins */

	  np = (long)RDB[ppw + PPW_NP];

	  /* Loop over energy groups and pins */
	  
	  for (m = 0; m < np; m++)
	    for (i = 0; i < nfg; i++)
	      {	     
		/* Score power */
		
		ptr = (long)RDB[gcu + GCU_RES_FG_PPW_POW];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		AddBuf(ppp[m*(nfg + 1) + i], 1.0, ptr, 0, -1, m, i);

		/* Pointer to coordinates */
		
		ptr = (long)RDB[gcu + GCU_RES_FG_PPW_XYZ];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

		/* Try group-wise */
		
		if ((div = ppp[m*(nfg + 1) + i]) > 0.0)
		  {
		    AddBuf(ppc[0*np*(nfg + 1) + m*(nfg + 1) + i]/div, 1.0, 
			   ptr, 0, -1, 0, m, i);
		    AddBuf(ppc[1*np*(nfg + 1) + m*(nfg + 1) + i]/div, 1.0, 
			   ptr, 0, -1, 1, m, i);
		    AddBuf(ppc[2*np*(nfg + 1) + m*(nfg + 1) + i]/div, 1.0, 
			   ptr, 0, -1, 2, m, i);
		  }
		
		/* Try total */
		
		else if ((div = ppp[m*(nfg + 1) + nfg]) > 0.0)
		  {
		    AddBuf(ppc[0*np*(nfg + 1) + m*(nfg + 1) + nfg]/div, 1.0, 
			   ptr, 0, -1, 0, m, i);
		    AddBuf(ppc[1*np*(nfg + 1) + m*(nfg + 1) + nfg]/div, 1.0, 
			   ptr, 0, -1, 1, m, i);
		    AddBuf(ppc[2*np*(nfg + 1) + m*(nfg + 1) + nfg]/div, 1.0, 
			   ptr, 0, -1, 2, m, i);
		  }
		else
		  {
		    /* Set value to -1E+18 to indicate zero scores */
		    /* (this is tested in homoflux.c). */

		    AddBuf(-1E+18, 1.0, ptr, 0, -1, 0, m, i);
		    AddBuf(-1E+18, 1.0, ptr, 0, -1, 1, m, i);
		    AddBuf(-1E+18, 1.0, ptr, 0, -1, 2, m, i);
		  }
	      }
	}

      /* Next universe */
      
      gcu = NextItem(gcu);
    }

  /****************************************************************************/

  /***** Currents for albedos and partial albedos *****************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check pointer */

      if ((alb = (long)RDB[gcu + GCU_PTR_ALB]) < VALID_PTR)
	{
	  /* Next universe */

	  gcu = NextItem(gcu);

	  /* Cycle loop */

	  continue;
	}

      /* Get current vectors */
      
      ptr = (long)RDB[gcu + GCU_MICRO_ALB_IN_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      JSin = &RES2[ptr];
      
      ptr = (long)RDB[gcu + GCU_MICRO_ALB_OUT_CURR];
      CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
      JSout = &RES2[ptr];

      /* Number of surfaces */

      ns = (long)RDB[alb + ALB_NSURF];
      CheckValue(FUNCTION_NAME, "ns", "", ns, 1, 8);

      /***********************************************************************/

      /***** Total albedos ***************************************************/

      /* Loop over incoming macro-group structure */

      for (i1 = 0; i1 < nfg; i1++)
	{
	  /* Calculate total inward current */

	  ptr = (long)RDB[gcu + GCU_MICRO_ALB_IN_CURR];
	  CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	  
	  div = 0.0;
	  for (m1 = 0; m1 < ns; m1++)
	    div = div + JSin[m1*nfg + i1];

	  /* Score albedos */

	  if (div > 0.0)
	    {
	      /* Loop over energy groups */

	      for (i2 = 0; i2 < nfg; i2++)
		{
		  val = 0.0;

		  /* Loop over surfaces */
	      
		  for (m1 = 0; m1 < ns; m1++)
		    for (m2 = 0; m2 < ns; m2++)
		      val = val + 
			JSout[m1*nfg*nfg*ns + i1*nfg*ns + m2*nfg + i2];

		  /* Total albedo */

		  ptr = (long)RDB[gcu + GCU_RES_FG_TOT_ALB];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  AddBuf(val/div, 1.0, ptr, 0, -1, i1, i2);		  
		}
	    }
	}

      /***********************************************************************/

      /***** Partial albedos ************************************************/

      /* Loop over incoming macro-group structure and surface */

      for (i1 = 0; i1 < nfg; i1++)
	for (m1 = 0; m1 < ns; m1++)
	  {
	    /* Get inward current */

	    ptr = (long)RDB[gcu + GCU_MICRO_ALB_IN_CURR];
	    CheckPointer(FUNCTION_NAME, "(ptr)", RES2_ARRAY, ptr);
	    div = JSin[m1*nfg + i1];

	    /* Score */

	    ptr = (long)RDB[gcu + GCU_RES_FG_ALB_IN_CURR];
	    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	    AddBuf(div, 1.0, ptr, 0, -1, m1, i1);

	    /* Score albedos */

	    if (div > 0.0)
	      {
		/* Loop over surfaces and energy groups */
		
		for (m2 = 0; m2 < ns; m2++)
		  for (i2 = 0; i2 < nfg; i2++)
		    {
		      val = JSout[m1*nfg*nfg*ns + i1*nfg*ns + m2*nfg + i2];

		      /* Outward current */
		
		      ptr = (long)RDB[gcu + GCU_RES_FG_ALB_OUT_CURR];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      AddBuf(val, 1.0, ptr, 0, -1, m1, i1, m2, i2);

		      /* Partial albedo */
		
		      ptr = (long)RDB[gcu + GCU_RES_FG_PART_ALB];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      AddBuf(val/div, 1.0, ptr, 0, -1, m1, i1, m2, i2);
		    }
	      }
	  }

      /***********************************************************************/

      /* Next universe */

      gcu = NextItem(gcu);
    }

  /****************************************************************************/

  /***** Macro-group calculations *********************************************/

  /* Reduce buffer */

  ReduceBuffer();

  /* 16.5.2014 */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      ApplyGCSymmetries(gcu);
      HomoFlux(gcu);
      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /***** Infinite spectrum calculations **************************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Loop over macro-group structure */

      for (i = 0; i < nfg; i++)
	{
	  /*******************************************************************/

	  /***** Few-group reaction cross sections ***************************/
	  
	  /* Flux */

	  ptr = (long)RDB[gcu + GCU_INF_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  flx = BufVal(ptr, i);
	  AddStat(norm*flx, ptr, i);

	  /* Check value */

	  if (flx < ZERO)
	    continue;

	  /* Total cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_TOT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Fission */
	  
	  ptr = (long)RDB[gcu + GCU_INF_FISS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Check */

	  if ((div = val) > 0.0)
	    {
	      /* Fission neutron production */
	      
	      ptr = (long)RDB[gcu + GCU_INF_NSF];;
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i);
	      AddStat(val/flx, ptr, i);

	      /* Fission nubar */
	      
	      ptr = (long)RDB[gcu + GCU_INF_NUBAR];;
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddStat(val/div, ptr, i);

	      /* Fission kappa */
	      
	      ptr = (long)RDB[gcu + GCU_INF_KAPPA];;
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i)/MEV;
	      AddStat(val/div, ptr, i);
	    }
	  
	  /* Capture */
	  
	  ptr = (long)RDB[gcu + GCU_INF_CAPT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Absorption */
	  
	  ptr = (long)RDB[gcu + GCU_INF_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* 1/v */
	  
	  ptr = (long)RDB[gcu + GCU_INF_INVV];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Reduced absorption */
	  
	  ptr = (long)RDB[gcu + GCU_INF_RABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Removal */
	  
	  ptr = (long)RDB[gcu + GCU_INF_REMXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Transport cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_TRANSPXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  if (val > 0.0)
	    AddStat(flx/val, ptr, i);

	  /* Diffusion coefficient */

	  ptr = (long)RDB[gcu + GCU_INF_DIFFCOEF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  AddStat(val/flx/3.0, ptr, i);

	  /* Total fission spectrum */
	  
	  ptr = (long)RDB[gcu + GCU_INF_CHIT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  

	  div = 0.0;
	  for (n = 0; n < nfg; n++)
	    div = div + BufVal(ptr, n);

	  if (div > 0.0)
	    {
	      val = BufVal(ptr, i);
	      AddStat(val/div, ptr, i);
	    }

	  /* Prompt fission spectrum */
	  
	  ptr = (long)RDB[gcu + GCU_INF_CHIP];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  

	  div = 0.0;
	  for (n = 0; n < nfg; n++)
	    div = div + BufVal(ptr, n);

	  if (div > 0.0)
	    {
	      val = BufVal(ptr, i);
	      AddStat(val/div, ptr, i);
	    }

	  /* Delayed fission spectrum */
	  
	  ptr = (long)RDB[gcu + GCU_INF_CHID];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  

	  div = 0.0;
	  for (n = 0; n < nfg; n++)
	    div = div + BufVal(ptr, n);

	  if (div > 0.0)
	    {
	      val = BufVal(ptr, i);
	      AddStat(val/div, ptr, i);
	    }

	  /*******************************************************************/
	  
	  /***** Poison cross sections ***************************************/

	  /* Check poison calculation */
	  
	  if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
	    {
	      /* Get fission rate */
	  
	      ptr = (long)RDB[gcu + GCU_INF_FISS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      div = BufVal(ptr, i);
	      
	      /* Check */

	      if (div > 0.0)
		{
		  /* I-135 yield */
		  
		  ptr = (long)RDB[gcu + GCU_INF_I135_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr); 
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Xe-135 yield */
		  
		  ptr = (long)RDB[gcu + GCU_INF_XE135_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Pm yield */
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM147_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM148_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM148M_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM149_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Sm-149 yield */
		  
		  ptr = (long)RDB[gcu + GCU_INF_SM149_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		}
	      
	      /* Get flux in fissile zones */
	  
	      ptr = (long)RDB[gcu + GCU_INF_FISS_FLX];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      val = BufVal(ptr, i);
	      AddStat(norm*val, ptr, i);

	      /* Divide with total flux and volume ratio */
	      /* if provided (2.1.26) */

	      if (RDB[DATA_POISON_XS_VOL_RAT] > 0.0)
		div = flx*RDB[DATA_POISON_XS_VOL_RAT];
	      else
		div = flx;

	      /* Check */

	      if (div > 0.0)
		{
		  /* I-135 absorption */
	      
		  ptr = (long)RDB[gcu + GCU_INF_I135_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Xe-135 absorption */
		  
		  ptr = (long)RDB[gcu + GCU_INF_XE135_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Pm absorption */
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM147_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM148_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM148M_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_INF_PM149_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Sm-149 absorption */
		  
		  ptr = (long)RDB[gcu + GCU_INF_SM149_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		}

	      /* Macroscopic Xe-135 absorption */
	      
	      ptr = (long)RDB[gcu + GCU_INF_XE135_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i);
	      AddStat(val/flx, ptr, i);
	      
	      /* Macroscopic Sm-149 absorption */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SM149_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i);
	      AddStat(val/flx, ptr, i);
	    }
	  
	  /*******************************************************************/

	  /***** Total scattering cross sections *****************************/

	  /* P0 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT0];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P0 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP0];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P1 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT1];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P1 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP1];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* P2 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT2];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P2 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP2];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P3 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT3];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P3 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP3];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P4 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT4];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P4 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP4];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P5 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT5];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P5 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP5];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P6 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT6];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P6 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP6];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P7 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATT7];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P7 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_INF_SCATTP7];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /*******************************************************************/
	  
	  /***** Scattering matrixes *****************************************/

	  /* Loop over groups */

	  for (j = 0; j < nfg; j++)
	    {
	      /* P0 scattering cross section */

	      ptr = (long)RDB[gcu + GCU_INF_S0];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P0 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP0];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P1 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_S1];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P1 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP1];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P2 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_S2];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P2 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP2];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P3 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_S3];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P3 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP3];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P4 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_S4];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P4 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP4];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P5 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_S5];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P5 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP5];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P6 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_S6];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P6 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP6];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P7 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_S7];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P7 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_INF_SP7];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	    }
	}
      
      /************************************************************************/

      /* Next universe */

      gcu = NextItem(gcu);

      /************************************************************************/
    }

  /***************************************************************************/

  /***** Critical spectrum calculations **************************************/

  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check convergence */

      if ((long)RDB[gcu + GCU_B1_CONV] == NO)
	{
      	  /* Not run or no convergence, get pointer to next universe */
	  
	  gcu = NextItem(gcu);
	  
	  /* Cycle loop */
	  
	  continue;
	}

      /* Loop over macro-group structure */

      for (i = 0; i < nfg; i++)
	{
	  /*******************************************************************/

	  /***** Few-group reaction cross sections ***************************/
	  
	  /* Flux */

	  ptr = (long)RDB[gcu + GCU_B1_FLX];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  flx = BufVal(ptr, i);
	  AddStat(norm*flx, ptr, i);

	  /* Check value */

	  if (flx < ZERO)
	    continue;

	  /* Total cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_TOT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Fission */
	  
	  ptr = (long)RDB[gcu + GCU_B1_FISS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Check */

	  if ((div = val) > 0.0)
	    {
	      /* Fission neutron production */
	      
	      ptr = (long)RDB[gcu + GCU_B1_NSF];;
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i);
	      AddStat(val/flx, ptr, i);

	      /* Fission nubar */
	      
	      ptr = (long)RDB[gcu + GCU_B1_NUBAR];;
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      AddStat(val/div, ptr, i);

	      /* Fission kappa */
	      
	      ptr = (long)RDB[gcu + GCU_B1_KAPPA];;
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i)/MEV;
	      AddStat(val/div, ptr, i);
	    }
	  
	  /* Capture */
	  
	  ptr = (long)RDB[gcu + GCU_B1_CAPT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Absorption */
	  
	  ptr = (long)RDB[gcu + GCU_B1_ABS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* 1/v */
	  
	  ptr = (long)RDB[gcu + GCU_B1_INVV];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Reduced absorption */
	  
	  ptr = (long)RDB[gcu + GCU_B1_RABSXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Removal */
	  
	  ptr = (long)RDB[gcu + GCU_B1_REMXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Diffusion coefficient */

	  ptr = (long)RDB[gcu + GCU_B1_DIFFCOEF];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* Transport cross section (muutettu 21.8.2014 / 2.1.22) */

	  /*
	  ptr = (long)RDB[gcu + GCU_B1_TRANSPXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  if (val > 0.0)
	    AddStat(flx/val/3.0, ptr, i);
	  */

	  ptr = (long)RDB[gcu + GCU_B1_TRANSPXS];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	  val = BufVal(ptr, i);
	  if (val > 0.0)
	    AddStat(flx/val, ptr, i);

	  /* Total fission spectrum */
	  
	  ptr = (long)RDB[gcu + GCU_B1_CHIT];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  

	  div = 0.0;
	  for (n = 0; n < nfg; n++)
	    div = div + BufVal(ptr, n);

	  if (div > 0.0)
	    {
	      val = BufVal(ptr, i);
	      AddStat(val/div, ptr, i);
	    }

	  /* Prompt fission spectrum */
	  
	  ptr = (long)RDB[gcu + GCU_B1_CHIP];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  

	  div = 0.0;
	  for (n = 0; n < nfg; n++)
	    div = div + BufVal(ptr, n);

	  if (div > 0.0)
	    {
	      val = BufVal(ptr, i);
	      AddStat(val/div, ptr, i);
	    }

	  /* Delayed fission spectrum */
	  
	  ptr = (long)RDB[gcu + GCU_B1_CHID];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  

	  div = 0.0;
	  for (n = 0; n < nfg; n++)
	    div = div + BufVal(ptr, n);

	  if (div > 0.0)
	    {
	      val = BufVal(ptr, i);
	      AddStat(val/div, ptr, i);
	    }

	  /*******************************************************************/
	  
	  /***** Poison cross sections ***************************************/

	  /* Check poison calculation */
	  
	  if ((long)RDB[DATA_OPTI_POISON_CALC] == YES)
	    {
	      /* Get fission rate */
	  
	      ptr = (long)RDB[gcu + GCU_B1_FISS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      div = BufVal(ptr, i);
	      
	      /* Check */

	      if (div > 0.0)
		{
		  /* I-135 yield */
	      
		  ptr = (long)RDB[gcu + GCU_B1_I135_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr); 
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Xe-135 yield */
		  
		  ptr = (long)RDB[gcu + GCU_B1_XE135_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Pm yield */
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM147_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM148_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM148M_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM149_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Sm-149 yield */
		  
		  ptr = (long)RDB[gcu + GCU_B1_SM149_YIELD];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		}

	      /* Get flux in fissile zones */
	  
	      ptr = (long)RDB[gcu + GCU_B1_FISS_FLX];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);	  
	      val = BufVal(ptr, i);
	      AddStat(norm*val, ptr, i);

	      /* Divide with total flux and volume ratio */
	      /* if provided (2.1.26) */

	      if (RDB[DATA_POISON_XS_VOL_RAT] > 0.0)
		div = flx*RDB[DATA_POISON_XS_VOL_RAT];
	      else
		div = flx;

	      /* Check */

	      if (div > 0.0)
		{	      
		  /* I-135 absorption */
		  
		  ptr = (long)RDB[gcu + GCU_B1_I135_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Xe-135 absorption */
		  
		  ptr = (long)RDB[gcu + GCU_B1_XE135_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Pm absorption */
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM147_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM148_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM148M_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  ptr = (long)RDB[gcu + GCU_B1_PM149_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		  
		  /* Sm-149 absorption */
		  
		  ptr = (long)RDB[gcu + GCU_B1_SM149_ABS];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, i);
		  AddStat(val/div, ptr, i);
		}

	      /* Macroscopic Xe-135 absorption */
	      
	      ptr = (long)RDB[gcu + GCU_B1_XE135_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i);
	      AddStat(val/flx, ptr, i);
	      
	      /* Macroscopic Sm-149 absorption */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SM149_MACRO_ABS];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      val = BufVal(ptr, i);
	      AddStat(val/flx, ptr, i);
	    }
	  
	  /*******************************************************************/

	  /***** Total scattering cross sections *****************************/

	  /* P0 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT0];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P0 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP0];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P1 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT1];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P1 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP1];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);

	  /* P2 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT2];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P2 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP2];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P3 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT3];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P3 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP3];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P4 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT4];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P4 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP4];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P5 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT5];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P5 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP5];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P6 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT6];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P6 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP6];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P7 scattering cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATT7];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /* P7 scattering production cross section */
	  
	  ptr = (long)RDB[gcu + GCU_B1_SCATTP7];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	  val = BufVal(ptr, i);
	  AddStat(val/flx, ptr, i);
	  
	  /*******************************************************************/
	  
	  /***** Scattering matrixes *****************************************/

	  /* Loop over groups */

	  for (j = 0; j < nfg; j++)
	    {
	      /* P0 scattering cross section */

	      ptr = (long)RDB[gcu + GCU_B1_S0];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P0 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP0];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P1 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_S1];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P1 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP1];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P2 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_S2];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P2 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP2];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P3 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_S3];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P3 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP3];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);

	      /* P4 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_S4];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P4 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP4];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P5 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_S5];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P5 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP5];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P6 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_S6];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P6 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP6];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P7 scattering cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_S7];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	      
	      /* P7 scattering production cross section */
	      
	      ptr = (long)RDB[gcu + GCU_B1_SP7];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);  
	      val = BufVal(ptr, i, j);
	      AddStat(val/flx, ptr, i, j);
	    }
	}
    
      /***********************************************************************/

      /* Next universe */

      gcu = NextItem(gcu);

      /***********************************************************************/
    }

  /***************************************************************************/

  /***** Discontinuity factors ***********************************************/
  
  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check pointer */

      if ((adf = (long)RDB[gcu + GCU_PTR_ADF]) > VALID_PTR)
	{
	  /* Get areas of surface, mid-point and corner segments */
	  /* (lengths in 2D calculation) */

	  ptr = (long)RDB[adf + ADF_PTR_SURF_AREA];
	  ds = &RDB[ptr];
	  
	  ptr = (long)RDB[adf + ADF_PTR_MID_AREA];
	  dm = &RDB[ptr];
	  
	  ptr = (long)RDB[adf + ADF_PTR_CORN_AREA];
	  dc = &RDB[ptr];

	  /* Loop over energy groups */

	  for (i = 0; i < nfg; i++)
	    {
	      /* Homogeneous volume flux */

	      ptr = RDB[gcu + GCU_RES_FG_DF_HOM_VOL_FLUX];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      flx = BufVal(ptr, i);
	      AddStat(norm*flx, ptr, i);

	      /* Heterogeneous volume flux */

	      ptr = RDB[gcu + GCU_RES_FG_DF_HET_VOL_FLUX];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      flx = BufVal(ptr, i);
	      AddStat(norm*flx, ptr, i);

	      /* Loop over vertices */

	      for (m = 0; m < (long)RDB[adf + ADF_NSURF]; m++)
		{
		  /* Currents */

		  ptr = RDB[gcu + GCU_RES_FG_DF_SURF_IN_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);

		  ptr = RDB[gcu + GCU_RES_FG_DF_SURF_OUT_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);

		  ptr = RDB[gcu + GCU_RES_FG_DF_SURF_NET_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);

		  ptr = RDB[gcu + GCU_RES_FG_DF_MID_IN_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);
		  
		  ptr = RDB[gcu + GCU_RES_FG_DF_MID_OUT_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);

		  ptr = RDB[gcu + GCU_RES_FG_DF_MID_NET_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);

		  /* Homogeneous surface flux */

		  ptr = RDB[gcu + GCU_RES_FG_DF_HOM_SURF_FLUX];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  flx = BufVal(ptr, m, i);
		  AddStat(norm*flx, ptr, m, i);

		  /* Heterogeneous surface flux */

		  ptr = RDB[gcu + GCU_RES_FG_DF_HET_SURF_FLUX];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i)/ds[m];
		  AddStat(norm*val, ptr, m, i);

		  /* Discontinuity factor */

		  if (flx > 0.0)
		    {
		      ptr = RDB[gcu + GCU_RES_FG_DF_SURF_DF];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      AddStat(val/flx, ptr, m, i);
		    }
		}
	      
	      /* Loop over corners */

	      for (m = 0; m < (long)RDB[adf + ADF_NCORN]; m++)
		{
		  /* Currents */

		  ptr = RDB[gcu + GCU_RES_FG_DF_CORN_IN_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);
		  
		  ptr = RDB[gcu + GCU_RES_FG_DF_CORN_OUT_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);

		  ptr = RDB[gcu + GCU_RES_FG_DF_CORN_NET_CURR];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i);
		  AddStat(norm*val, ptr, m, i);

		  /* Homogeneous corner flux */

		  ptr = RDB[gcu + GCU_RES_FG_DF_HOM_CORN_FLUX];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  flx = BufVal(ptr, m, i);
		  AddStat(norm*flx, ptr, m, i);

		  /* Heterogeneous corner flux */

		  ptr = RDB[gcu + GCU_RES_FG_DF_HET_CORN_FLUX];
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  val = BufVal(ptr, m, i)/dc[m];
		  AddStat(norm*val, ptr, m, i);

		  /* Discontinuity factor */

		  if (flx > 0.0)
		    {
		      ptr = RDB[gcu + GCU_RES_FG_DF_CORN_DF];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      AddStat(val/flx, ptr, m, i);
		    }
		}
	    }
	}
      
      /* Next universe */
      
      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /***** Pin-power distribution form factors *********************************/
  
  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check pointer */

      if ((ppw = (long)RDB[gcu + GCU_PTR_PPW]) > VALID_PTR)
	{
	  /* Reset total power */

	  div = 0.0;
	  
	  /* Get number of pins */

	  np = (long)RDB[ppw + PPW_NP];

	  /* Loop over energy groups and pins */
	  
	  for (i = 0; i < nfg; i++)
	    for (m = 0; m < np; m++)
	      {
		/* Homogeneous flux */

		ptr = RDB[gcu + GCU_RES_FG_PPW_HOM_FLUX];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		flx = BufVal(ptr, m, i);
		AddStat(norm*flx, ptr, m, i);

		/* Power */
		    
		ptr = (long)RDB[gcu + GCU_RES_FG_PPW_POW];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		val = BufVal(ptr, m, i);
		
		/* Check value */
		
		if (flx != 0.0)
		  {
		    /* Store form factor */
		    
		    ptr = RDB[gcu + GCU_RES_FG_PPW_FF];
		    CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		    AddStat(val/flx, ptr, m, i);
		  }

		/* Add to total */
		    
		div = div + val;
	      }
	  
	  /* Fractional power */
	  
	  ptr = (long)RDB[gcu + GCU_RES_FG_PPW_POW];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	  for (i = 0; i < nfg; i++)
	    for (m = 0; m < np; m++)
	      {
		val = BufVal(ptr, m, i);
		AddStat(val/div, ptr, m, i);
	      }
	}

      /* Next universe */
      
      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /***** Albedos and partial albedos *****************************************/
  
  /* Loop over universes */

  gcu = (long)RDB[DATA_PTR_GCU0];
  while (gcu > VALID_PTR)
    {
      /* Check pointer */

      if ((alb = (long)RDB[gcu + GCU_PTR_ALB]) > VALID_PTR)
	{
	  /* Number of surfaces */

	  ns = (long)RDB[alb + ALB_NSURF];
	  CheckValue(FUNCTION_NAME, "ns", "", ns, 1, 8);

	  /* Total albedos */

	  ptr = RDB[gcu + GCU_RES_FG_TOT_ALB];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Loop over energy groups */

	  for (i1 = 0; i1 < nfg; i1++)
	    for (i2 = 0; i2 < nfg; i2++)
	      {
		val = BufVal(ptr, i1, i2);
		AddStat(val, ptr, i1, i2);
	      }

	  /* Currents and partial albedos */

	  for (m1 = 0; m1 < ns; m1++)
	    for (i1 = 0; i1 < nfg; i1++)
	      {
		/* Inward current */

		ptr = RDB[gcu + GCU_RES_FG_ALB_IN_CURR];
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		val = BufVal(ptr, m1, i1);
		AddStat(norm*val, ptr, m1, i1);

		/* Loop over second surface and group */

		for (m2 = 0; m2 < ns; m2++)
		  for (i2 = 0; i2 < nfg; i2++)
		    {
		      /* Outward current */
		      
		      ptr = RDB[gcu + GCU_RES_FG_ALB_OUT_CURR];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      val = BufVal(ptr, m1, i1, m2, i2);
		      AddStat(norm*val, ptr, m1, i1, m2, i2);

		      /* Partial albedo */

		      ptr = RDB[gcu + GCU_RES_FG_PART_ALB];
		      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		      val = BufVal(ptr, m1, i1, m2, i2);
		      AddStat(val, ptr, m1, i1, m2, i2);
		    }
	      }
	}

      /* Next universe */
      
      gcu = NextItem(gcu);
    }

  /***************************************************************************/

  /* For EDo (18.1.2014) */

  DiffCoefED(4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  /****************************************************************************/
}

/*****************************************************************************/

#ifdef __cplusplus 
} 
#endif 
