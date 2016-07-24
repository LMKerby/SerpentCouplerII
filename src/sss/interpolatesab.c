#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : interpolatesab.c                               */
/*                                                                           */
/* Created:       2015/02/12 (TVi)                                           */
/* Last modified: 2016/02/19 (TVi)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Interpolates S(a,b) cross sections and secondary particle    */
/*              distributions in temperature. Returns interpolated nuclide.  */
/*                                                                           */
/* Comments:    - Methodology is the same as in makxsf                       */
/*              - This mode is activated with "therm",                       */
/*                "thermstoch" uses stochastic mixing approach for           */
/*                 interpolation                                             */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InterpolateSab:"

long InterpolateSab(long ptr)
{
  long nuc, nuc1, nuc2, rea, rea1, rea2, ne1, ne2, erg1, erg2;
  long ace, ace1, ace2, i, j, k, xsd1, xsd2, loc, loc0, idx;
  long NXS1[16], JXS1[32], NXS2[16], JXS2[32], ptr1, ptr2, mt, eidx;
  double T, T1, T2, *E1, *E2, *XSS, *XSS1, *XSS2, f, fe, *xs, *xs1, *xs2, *dat, *dat1, *dat2;
  char newname[50];

  /* Avoid compiler warning */

  dat = NULL;
  dat1 = NULL;
  dat2 = NULL;
  
  /* Pointers to nuclides */  

  loc = (long)RDB[ptr + THERM_PTR_SAB];
  if(loc < VALID_PTR)
    Die(FUNCTION_NAME, "S(a,b) data not found for interpolation");
   
  nuc1 = (long)RDB[loc + SAB_PTR_ISO];

  loc = NextItem(loc);
  if(loc < VALID_PTR)
    Die(FUNCTION_NAME, "S(a,b) data not found for interpolation");
    
  nuc2 = (long)RDB[loc + SAB_PTR_ISO];

  CheckPointer(FUNCTION_NAME, "(nuc1)", DATA_ARRAY, nuc1);
  CheckPointer(FUNCTION_NAME, "(nuc2)", DATA_ARRAY, nuc2);

  /* Temperatures */

  T = RDB[ptr + THERM_T];
  T1 = RDB[nuc1 + NUCLIDE_TEMP];
  T2 = RDB[nuc2 + NUCLIDE_TEMP];

  CheckValue(FUNCTION_NAME, "SAB data T1","", T1, 0, 5000.00);
  CheckValue(FUNCTION_NAME, "SAB data T2","", T2, 0, 5000.00);

  /* Check if interpolation is unnecessary */
  if(T == T1)
    return nuc1;
  else if(T == T2)
    return nuc2;

  /* Check that interpolation is possible */

  if(T1 > T) 
    Error(ptr, "Temperature %1.1fK is below minimum %1.1fK", T, T1);
  else if (T2 < T) 
    Error(ptr, "Temperature %1.1fK is above maximum %1.1fK", T, T2);

  /* Message */
  
  fprintf(out,"Interpolating to %.2f K between S(a,b) nuclides %s @ %.2f K and %s @ %.2f K\n", 
	 T, GetText(nuc1 + NUCLIDE_PTR_NAME), T1, 
	 GetText(nuc2 + NUCLIDE_PTR_NAME), T2);

  /************* Read & init relevant data arrays *************/

  /* Pointer to ACE data of nuc1 */
  
  ace1 = (long)RDB[nuc1 + NUCLIDE_PTR_ACE];
  ace2 = (long)RDB[nuc2 + NUCLIDE_PTR_ACE];
   
  CheckPointer(FUNCTION_NAME, "(ace)", ACE_ARRAY, ace1);
  CheckPointer(FUNCTION_NAME, "(ace)", ACE_ARRAY, ace2);

  /* Read in NXS and JXS arrays */
 
  ptr1 = (long)ACE[ace1 + ACE_PTR_NXS];
  ptr2 = (long)ACE[ace2 + ACE_PTR_NXS];

  CheckPointer(FUNCTION_NAME, "(nxs)", ACE_ARRAY, ptr1);
  CheckPointer(FUNCTION_NAME, "(nxs)", ACE_ARRAY, ptr2);
 
  for(i=0; i<16; i++){
    NXS1[i] = (long)ACE[ptr1++];    
    NXS2[i] = (long)ACE[ptr2++];

    /* Sanity check: if the S(a,b) data are for the same nuclide,
       the NXS data should be the same excluding the number of 
       energy grid points which may vary (at least H in HZr) */ 

    if(i > 0 && NXS1[i] != NXS2[i])
      Die(FUNCTION_NAME, 
	   "Large discrepancies in S(a,b) data. (%s %s) (NXS1[%ld] = %ld != %ld = NXS2[%ld])", 
	   GetText(nuc1 + NUCLIDE_PTR_NAME), GetText(nuc2 + NUCLIDE_PTR_NAME), 
	   i, NXS1[i], NXS2[i], i);

  }

  ptr1 = (long)ACE[ace1 + ACE_PTR_JXS];
  ptr2 = (long)ACE[ace2 + ACE_PTR_JXS];

  CheckPointer(FUNCTION_NAME, "(jxs)", ACE_ARRAY, ptr1);
  CheckPointer(FUNCTION_NAME, "(jxs)", ACE_ARRAY, ptr2); 

  for(i=0; i<32; i++){
    JXS1[i] = (long)ACE[ptr1++];
    JXS2[i] = (long)ACE[ptr2++];    

    /* Sanity check: if the S(a,b) data are for the same nuclide,
       the JXS data should be the same excluding the number of 
       energy grid points which may vary (at least for H in HZr) */ 
        
    if(i != 4 && i != 5 && JXS1[i] != JXS2[i])
      Die(FUNCTION_NAME, 
	  "Large discrepancies in S(a,b) data. (%s %s) (JXS1[%ld] = %ld != %ld = JXS2[%ld])", 
	   GetText(nuc1 + NUCLIDE_PTR_NAME), GetText(nuc2 + NUCLIDE_PTR_NAME), 
	  i, JXS1[i], JXS2[i], i);
  }
  
  /************ Make a duplicate of nuc1 (interpolated nuclide) *************/

  /* Pitää tehdä kokonaan uusi nuklidi ja kopioida dataa, koska voi olla 
     että usammassa interpoloinnissa käytetään samaa S(a,b) -dataa, joten
     Sab-nuklideita ei voi ylikirjoittaa */

  /* Name of new nuclide is obtained by adding 'i' (=interpolated) 
     before the last character (='t') */

  /* Onkohan tää tarpeellista? Tunnistetaan kuitenkin myös T:n perusteella */
  
  sprintf(newname, "%st", GetText(nuc1 + NUCLIDE_PTR_NAME));
  newname[strlen(newname)-2]='i';
  
  /* Create new nuclide & copy data from nuc1 */

  nuc = NewItem(DATA_PTR_NUC0, NUCLIDE_BLOCK_SIZE);
  memcpy(&WDB[nuc + LIST_DATA_SIZE], &RDB[nuc1 + LIST_DATA_SIZE],
	 (NUCLIDE_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));

  /* Copy name */

  WDB[nuc + NUCLIDE_PTR_NAME] = PutText(newname);
  
  /* Allocate memory for previous collision velocity and relative energy */
  
  WDB[nuc + NUCLIDE_PREV_COL_Z2] = NULLPTR;
  WDB[nuc + NUCLIDE_PREV_COL_COS] = NULLPTR;
  WDB[nuc + NUCLIDE_PREV_COL_ER] = NULLPTR;
  
  AllocValuePair(nuc + NUCLIDE_PREV_COL_Z2);
  AllocValuePair(nuc + NUCLIDE_PREV_COL_COS);
  AllocValuePair(nuc + NUCLIDE_PREV_COL_ER);
  
  /* Nuclide temperature (to be interpolated) */
  
  WDB[nuc + NUCLIDE_TEMP] = T;

  /* Lib ID */
  /* Pidetään samana vaikka toi nimi muuttuukin */

  WDB[nuc + NUCLIDE_PTR_LIB_ID] = (double)PutText(GetText(nuc1 + NUCLIDE_PTR_LIB_ID));

  /* Make a copy of ACE array and update pointer */

  ace = ReallocMem(ACE_ARRAY, ACE_BLOCK_SIZE);
  memcpy(&ACE[ace], &ACE[ace1], ACE_BLOCK_SIZE*sizeof(double));
  
  WDB[nuc + NUCLIDE_PTR_ACE] = (double)ace;

  /* Put new name and alias  */
	  
  ACE[ace + ACE_PTR_NAME] = PutText(newname);
  ACE[ace + ACE_PTR_ALIAS] = PutText(newname);

  /*** Copy XSS array ******************************************************/

  /* Pointer to XSS array of old nuclide */

  loc = (long)ACE[ace1 + ACE_PTR_XSS];
  CheckPointer(FUNCTION_NAME, "(xss)", ACE_ARRAY, loc);
  
  /* Alloc mem for new XSS array & copy */

  loc0 = ReallocMem(ACE_ARRAY, NXS1[0]);
  memcpy(&ACE[loc0], &ACE[loc], NXS1[0]*sizeof(double));
    
  /* Pitäisikö myös NXS ja JXS -arraystä tehdä kopio??? */

  ACE[ace + ACE_PTR_XSS] = (double)loc0;

  /****** Init XSS pointers ************************************************/
  /* Pointers to XSS data (new, old1, old2). 
     Toi ace1 -pointteri tarkistettiin just */

  CheckPointer(FUNCTION_NAME, "(xss0)", ACE_ARRAY, loc0);  
  CheckPointer(FUNCTION_NAME, "(xss1)", ACE_ARRAY, (long)ACE[ace1 + ACE_PTR_XSS]);  
  CheckPointer(FUNCTION_NAME, "(xss2)", ACE_ARRAY, (long)ACE[ace2 + ACE_PTR_XSS]);
   
  XSS = &ACE[loc0];
  XSS1 = &ACE[(long)ACE[ace1 + ACE_PTR_XSS]];
  XSS2 = &ACE[(long)ACE[ace2 + ACE_PTR_XSS]];
  
  /************* Interpolate all reactions *********************************/

  /* Reset REA list of new nuclide */

  WDB[nuc + NUCLIDE_PTR_REA] = NULLPTR;

  /* Temperature interpolation factor */
  
  f = (T-T1)/(T2-T1);

  rea1 = (long)RDB[nuc1 + NUCLIDE_PTR_REA];

  while(rea1 > VALID_PTR){

    /* Find  corresponding reaction for nuc2 */

    rea2=(long)RDB[nuc2 + NUCLIDE_PTR_REA];    
    mt = (long)RDB[rea1 + REACTION_MT];

    while(rea2 > VALID_PTR){
      
      if( (long)RDB[rea2 + REACTION_MT] == mt)
	break;
      
      rea2=NextItem(rea2);
    }

    if(rea2 < VALID_PTR)
      Die(FUNCTION_NAME, "mt = %ld not found for S(a,b) nuclide %s", 
	  mt, GetText(nuc2 + NUCLIDE_PTR_NAME));
         
    /* Get pointers to energy grids and numbers of grid points */
    /* (Pointers in XSS array) */

    ne1 = (long)RDB[rea1 + REACTION_XS_NE];
    erg1 = (long)RDB[rea1 + REACTION_PTR_EGRID];
    xsd1 = (long)RDB[rea1 + REACTION_PTR_XS];

    ne2 = (long)RDB[rea2 + REACTION_XS_NE];
    erg2 = (long)RDB[rea2 + REACTION_PTR_EGRID];
    xsd2 = (long)RDB[rea2 + REACTION_PTR_XS];
   
    /* Init new reaction based on rea1 */
   
    rea = NewItem(nuc + NUCLIDE_PTR_REA, REACTION_BLOCK_SIZE);
    memcpy(&WDB[rea + LIST_DATA_SIZE], &RDB[rea1 + LIST_DATA_SIZE],
	   (REACTION_BLOCK_SIZE - LIST_DATA_SIZE)*sizeof(double));
    
    WDB[rea + REACTION_PTR_NUCLIDE] = nuc;
    
    /* Incident energy arrays */
       
    E1 = &XSS1[erg1];
    E2 = &XSS2[erg2];

    /* Cross section arrays */

    xs  = &XSS[xsd1];
    xs1 = &XSS1[xsd1];
    xs2 = &XSS2[xsd2];
    
    /* Data arrays for secondary energies and cosines */

    if(mt == 1004){
      dat = &XSS[JXS1[2]-1];
      dat1 = &XSS1[JXS1[2]-1];
      dat2 = &XSS2[JXS2[2]-1];
    }
    else if(mt == 1002){
      dat  = &XSS[JXS1[5]-1];
      dat1 = &XSS1[JXS1[5]-1];
      dat2 = &XSS2[JXS2[5]-1];
    }
    else
      Die(FUNCTION_NAME, "unknown S(a,b) reaction");

    /* Reset index over energy grid of nuclide nuc2 (in case of 
       different energy grids at different temperatures) */ 

    eidx = 0;

    /* Reset index over secondary particle data */
    idx = 0;

    /********* Interpolate XS and outgoing energy & cosine tables ************/

    /* Pienellä käänteisinsinööreilyllä oikeat makxsf:n käyttämät 
       interpolointitavat löytyivät:

       - vaikutusalat interpoloidaan lineaarisesti.
       - epäelastisen sironnan tapauksessa energiat interpoloidaan samalla
       tavalla kuin rinnan olevat vastukset.
       - epäelastisen ja elastisen sironnan tapauksessa kosinitaulukot 
       interpoloidaan lineaarisesti. 
       
       Pientä (ilmeisesti numeerista) virhettä tulee suhteessa makxsf:n,
       mutta sekin esiintyy vain erittäin pienillä kosineilla */ 

    /* Loop over XS energy grid (Ein values) */
    
    for(i=0; i<ne1; i++){
      
      /* Usein energiagridit on samat eri lämpötiloissa, mutta esim H in HZr:llä 
	 se näyttää vaihtelevan. Eri gridien tapauksessa pitää interpoloida myös
	 energian suhteen. Tallennetaan lopputulos nuklidin 1 energiagridiä käyttäen. */

      /* Check energies */

      if( (i > ne2-1) || (E1[i] != E2[i]) ){
	
	/* Find correct index from nuclide 2 grid */
	
	while( (eidx < ne2-2) && (E2[eidx+1] < E1[i]) )
	  eidx++;
	
	/* Check if last energy point was passed */
	
	if( (eidx == ne2-2) && E2[eidx+1] < E1[i] ){
	  /* Exceeded, use last value */
	  eidx = ne2-2;
	  fe = 1.0;
	}
	else
	  {
	    /* Calculate energy interpolation factor */
	    
	    fe = (E1[i]-E2[eidx])/(E2[eidx+1]-E2[eidx]);
	  }
	
	/********* Linear interpolation in both temperature and energy *****/
	
	xs[i] = (1-f)*xs1[i] + f*( (1-fe)*xs2[eidx] + fe*xs2[eidx+1]); 	
	
	/***** Interpolation of secondary particle tables in temperature *****/
	
	/*** Inelastic mode */
	
	if(mt == 1004){
	  
	  /* Loop over outgoing energies */
	  
	  for(k=0; k<NXS1[3]; k++){
	    
	    /*** In case of inelastic scattering mode, interpolate both energies 
		 and cosines */
	    
	    /* Interpolate outgoing energy using reciprocal law */
	    /* outgoing energy at E_in E1[i] and temperature T2 is interpolated 
	       linearly in energy */
	    
	    dat[idx] = 1/( (1-f)/dat1[idx] + f/( (1-fe)*dat2[eidx]+fe*dat2[eidx+1] ) );
	    idx++;
	    
	    /* Interpolate cosines linearly */

	    /* outgoing cosine at E_in E1[i] and temperature T2 is interpolated 
	       linearly in energy */
	    
	    for(j=0; j<NXS1[2]+1; j++){
	      dat[idx] = (1-f)*dat1[idx] + f*( (1-fe)*dat2[eidx] + fe*dat2[eidx+1]);
	      idx++;
	    }
	    
	  }
	}
	
	/* Elastic mode */

	else if(mt == 1002){

	  /* Normal treatment. In case of exact treatment (NXS1[4] == 4)
	     only xs needs to be interpolated*/

	  if(NXS1[4] != 4 ){

	    for(j=0; j<NXS1[5]+1; j++){

	      /* In case of elastic scattering mode, interpolate only cosines */
	      
	      dat[idx] = (1-f)*dat1[idx] + f*( (1-fe)*dat2[eidx] + fe*dat2[eidx+1]);
	      idx++;

	    }
	  }
	  
	}
      }         
      else {

	/* Tää on periaatteessa turha koska tuolla eriävien energioiden 
	   tapauksella voidaan hoitaa myös sama gridi -keissi, mutta jätin 
	   optimointisyistä paikoilleen (yhtenevät energiagridit on 
	   yleisin tapaus) */

	/** Simple, linear interpolation of cross sections in temperature **/

      	xs[i] = (1-f)*xs1[i] + f*xs2[i];
	
	/***** Interpolation of secondary particle tables in temperature *****/

	/*** Inelastic mode ***/
	
	if(mt == 1004){
	  
	  /* Loop over outgoing energies */
	  for(k=0; k<NXS1[3]; k++){
	    
	    /* In case of inelastic scattering mode, interpolate both energies 
	       and cosines */
	    
	    /* Interpolate outgoing energy using reciprocal law */

	    dat[idx] =1/( (1-f)/dat1[idx] + f/dat2[idx] ); 
	    idx++;
	    
	    /* Interpolate cosines linearly */
	    
	    for(j=0; j<NXS1[2]+1; j++){
	      dat[idx] = (1-f)*dat1[idx] + f*dat2[idx];	      	      
	      idx++;	      
	    }
	  }
	}
	/*** Elastic mode ****/
	
	/* In case of elastic scattering mode, interpolate only cosines */

	else if(mt == 1002){
	  
	  /* Normal treatment. (In case of exact treatment, do nothing 
	    since nothing to interpolate in data) */

	  if(NXS1[4] != 4 ){

	    for(j=0; j<NXS1[5]+1; j++){
	     
	      /* Simple linear interpolation */

	      dat[idx] = (1-f)*dat1[idx] + f*dat2[idx];            
	      idx++;

	    }
	  }  
	}
      }
    	        
    }
   
    rea1 = NextItem(rea1);
  }

  /* For checking purposes obtain pointer to nuclide using addnuclide() */
  
  if( (nuc = AddNuclide(newname, RDB[nuc1 + NUCLIDE_ZAI], 
			GetText(nuc1 + NUCLIDE_PTR_LIB_ID), T, NUCLIDE_TYPE_SAB, NO)) == 0){
    Die(FUNCTION_NAME, "Interpolated S(a,b) nuclide %s not found", 
	newname);
  }

  return nuc;
}
#ifdef __cplusplus 
} 
#endif 
