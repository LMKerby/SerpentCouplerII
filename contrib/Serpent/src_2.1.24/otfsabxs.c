/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : otfsabxs.c                                     */
/*                                                                           */
/* Created:       2015/03/20 (TVi)                                           */
/* Last modified: 2015/05/19 (TVi)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Returns interpolated cross sections for on-the-fly S(a,b)    */
/*              treatment.                                                   */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "OTFSabXS:"

double OTFSabXS(long rea, double E, double T, long id){

  long nuc, sab1, sab2, iso1, iso2, mt, rea2, sab0, ptr, ncol; 
  double f, xs1, xs2, T1, T2, xs;
 
  /* Avoid compiler warning */

  sab1 = -1;

  /* Get reaction mt and rea nuclide */

  mt = (long)RDB[rea + REACTION_MT];
  nuc = (long)RDB[rea + REACTION_PTR_NUCLIDE];

  if(E > RDB[nuc + NUCLIDE_SAB_EMAX]){
    return MicroXS(rea, E, id);
  }

  /* Check that data exists */

  if( (sab2 = (long)RDB[nuc + NUCLIDE_PTR_SAB]) < VALID_PTR ){
    Die(FUNCTION_NAME, "S(a,b) data not available for nuclide %s", GetText(nuc + NUCLIDE_PTR_NAME));
  }
  
  /* Get collision number */
  
  ptr = (long)RDB[DATA_PTR_COLLISION_COUNT];
  ncol = (long)GetPrivateData(ptr, id);

  /* Find correct temperature */
  
  while(sab2 > VALID_PTR){

    if(RDB[sab2 + SAB_T] > T){
      sab1 = PrevItem(sab2);
      break;
    }

    sab2 = NextItem(sab2);
  }
  
  /* Check that sab was found */

  if(sab2 < VALID_PTR)
    Die(FUNCTION_NAME, "S(a,b) OTF nuclide not found for %s", GetText(nuc + NUCLIDE_PTR_NAME));

  /* Temperatures */

  T1 = RDB[sab1 + SAB_T];
  T2 = RDB[sab2 + SAB_T];

  /* Pointers to S(a,b) nuclides */

  iso1=RDB[sab1 + SAB_PTR_ISO];
  iso2=RDB[sab2 + SAB_PTR_ISO];

  CheckPointer(FUNCTION_NAME, "(iso1)", DATA_ARRAY, iso1);
  CheckPointer(FUNCTION_NAME, "(iso2)", DATA_ARRAY, iso2);

  /* Find reaction pointers in S(a,b) data */
  
  /* Elastic scattering xs = total xs in case of S(a,b) nuclides 
     (which only have 1004 and 1002 reactions and total is 
     calculated over them) */
  
  if(mt == 1 || mt == 2 )
    rea2 = (long)RDB[iso1 + NUCLIDE_PTR_TOTXS];
  
  /* Find reaction at first temperature */

  else{

    rea2 = (long)RDB[iso1 + NUCLIDE_PTR_REA];
    
    while(rea2 > VALID_PTR){
      
      if((long)RDB[rea2 + REACTION_MT] == mt-1000)
	break;
      
      rea2 = NextItem(rea2);
    }
    
    if(rea2 < VALID_PTR)
      Die(FUNCTION_NAME, "mt %ld not found for S(a,b) nuclide %s", mt, 
	  GetText(iso1 + NUCLIDE_PTR_NAME));
        
  }

  xs1 = MicroXS(rea2, E, id);

  if(mt == 1 || mt == 2 )
    rea2 = (long)RDB[iso2 + NUCLIDE_PTR_TOTXS];

  else{
    /* Find reaction at second temperature */
    
    rea2 = (long)RDB[iso2 + NUCLIDE_PTR_REA];
    
    while(rea2 > VALID_PTR){
      
      if((long)RDB[rea2 + REACTION_MT] == mt-1000)
	break;
      
      rea2 = NextItem(rea2);
    }
    
    if(rea2 < VALID_PTR)
      Die(FUNCTION_NAME, "mt %ld not found for S(a,b) nuclide %s", mt, 
	  GetText(iso1 + NUCLIDE_PTR_NAME));
  
  }

  xs2 = MicroXS(rea2, E, id);

  f = (T-T1)/(T2-T1);

  /* Check values */

  CheckValue(FUNCTION_NAME, "f", "", f, 0.0, 1.0);
  CheckValue(FUNCTION_NAME, "xs1", "", f, 0.0, 1.0E10);
  CheckValue(FUNCTION_NAME, "xs2", "", f, 0.0, 1.0E10);

  /* Pointer to first item in sab list */

  sab0 = (long)RDB[nuc + NUCLIDE_PTR_SAB];

  /* Store values */

  StoreValuePair(sab0 + SAB_PTR_PREV_FRAC, ncol, f, id);
  StoreValuePair(sab0 + SAB_PTR_PREV_SAB1, ncol, (double)sab1, id);

  xs = xs1 + f*(xs2-xs1);
  
  StoreValuePair(rea + REACTION_PTR_PREV_XS, E, xs, id);
  StoreValuePair(rea + REACTION_PTR_PREV_XS0, E, xs, id);  

  /* Näillä voisi varmaan optimoida myös tuota rutiinin alkupäätä */

  return xs;
}
