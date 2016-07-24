#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : addbranching.c                                 */
/*                                                                           */
/* Created:       2010/09/11 (JLe)                                           */
/* Last modified: 2015/06/09 (JLe)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Creates branches to isomeric states, secondary products      */
/*              and fission product yields at different energies             */
/*                                                                           */
/* Comments: - Noiden fixed-arvojen lisäksi data pitää lukea käyttäjältä     */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "AddBranching:"

/*****************************************************************************/

void AddBranching(long nuc)
{
  long mt, mt1, rea, loc0, loc1, ptr, yld, n;
  double br, Q;
    
  /* Check burnup mode */

  if (((long)RDB[DATA_BURNUP_CALCULATION_MODE] == NO) &&
      ((long)RDB[DATA_XENON_EQUILIBRIUM_MODE] == -1) &&
      ((long)RDB[DATA_SAMARIUM_EQUILIBRIUM_MODE] == -1))
    return;

  /* Check DBRC flag */

  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_DBRC)
    return;

  /***************************************************************************/

  /***** Energy-dependent branching ratios ***********************************/

  /* Check type */

  if ((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT)
    {
      /* Loop over data */

      loc0 = (long)RDB[DATA_PTR_BRA_LIST];
      while (loc0 > VALID_PTR)
	{
	  /* Compare ZAI */

	  if ((long)RDB[loc0 + BRA_LIST_ZAI] != (long)RDB[nuc + NUCLIDE_ZAI])
	    {
	      /* Next in list */

	      loc0 = NextItem(loc0);

	      /* Cycle loop */

	      continue;
	    }

	  /* Find reaction */

	  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
	  while (rea > VALID_PTR)
	    {
	      /* Compare mt */

	      if ((long)RDB[loc0 + BRA_LIST_MT] == 
		  (long)RDB[rea + REACTION_MT])
		break;
	      
	      /* Next reaction */
	      
	      rea = NextItem(rea);
	    }

	  /* Check if found */

	  if (rea < VALID_PTR)
	    {
	      /* Next in list */

	      loc0 = NextItem(loc0);

	      /* Cycle loop */

	      continue;
	    }
	  
	  /* Set flag and type */

	  SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_BRA_DATA);
	  WDB[nuc + NUCLIDE_BRA_TYPE] = (double)BRA_TYPE_ENE;

	  /* Get pointer to first state */

	  loc1 = (long)RDB[loc0 + BRA_LIST_PTR_STATES];
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  /* Check and put pointer */

	  if ((long)RDB[loc1 + BRA_STATE_LFS] != 0)
	    Die(FUNCTION_NAME, "Not a ground state");
	  else
	    WDB[rea + REACTION_PTR_BRA_STATE] = (double)loc1;

	  /* Get pointer to second state */

	  loc1 = NextItem(loc1);
	  CheckPointer(FUNCTION_NAME, "(loc1)", DATA_ARRAY, loc1);

	  /* Loop over remaining */

	  while (loc1 > VALID_PTR)
	    {
	      /* Create duplicate */
	  
	      ptr = DuplicateItem(rea);
	  
	      /* Check state */

	      if ((long)RDB[loc1 + BRA_STATE_LFS] == 0)
		Die(FUNCTION_NAME, "Ground state");
	  
	      /* Put state and pointer to data */
	      /*
	      WDB[ptr + REACTION_RFS] = RDB[loc1 + BRA_STATE_LFS];
	      */
	      WDB[ptr + REACTION_RFS] = 1.0;
	      WDB[ptr + REACTION_PTR_BRA_STATE] = (double)loc1;	      
	  
	      /* Set pointer to parent reaction */
	      /*
	      WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
	      */
	      /* Set type */
	      /*
	      WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
	      */
	      /* Set branch mt */
	      /*
	      WDB[ptr + REACTION_BRANCH_MT] = RDB[ptr + REACTION_MT];
	      */
	      /* Next state */

	      loc1 = NextItem(loc1);
	    }

	  /* Next in list */

	  loc0 = NextItem(loc0);
	}
    }

  /***************************************************************************/

  /***** Fixed branching ratios to isomeric states ***************************/

  /* Check type and if flag is already set */

  if (((long)RDB[nuc + NUCLIDE_TYPE] == NUCLIDE_TYPE_TRANSPORT) &&
      (!((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_BRA_DATA)))
    {
      /* Check fixed values */
      
      switch((long)RDB[nuc + NUCLIDE_ZAI])
	{

#ifdef OLD_FIXED_BRANCHING_RATIOS

	  /* These values are taken from Serpent 1 and originate from some  */
	  /* long-forgotten source. */
	  
	case 952410:
	  {
	    br = 0.885;
	    mt = 102;
	    break;
	  }
	case 952430:
	  {
	    br = 0.050;
	    mt = 102;
	    break;
	  }
	case 611470:
	  {
	    br = 0.530;
	    mt = 102;
	    break;
	  }
	case 471090:
	  {
	    br = 0.941;
	    mt = 102;
	    break;
	  }

#else
	  /* These values are calculated from energy-dependent branching  */
	  /* ratios in the JEFF-3.1 activation file in PWR flux spectrum. */

	case 110230:
	  {
	    br = 0.2320;
	    mt = 102;
	    break;
	  }
	case 170370:
	  {
	    br = 0.8809;
	    mt = 102;
	    break;
	  }
	case 210450:
	  {
	    br = 0.5560;
	    mt = 102;
	    break;
	  }
	case 270590:
	  {
	    br = 0.4440;
	    mt = 102;
	    break;
	  }
	case 320720:
	  {
	    br = 0.5012;
	    mt = 102;
	    break;
	  }
	case 320740:
	  {
	    br = 0.6660;
	    mt = 102;
	    break;
	  }
	case 320760:
	  {
	    br = 0.4005;
	    mt = 102;
	    break;
	  }
	case 340760:
	  {
	    br = 0.7409;
	    mt = 102;
	    break;
	  }
	case 340780:
	  {
	    br = 0.1178;
	    mt = 102;
	    break;
	  }
	case 340800:
	  {
	    br = 0.8454;
	    mt = 102;
	    break;
	  }
	case 340820:
	  {
	    br = 0.1402;
	    mt = 102;
	    break;
	  }
	case 350790:
	  {
	    br = 0.7687;
	    mt = 102;
	    break;
	  }
	case 350810:
	  {
	    br = 0.0914;
	    mt = 102;
	    break;
	  }
	case 360780:
	  {
	    br = 0.9704;
	    mt = 102;
	    break;
	  }
	case 360800:
	  {
	    br = 0.6031;
	    mt = 102;
	    break;
	  }
	case 360820:
	  {
	    br = 0.3330;
	    mt = 102;
	    break;
	  }
	case 360840:
	  {
	    br = 0.1839;
	    mt = 102;
	    break;
	  }
	case 370850:
	  {
	    br = 0.8791;
	    mt = 102;
	    break;
	  }
	case 380840:
	  {
	    br = 0.2530;
	    mt = 102;
	    break;
	  }
	case 380860:
	  {
	    br = 0.1988;
	    mt = 102;
	    break;
	  }
	case 390890:
	  {
	    br = 0.9979;
	    mt = 102;
	    break;
	  }
	case 390900:
	  {
	    br = 0.7496;
	    mt = 102;
	    break;
	  }
	case 410930:
	  {
	    br = 0.3101;
	    mt = 102;
	    break;
	  }
	case 410940:
	  {
	    br = 0.9610;
	    mt = 102;
	    break;
	  }
	case 420920:
	  {
	    br = 0.9978;
	    mt = 102;
	    break;
	  }
	case 451030:
	  {
	    br = 0.9240;
	    mt = 102;
	    break;
	  }
	case 451050:
	  {
	    br = 0.9040;
	    mt = 102;
	    break;
	  }
	case 461060:
	  {
	    br = 0.9527;
	    mt = 102;
	    break;
	  }
	case 461080:
	  {
	    br = 0.9779;
	    mt = 102;
	    break;
	  }
	case 461100:
	  {
	    br = 0.8500;
	    mt = 102;
	    break;
	  }
	case 471070:
	  {
	    br = 0.9898;
	    mt = 102;
	    break;
	  }
	case 471090:
	  {
	    br = 0.9540;
	    mt = 102;
	    break;
	  }
	case 481100:
	  {
	    br = 0.9945;
	    mt = 102;
	    break;
	  }
	case 481120:
	  {
	    br = 0.8685;
	    mt = 102;
	    break;
	  }
	case 481140:
	  {
	    br = 0.8812;
	    mt = 102;
	    break;
	  }
	case 481160:
	  {
	    br = 0.6660;
	    mt = 102;
	    break;
	  }
	case 491130:
	  {
	    br = 0.4191;
	    mt = 102;
	    break;
	  }
	case 501120:
	  {
	    br = 0.7253;
	    mt = 102;
	    break;
	  }
	case 501160:
	  {
	    br = 0.9568;
	    mt = 102;
	    break;
	  }
	case 501180:
	  {
	    br = 0.9794;
	    mt = 102;
	    break;
	  }
	case 501200:
	  {
	    br = 0.9875;
	    mt = 102;
	    break;
	  }
	case 501220:
	  {
	    br = 0.0112;
	    mt = 102;
	    break;
	  }
	case 501240:
	  {
	    br = 0.0375;
	    mt = 102;
	    break;
	  }
	case 501260:
	  {
	    br = 0.3018;
	    mt = 102;
	    break;
	  }
	case 511210:
	  {
	    br = 0.9369;
	    mt = 102;
	    break;
	  }
	case 521200:
	  {
	    br = 0.8871;
	    mt = 102;
	    break;
	  }
	case 521220:
	  {
	    br = 0.6448;
	    mt = 102;
	    break;
	  }
	case 521240:
	  {
	    br = 0.9912;
	    mt = 102;
	    break;
	  }
	case 521260:
	  {
	    br = 0.8689;
	    mt = 102;
	    break;
	  }
	case 521280:
	  {
	    br = 0.9245;
	    mt = 102;
	    break;
	  }
	case 521300:
	  {
	    br = 0.8559;
	    mt = 102;
	    break;
	  }
	case 521320:
	  {
	    br = 0.8517;
	    mt = 102;
	    break;
	  }
	case 531290:
	  {
	    br = 0.4130;
	    mt = 102;
	    break;
	  }
	case 531310:
	  {
	    br = 0.9839;
	    mt = 102;
	    break;
	  }
	case 541240:
	  {
	    br = 0.8300;
	    mt = 102;
	    break;
	  }
	case 541260:
	  {
	    br = 0.8691;
	    mt = 102;
	    break;
	  }
	case 541280:
	  {
	    br = 0.8923;
	    mt = 102;
	    break;
	  }
	case 541300:
	  {
	    br = 0.9164;
	    mt = 102;
	    break;
	  }
	case 541320:
	  {
	    br = 0.8867;
	    mt = 102;
	    break;
	  }
	case 541330:
	  {
	    br = 0.9600;
	    mt = 102;
	    break;
	  }
	case 541340:
	  {
	    br = 0.9853;
	    mt = 102;
	    break;
	  }
	case 551330:
	  {
	    br = 0.9070;
	    mt = 102;
	    break;
	  }
	case 551340:
	  {
	    br = 0.9960;
	    mt = 102;
	    break;
	  }
	case 551350:
	  {
	    br = 0.9840;
	    mt = 102;
	    break;
	  }
	case 551370:
	  {
	    br = 0.9021;
	    mt = 102;
	    break;
	  }
	case 561300:
	  {
	    br = 0.8871;
	    mt = 102;
	    break;
	  }
	case 561320:
	  {
	    br = 0.9175;
	    mt = 102;
	    break;
	  }
	case 561340:
	  {
	    br = 0.9263;
	    mt = 102;
	    break;
	  }
	case 561350:
	  {
	    br = 0.9978;
	    mt = 102;
	    break;
	  }
	case 561360:
	  {
	    br = 0.9731;
	    mt = 102;
	    break;
	  }
	case 581360:
	  {
	    br = 0.8662;
	    mt = 102;
	    break;
	  }
	case 581380:
	  {
	    br = 0.9787;
	    mt = 102;
	    break;
	  }
	case 591410:
	  {
	    br = 0.6519;
	    mt = 102;
	    break;
	  }
	case 591430:
	  {
	    br = 0.3100;
	    mt = 102;
	    break;
	  }
	case 611470:
	  {
	    br = 0.5330;
	    mt = 102;
	    break;
	  }
	case 631530:
	  {
	    br = 0.9840;
	    mt = 102;
	    break;
	  }
	case 661640:
	  {
	    br = 0.3700;
	    mt = 102;
	    break;
	  }
	case 671650:
	  {
	    br = 0.9490;
	    mt = 102;
	    break;
	  }
	case 681660:
	  {
	    br = 0.2503;
	    mt = 102;
	    break;
	  }
	case 711750:
	  {
	    br = 0.3331;
	    mt = 102;
	    break;
	  }
	case 711760:
	  {
	    br = 0.9990;
	    mt = 102;
	    break;
	  }
	case 721790:
	  {
	    br = 0.9910;
	    mt = 102;
	    break;
	  }
	case 741820:
	  {
	    br = 0.8699;
	    mt = 102;
	    break;
	  }
	case 741840:
	  {
	    br = 0.9983;
	    mt = 102;
	    break;
	  }
	case 751850:
	  {
	    br = 0.9990;
	    mt = 102;
	    break;
	  }
	case 751870:
	  {
	    br = 0.9729;
	    mt = 102;
	    break;
	  }
	case 791970:
	  {
	    br = 0.9990;
	    mt = 102;
	    break;
	  }
	case 801960:
	  {
	    br = 0.9660;
	    mt = 102;
	    break;
	  }
	case 801980:
	  {
	    br = 0.9918;
	    mt = 102;
	    break;
	  }
	case 822060:
	  {
	    br = 0.9783;
	    mt = 102;
	    break;
	  }
	case 832090:
	  {
	    br = 0.6791;
	    mt = 102;
	    break;
	  }
	case 912330:
	  {
	    br = 0.4871;
	    mt = 102;
	    break;
	  }
	case 922340:
	  {
	    br = 0.5000;
	    mt = 102;
	    break;
	  }
	case 932350:
	  {
	    br = 0.4000;
	    mt = 102;
	    break;
	  }
	case 932390:
	  {
	    br = 0.3573;
	    mt = 102;
	    break;
	  }
	case 942360:
	  {
	    br = 0.5001;
	    mt = 102;
	    break;
	  }
	case 952410:
	  {
	    br = 0.9190;
	    mt = 102;
	    break;
	  }
	case 952430:
	  {
	    br = 0.0626;
	    mt = 102;
	    break;
	  }
	case 972470:
	  {
	    br = 0.4000;
	    mt = 102;
	    break;
	  }
	case 992530:
	  {
	    br = 0.0320;
	    mt = 102;
	    break;
	  }
	case 992550:
	  {
	    br = 0.9840;
	    mt = 102;
	    break;
	  }
	  
#endif
	  
	default:
	  {
	    br = -1.0;
	    mt = -1;
	    break;
	  }
	}
      
      /* Check ratio */
      
      if (br != -1.0)
	{
	  /* Loop over reaction channels to find correct mt */
	  
	  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
	  while (rea > VALID_PTR)
	    {
	      /* Check MT */
	      
	      if ((long)RDB[rea + REACTION_MT] == mt)
		break;
	      
	      /* Next */
	      
	      rea = NextItem(rea);
	    }
	  
	  /* Check if found */
	  
	  if (rea < 0)
	    Error(0, "No isomeric branching reaction %ld for nuclide %ld",
		  mt, (long)RDB[nuc + NUCLIDE_ZAI]);
	  
	  /* Set flag and type */
	  
	  SetOption(nuc + NUCLIDE_TYPE_FLAGS, NUCLIDE_FLAG_BRA_DATA);
	  WDB[nuc + NUCLIDE_BRA_TYPE] = (double)BRA_TYPE_FIX;
	  
	  /* Set branching ratio to ground state */
	  
	  WDB[rea + REACTION_BR] = (double)br;
	  
	  /* Create duplicate */
	  
	  ptr = DuplicateItem(rea);
	  
	  /* Set branching ratio to isomeric state */
	  
	  WDB[ptr + REACTION_BR] = 1.0 - br;
	  WDB[ptr + REACTION_RFS] = 1.0;
	  
	  /* Set pointer to parent reaction */
	  
	  WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
	  
	  /* Set type */
	  
	  WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
	  
	  /* Set branch mt */
	  
	  WDB[ptr + REACTION_BRANCH_MT] = RDB[ptr + REACTION_MT];
	}
    }
  
  /***************************************************************************/
  
  /***** Branching to secondary products *************************************/
  
  /* Loop over reactions */
  
  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check type */
      
      if (((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_TRA_BRANCH) &&
	  ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_SPECIAL))
	{
	  /* Loop over successive decay modes */
	  
	  for (n = 0; n < 5; n++)
	    {
	      /* Avoid compiler warning */
	      
	      mt = -1;
	      
	      /* Get mt */
	      
	      if (n == 0)
		mt = (long)RDB[rea + REACTION_MT];
	      else if (n == 1)
		mt = (long)RDB[rea + REACTION_RTYP2];
	      else if (n == 2)
		mt = (long)RDB[rea + REACTION_RTYP3];
	      else if (n == 3)
		mt = (long)RDB[rea + REACTION_RTYP4];
	      else if (n == 4)
		mt = (long)RDB[rea + REACTION_RTYP5];
	      else
		Die(FUNCTION_NAME, "Overflow");
	      
	      /* Break if no secondary modes */
	      
	      if ((mt == 10000) || (mt == 0))
		break;
	      
	      /* Put mt for partial modes */
	      
	      if ((mt > 599) && (mt < 650))
		mt1 = 103;
	      else if ((mt > 649) && (mt < 700))
		mt1 = 104;
	      else if ((mt > 699) && (mt < 750))
		mt1 = 105;
	      else if ((mt > 749) && (mt < 800))
		mt1 = 106;
	      else if ((mt > 799) && (mt < 850))
		mt1 = 107;
	      else
		mt1 = mt;
	      
	      /* Neutron reactions */
	      
	      switch (mt1)
		{
		case 22:
		  {
		    /*(n,na) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 23:
		  {
		    /* (n, n3a) */
		
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 3.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 24:
		  {
		    /* (n, 2na) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 25:
		  {
		    /* (n, 3na) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 28:
		  {
		    /* (n, np) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 29:
		  {
		    /* (n, n2a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 30:
		  {
		    /* (n, 2n2a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;

		    break;
		  }
		case 32:
		  {
		    /* (n, nd) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1002;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 33:
		  {
		    /* (n, nt) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1003;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 34:
		  {
		    /* (n, nHe-3) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2003;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 35:
		  {
		    /* (n, nd2a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1002;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 36:
		  {
		    /* (n, nt2a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1003;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;		
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 41:
		  {
		    /* (n,2np) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 42:
		  {
		    /* (n,3np) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 44:
		  {
		    /* (n,n2p) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 45:
		  {
		    /* (n, npa) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 103:
		  {
		    /* (n,p) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 104:
		  {
		    /* (n,d) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1002;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 105:
		  {
		    /* (n, t) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1003;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 106:
		  {
		    /* (n, He-3) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2003;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 107:
		  {
		    /* (n, a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 108:
		  {
		    /* (n, 2a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 109:
		  {
		    /* (n, 3a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 3.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 111:
		  {
		    /* (n,2p) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 112:
		  {
		    /* (n, pa) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;		
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 113:
		  {
		    /* (n, t2a) - Special treatment for B-10 */
		    
		    if ((long)RDB[nuc + NUCLIDE_ZAI] != 50100)
		      {
			ptr = DuplicateItem(rea);	    
			WDB[ptr + REACTION_MT] = 20000 + 1003;
			WDB[ptr + REACTION_RFS] = 0.0;
			WDB[ptr + REACTION_TYPE] = 
			  (double)REACTION_TYPE_TRA_BRANCH;
			WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
			WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
			WDB[ptr + REACTION_RTYP2] = 0.0;
			WDB[ptr + REACTION_RTYP3] = 0.0;
			WDB[ptr + REACTION_RTYP4] = 0.0;
			WDB[ptr + REACTION_RTYP5] = 0.0;
		      }
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 114:
		  {
		    /* (n, d2a) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1002;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;		
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BR] = 2.0*RDB[rea + REACTION_BR];
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 115:
		  {
		    /* (n, pd) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1002;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 116:
		  {
		    /* (n, pt) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1003;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 117:
		  {
		    /* (n, da) */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1002;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		}
	      
	      /* Decay reactions */
	      
	      switch (mt - 10000)
		{
		case 4:
		  {
		    /* Alpha decay */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 2004;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_DEC_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		case 7: 
		  {
		    /* Proton emission */
		    
		    ptr = DuplicateItem(rea);	    
		    WDB[ptr + REACTION_MT] = 20000 + 1001;
		    WDB[ptr + REACTION_RFS] = 0.0;
		    WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_DEC_BRANCH;
		    WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;
		    WDB[ptr + REACTION_BRANCH_MT] = (double)mt;
		    WDB[ptr + REACTION_RTYP2] = 0.0;
		    WDB[ptr + REACTION_RTYP3] = 0.0;
		    WDB[ptr + REACTION_RTYP4] = 0.0;
		    WDB[ptr + REACTION_RTYP5] = 0.0;
		    
		    break;
		  }
		}
	    }
	}
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /* Reset branching state pointers in duplicates */
  
  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check type and reset */

      if ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_PARTIAL)
	WDB[rea + REACTION_PTR_BRA_STATE] = NULLPTR;
      
      /* Next reaction */
      
      rea = NextItem(rea);
    }

  /***************************************************************************/

  /***** Fission yields ******************************************************/

  /* Create duplicates */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check type (NOTE: spontaneous fissiond shouldn't have branches) */

      if (((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_TRA_BRANCH) && 
	  ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_DEC_BRANCH))
	{
	  /* Pointer to first yield */
	  
	  if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
	    {
	      /* Put yield energies */

	      WDB[rea + REACTION_FISSY_IE0] = -INFTY;
	      WDB[rea + REACTION_FISSY_IE1] = RDB[yld + FISSION_YIELD_E];
	      WDB[rea + REACTION_FISSY_IE2] = INFTY;

	      /* Loop over remaining */
	      
	      yld = NextItem(yld);
	      
	      while(yld > VALID_PTR)
		{
		  /* Check type */

		  if ((long)RDB[rea + REACTION_TYPE] != REACTION_TYPE_PARTIAL)
		    Die(FUNCTION_NAME, "Invalid reaction type");

		  /* Duplicate block */
		  
		  ptr = DuplicateItem(rea);	    

		  /* Set yield pointer */
		  
		  WDB[ptr + REACTION_PTR_FISSY] = (long)yld;

		  /* Put yield energies */

		  WDB[ptr + REACTION_FISSY_IE0] = -INFTY;
		  WDB[ptr + REACTION_FISSY_IE1] = RDB[yld + FISSION_YIELD_E];
		  WDB[ptr + REACTION_FISSY_IE2] = INFTY;

		  /* Set type */
		  
		  WDB[ptr + REACTION_TYPE] = (double)REACTION_TYPE_TRA_BRANCH;
		  
		  /* Set pointer to parent */
		  
		  WDB[ptr + REACTION_PTR_BRANCH_PARENT] = (double)rea;

		  /* Set branch mt */

		  WDB[ptr + REACTION_BRANCH_MT] = RDB[rea + REACTION_MT];

		  /* Reset secondary decay modes */

		  WDB[ptr + REACTION_RTYP2] = 0.0;
		  WDB[ptr + REACTION_RTYP3] = 0.0;
		  WDB[ptr + REACTION_RTYP4] = 0.0;
		  WDB[ptr + REACTION_RTYP5] = 0.0;
		  
		  /* Next yield */
		  
		  yld = NextItem(yld);
		}
	    }
	}
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /* Put boundaries (tän voisi tehdä myös tossa edellisessä luupissa) */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Pointer to yield */
	  
      if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)
	{
	  /* Check previous yield */

	  if ((ptr = PrevItem(yld)) > VALID_PTR)
	    WDB[rea + REACTION_FISSY_IE0] = RDB[ptr + FISSION_YIELD_E];
	      
	  /* Check next yield */

	  if ((ptr = NextItem(yld)) > VALID_PTR)
	    WDB[rea + REACTION_FISSY_IE2] = RDB[ptr + FISSION_YIELD_E];
	}
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /* Check data */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Pointer to first yield */
      
      if ((yld = (long)RDB[rea + REACTION_PTR_FISSY]) > VALID_PTR)  
	{
	  /* Check that energy matches yield energy */

	  if (RDB[yld + FISSION_YIELD_E] != RDB[rea + REACTION_FISSY_IE1])
	    Die(FUNCTION_NAME, "Mismatch in energy");
	  
	  /* Check boundaries */

	  if (RDB[rea + REACTION_FISSY_IE0] >= RDB[rea + REACTION_FISSY_IE1])
	    Die(FUNCTION_NAME, "Error in boundaries");

	  if (RDB[rea + REACTION_FISSY_IE1] >= RDB[rea + REACTION_FISSY_IE2])
	    Die(FUNCTION_NAME, "Error in boundaries");

	  /* Check zeros (joissain jakaumissa on nollaenergioita) */
	  /*
	  if (RDB[rea + REACTION_FISSY_IE0] == 0.0)
	    Die(FUNCTION_NAME, "Zero yield energy");

	  if (RDB[rea + REACTION_FISSY_IE1] == 0.0)
	    Die(FUNCTION_NAME, "Zero yield energy");

	  if (RDB[rea + REACTION_FISSY_IE2] == 0.0)
	    Die(FUNCTION_NAME, "Zero yield energy");
	  */
	}

      /* Next */
      
      rea = NextItem(rea);
    }

  /***************************************************************************/

  /***** Allocate memory one-group transmutation reactions *******************/

  /* TODO: Tämä omaksi aliohjelmaksi */

  /* Transmutation and fission branching reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Get mt */

      mt = (long)RDB[rea + REACTION_MT];

      /* Allocate memory transmutation reactions */

      if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_PARTIAL) &&
	  (((mt > 15) && (mt < 50)) || ((mt > 101) && (mt < 200)) ||
	   ((mt > 599) && (mt < 900))))
	AllocValuePair(rea + REACTION_PTR_TRANSMUXS);

      /* Allocate memory for fission branch reactions */

      if ((((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_TRA_BRANCH) &&
	   ((long)RDB[rea + REACTION_PTR_FISSY] > VALID_PTR)))
	AllocValuePair(rea + REACTION_PTR_TRANSMUXS);
      
      /* Next */
      
      rea = NextItem(rea);
    }

  /* Link pointer for other branch reactions */

  rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
  while (rea > VALID_PTR)
    {
      /* Check type and allocate memory (NOTE: pointteri PRIVA-blokkiin)*/

      if (((long)RDB[rea + REACTION_TYPE] == REACTION_TYPE_TRA_BRANCH) &&
	  ((long)RDB[rea + REACTION_PTR_TRANSMUXS] < 1))
	{
	  /* Pointer to parent reaction */

	  ptr = (long)RDB[rea + REACTION_PTR_BRANCH_PARENT];
	  
	  /* Check pointer */

	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Link pointer */

	  WDB[rea + REACTION_PTR_TRANSMUXS] = 
	    RDB[ptr + REACTION_PTR_TRANSMUXS];
	}

      /* Next */
      
      rea = NextItem(rea);
    }

  /* Add reaction structure for total fission */
  
  if ((long)RDB[nuc + NUCLIDE_TYPE_FLAGS] & NUCLIDE_FLAG_FISSILE)
    {
      /* Reset Q-value */

      Q = -1.0;

      /* Loop over reactions to get Q-value */
      
      rea = (long)RDB[nuc + NUCLIDE_PTR_REA];
      while (rea > VALID_PTR)
	{
	  /* Check mt */

	  if (((long)RDB[rea + REACTION_MT] == 18) ||
	      ((long)RDB[rea + REACTION_MT] == 19))
	    {
	      /* Get value */

	      Q = RDB[rea + REACTION_Q];
	      
	      /* Break loop */

	      break;
	    }
	  
	  /* Next reaction */
	  
	  rea = NextItem(rea);
	}

      /* Check Q-value */

      if (Q < ZERO)
	Die(FUNCTION_NAME, "Nuclide %s has no fission channel",
	    GetText(nuc + NUCLIDE_PTR_NAME));

      /* Add reaction */

      rea = NewItem(nuc + NUCLIDE_PTR_TOTFISS_REA, REACTION_BLOCK_SIZE);

      /* Set mt and Q-value */

      WDB[rea + REACTION_MT] = -INFTY;
      WDB[rea + REACTION_Q] = Q;

      /* Set nuclide pointer */

      WDB[rea + REACTION_PTR_NUCLIDE] = (double)nuc;

      /* Allocate memory for one-group transmutation xs */
      
      AllocValuePair(rea + REACTION_PTR_TRANSMUXS);
    }

  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
