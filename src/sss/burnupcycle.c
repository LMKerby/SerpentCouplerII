/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : burnupcycle.c                                  */
/*                                                                           */
/* Created:       2011/05/23 (JLe)                                           */
/* Last modified: 2016/04/04 (VVa)                                           */
/* Version:       2.1.26                                                     */
/*                                                                           */
/* Description: Loop over burnup history                                     */
/*                                                                           */
/* Comments: the "step" variable indexes steps within current burnup interval*/
/*           while DATA_BURN_STEP  indexes total steps.                      */
/*                                                                           */
/*           See:                                                            */
/*                                                                           */
/*           A.E Isotalo & P.A. Aarnio "Higher Order methods for burnup      */
/*           calculations with Bateman solutions" Ann.Nucl.Energy 38 (2011)  */
/*           pp. 1987-1995.                                                  */
/*                                                                           */
/*           and                                                             */
/*                                                                           */
/*           A.E. Isotalo & P.A. Aarnio "Substep methods for burnup          */
/*           calculations with Bateman solutions" Ann.Nucl.Energy            */
/*           (Submitted)                                                     */
/*                                                                           */
/* for description of the burnup calculation methods used.                   */
/*                                                                           */
/* TODO? tätä vois ehkä  muuttaa niin että asteluvut, aliaskeleet ja tuleeko */
/*       correctori määritellään askeleen aluksi erillisessä funktiossa      */
/*       niin ettei tarkastuksia olisi siellä täällä                         */
/*                                                                           */
/*      -Voisi vaihtaa niin, että viimeinen correctori lasketaan             */
/*       paremmilla statistiikoilla, ehkä                                    */
/*      -Ton batchikoon vaihtamisen SIE:n ym. tekniikoiden kanssa vois       */
/*       siisteyden vuoksi laittaa johonkin omaan aliohjelmaansa tms.        */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "BurnupCycle:"

void BleedFluxes();

/*****************************************************************************/
void BurnupCycle()
{
  long dep, type, step, steps, pc;
  char tmpstr[MAX_STR];
  long nbatch, cycles, skip;
  
  /* avoid compiler warnings */
  nbatch=0; cycles=0; skip=0;

  /* Remove binary restart file */

  if ((long)RDB[DATA_WRITE_RESTART_FILE] == YES)
    {
      /* Get file name */

      if ((long)RDB[DATA_RESTART_WRITE_PTR_FNAME] > VALID_PTR)
	sprintf(tmpstr, "%s", GetText(DATA_RESTART_WRITE_PTR_FNAME));
      else
	sprintf(tmpstr, "%s.wrk", GetText(DATA_PTR_INPUT_FNAME));
      
      /* Remove */

      remove(tmpstr);
    }

  /* First loop is over intervals */

  dep = (long)RDB[DATA_BURN_PTR_DEP];
  while (dep > VALID_PTR)
    {
      /* Set normalization */

      SetNormalization(dep);

      /* Get step type and number of steps */

      type = (long)RDB[dep + DEP_HIS_STEP_TYPE];
      steps = (long)RDB[dep + DEP_HIS_N_STEPS];

      /* Put type */

      WDB[DATA_BURN_STEP_TYPE] = (double)type;

      /* Add final step for last interval */

      if (NextItem(dep) < VALID_PTR)
	steps++;

      /* Reprocess */

      Reprocess(dep);
      
      /* Second loop is over steps */

      for (step = 0; step < steps; step++)
	{
	  /* Initialize corrector iteration index */

	  WDB[DATA_BURN_CI_I] = (double)0;

	  /* Reset CI stopping flag */

	  if(RDB[DATA_BURN_SIE] == (double)YES)
	    WDB[DATA_BURN_CI_LAST] = (double)NO;
	  else
	    WDB[DATA_BURN_CI_LAST] = (double)YES;

	  /* Predictor-corrector -loop */

          for (pc = 0; pc < 2; pc++)
            {
             /* Set the predictor/corrector status. pc is not used directly */

              if (pc == 0)
                WDB[DATA_BURN_STEP_PC] = PREDICTOR_STEP;
              else
		{
		  WDB[DATA_BURN_STEP_PC] = CORRECTOR_STEP;
		
		  /* Signal the external program about moving */
		  /* to next burnup point */

		  if(RDB[DATA_BURN_CI_I] == 0)
		    SignalExternal(SIGUSR2);
		}

	      /* Alternate neutron population params for ci */
	      /* Tästä järkevämpi (VVa) */

	      if (((long)RDB[DATA_BURN_CI_NBATCH] > 0) &&
		  ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
		{
		  nbatch=(long)RDB[DATA_CRIT_POP];
		  cycles=(long)RDB[DATA_CRIT_CYCLES];
		  skip=(long)RDB[DATA_CRIT_SKIP];
		  
		  WDB[DATA_CRIT_POP]=RDB[DATA_BURN_CI_NBATCH];
		  WDB[DATA_CRIT_CYCLES]=RDB[DATA_BURN_CI_CYCLES];
		  WDB[DATA_CRIT_SKIP]=RDB[DATA_BURN_CI_SKIP];
		}

              /* Prepare transport cycle */
              
              PrepareTransportCycle();

              /* Transport calculation cycle if not decay step */
              
	      if ((type == DEP_STEP_DEC_STEP) || (type == DEP_STEP_DEC_TOT) ||
		  (type == DEP_STEP_ACT_STEP) || (type == DEP_STEP_ACT_TOT))
		{
		  /* Decay or activation step, set completed flag */

		  WDB[DATA_SIMULATION_COMPLETED] = (double)YES;

		  /* Calculate activities */
		  
		  CalculateActivities();

		  /* Print output */
		  
		  MatlabOutput();
		}
	      else
		{
		  /* Run transportcycle(s) */

		  do
		    {

		      /* Prepare coupled calculation iteration if needed */

		      PrepareCCIter();

		      /* Run transport cycle */

		      TransportCycle();

		      /* Iterate coupled calculation routines */

		      IterateCC();

		      /* Repeat if needed */
		    }
		  while(RDB[DATA_ITERATE] == (double)YES);
		}

	      /* Print material compositions (tää siirrettiin tuolta */
	      /* transportcyclen edestä tähän 5.12.2012 / 2.1.10)    */

	      PrintCompositions((long)RDB[DATA_BURN_STEP]);

              /***/
	      /* Restore original population parameters */
	      /* Tästä järkevämpi */
	      /***/

	      if (((long)RDB[DATA_BURN_CI_NBATCH] > 0) &&
		  ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP))
		{
		  WDB[DATA_CRIT_POP]=(double)nbatch;
		  WDB[DATA_CRIT_CYCLES]=(double)cycles;
		  WDB[DATA_CRIT_SKIP]=(double)skip;
		}

              /* Start burnup timers */
              
              ResetTimer(TIMER_BURNUP);
              StartTimer(TIMER_BURNUP);
              StartTimer(TIMER_BURNUP_TOTAL);
              
	      /* output only at predictor (corresponding to BOS and earlier) */
              
              if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP) 
                {
		  /* Write binary depletion file */

		  WriteDepFile();       

		  /* Print depletion output (print only final output in */
		  /* in decay mode) */

		  if ((((long)RDB[DATA_BURN_DECAY_CALC] == NO) &&
		       ((long)RDB[DATA_BURN_PRINT_INTERMEDIATE] == YES)) ||
		      (step == RDB[dep + DEP_HIS_N_STEPS]))
		    {
		      fprintf(out, "Writing depletion output...\n");
		      
		      PrintDepOutput();     
		      
		      fprintf(out, "OK.\n\n");
		    }
		}
              
	      /* Add step to SIE counter */

	      WDB[DATA_BURN_STEP_TOT] = RDB[DATA_BURN_STEP_TOT] + 1.0;

              /* Break here if final step */
              
              if (step == RDB[dep + DEP_HIS_N_STEPS])
		{
		  /* Stop burnup timers */
              
		  StopTimer(TIMER_BURNUP);
		  StopTimer(TIMER_BURNUP_TOTAL);

		  /* Break loop */

		  break;
		}

              /* Calculate coefficients for the fit to xs/flux/power */
              
              DepletionPolyFit(dep, step);
              
              /* Set depletion step size */
              
              SetDepStepSize(dep, step);
              
              /* Burnup calculation */
              
              BurnMaterials(dep, step);
              
	      /* Collect material compositions from MPI parallel tasks */

	      CollectBurnData();

	      /* Copy compositions to parent materials (tää siirrettiin */
	      /* burnmaterial.c:n lopusta tähän 15.7.2013 / 2.1.15 että */
	      /* kaikilla MPI taskeilla olisi käytössään sama data). */

	      SumDivCompositions();

              /* Stop burnup timers */
              
              StopTimer(TIMER_BURNUP);
              StopTimer(TIMER_BURNUP_TOTAL);
             
	      /* Add to number of predictor and corrector cycles */

	      if ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP)
		WDB[DATA_BURN_PRED_STEP] = RDB[DATA_BURN_PRED_STEP] + 1.0;
	      else
		WDB[DATA_BURN_CORR_STEP] = RDB[DATA_BURN_CORR_STEP] + 1.0;

	      /* BLEED fluxes*/
	      /*
	      BleedFluxes();
	      */

	      /* Check for iterating the corrector */

	      /* Update CI stopping criterion*/
	      if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) && ((long)RDB[DATA_BURN_SIE] == YES))
		StopCI();

	      /* If further iterations are needed */
	      if(((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) && 
		 ((long)RDB[DATA_BURN_CI_LAST] == NO) && 
                 ((long)RDB[DATA_BURN_CI_MAXI] > 1.0) )
		{
                  /* increment iteration count (indexed 0,1,2...) */

                  WDB[DATA_BURN_CI_I] = RDB[DATA_BURN_CI_I] + 1.0;

		  /* Repeat corrector */
		  pc--;

		}

              /* No corrector if: -User requested none,                      */
	      /*                  -On decay step (flux is zero, CE is exact) */
              /*                  -Activation step (26.9.2015 / JLe),        */
              /*                  -When flux is zero (from normalization),   */
              /*                   making CE exact                           */

              if (((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE) ||
                  (type == DEP_STEP_DEC_STEP) || (type == DEP_STEP_DEC_TOT) ||
                  (type == DEP_STEP_ACT_STEP) || (type == DEP_STEP_ACT_TOT) ||
                  (Mean((long)RDB[RES_TOT_NEUTRON_FLUX], 0) <= 0.0))
		{
		  /* Signal the external program about moving */
		  /* to next burnup point */

		  if(RDB[DATA_BURN_CI_I] == 0)
		    SignalExternal(SIGUSR2);

		  /* Break from PC-loop */

		  break;
		}
            }

	  /* Update cumulative burnup and time */

	  WDB[DATA_BURN_CUM_BURNTIME] = RDB[DATA_BURN_CUM_BURNTIME] 
	    + RDB[DATA_BURN_TIME_INTERVAL];
	  WDB[DATA_BURN_CUM_BURNUP] = RDB[DATA_BURN_CUM_BURNUP]
	    +  RDB[DATA_BURN_BURNUP_INTERVAL];

	  /* Update burnup step (total, not in this interval) */
       
	  WDB[DATA_BURN_STEP] = RDB[DATA_BURN_STEP] + 1.0;

	}
      
      /* Next interval */

      dep = NextItem(dep);
    }

  /* Signal externally coupled program to end calculation */

  SignalExternal(SIGTERM);

  /* Check total time */

  if (RDB[DATA_BURN_CUM_BURNTIME] == 0.0)
    Die(FUNCTION_NAME, "No burnup calculation performed");  
}

/*****************************************************************************/

/***** BLEED FLUX OUTPUT *****************************************************/
/* This outputs EOS fluxes in burnable material regions to <input>_fluxes.m  */
/* This is done to extract the averaged ones from Dufek's method.            */
/* There is no simple alternative                                            */


void BleedFluxes()
{
  long mat;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Only applicable with PC methods */

  if ((long)RDB[DATA_BURN_CORR_TYPE] == CORR_TYPE_NONE)
    return;

  sprintf(tmpstr, "%s_fluxes.m", GetText(DATA_PTR_INPUT_FNAME));

  
  if (((long)RDB[DATA_BURN_STEP] == 0) &&
      ((long)RDB[DATA_BURN_STEP_PC] == PREDICTOR_STEP))
    {

      /* first predictor => flux(1) for the initial composition */

      if ((fp = fopen(tmpstr, "w")) == NULL) 
        Die(FUNCTION_NAME, "Unable to open file for writing");

      /* Write the header info */

      fprintf(fp,"%% flux densities, [cm^-2s^-1]\n");
      fprintf(fp,"%% _flux(1) is solved with initial compositions\n");
      fprintf(fp,"%% _flux(n>1) are final predicted EOS fluxes for\n");
      fprintf(fp,"%%            the t(n-1) when t(0)=0\n");
      fprintf(fp,"%% These values are valid results for Dufek's method");
      fprintf(fp," ONLY\n\n");
      
      /* Write "definitions" for each material */

      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {

          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            fprintf(fp,"%s_flux=zeros(%ld,1);\n",GetText(mat+MATERIAL_PTR_NAME),
                    (long)RDB[DATA_BURN_TOT_STEPS] + 1);

          mat = NextItem(mat);
        }
      fprintf(fp,"\n");

      /* write fluxes */
            
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            fprintf(fp,"%s_flux(%ld)=%.12e;\n",GetText(mat+MATERIAL_PTR_NAME),
                    (long)RDB[DATA_BURN_STEP] + 1,
                    WDB[mat + MATERIAL_BURN_FLUX_BOS]);
          
          mat = NextItem(mat);
        }
      fprintf(fp,"\n");
      

    }
  else if (((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) &&
           ((long)WDB[DATA_BURN_CI_LAST] == YES))
    {
      /* last corector => flux(step+2) */

      if ((fp = fopen(tmpstr, "a")) == NULL) 
        Die(FUNCTION_NAME, "Unable to open file for appending");


      /* write fluxes */
      
      
      mat = (long)RDB[DATA_PTR_M0];
      while (mat > VALID_PTR)
        {
          if ((long)RDB[mat + MATERIAL_OPTIONS] & OPT_BURN_MAT)
            fprintf(fp,"%s_flux(%ld)=%.12e;\n",GetText(mat+MATERIAL_PTR_NAME),
                    (long)RDB[DATA_BURN_STEP] + 2,
                    WDB[mat + MATERIAL_BURN_FLUX_EOS]);
          
          mat = NextItem(mat);
        }
      fprintf(fp,"\n");

    }
  else
    {
      return;
    }
  
  
  
   fclose(fp);


}

