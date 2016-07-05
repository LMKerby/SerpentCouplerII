/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writefinixifc.c                                */
/*                                                                           */
/* Created:       2013/03/11 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Creates an IFC-template for the ReadInterface()              */
/*              for the initialization of the FINIX interface                */
/*                                                                           */
/* Comments:   -Tän vois periaatteessa kirjoittaa suoraan muistiinkin        */
/*              mutta lienee järkevintä luoda kaikki rajapinnat vasta        */
/*              ReadInterface():ssa                                          */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "WriteFinixIFC:"

/*****************************************************************************/

void WriteFinixIFC()
{
  long fib, axi, i, j, k;
  double zmin, zmax, nt;
  char tmpstr[MAX_STR];
  int *nnodes;
  double **T;
  double **r;
  double **r_cold;
  FILE *fp;

  /* Print interface file name */

  sprintf(tmpstr,"./Finix.ifc");

  /* Try to open the interface file for writing */

  if ((fp = fopen(tmpstr, "w")) == NULL)
    Die(FUNCTION_NAME, "Could not open file \"%s\" for writing", tmpstr);

  fib = (long)RDB[DATA_PTR_FIN0];  

  /* Write the interface file header */

  /* Output filename*/
  sprintf(tmpstr,"Finixifcout.m");

  /* Interface type, outputfile, number of pins */

  fprintf(fp,"6 %s %ld\n", tmpstr, (long)RDB[fib + FINIX_N_PIN]);

  /* Loop over pins */

  while(fib > VALID_PTR)
    {
      /* Get the Finix pointers */

      nnodes = (int*)((long)RDB[fib + FINIX_PTR_NNODES]);
      T = (double **)((long)RDB[fib + FINIX_PTR_T]);
      r = (double**)((long)RDB[fib + FINIX_PTR_R]);
      r_cold = (double**)((long)RDB[fib + FINIX_PTR_R_COLD]);


      /* Pin universe */

      fprintf(fp,"%s\n",GetText(fib + FINIX_PTR_UNI_NAME));

      /* Power tally */

      fprintf(fp,"%ld %f %f 1 0 360 %ld %f %f\n", (long)RDB[fib + FINIX_NZ_POW], RDB[fib + FINIX_ZMIN], 
	      RDB[fib + FINIX_ZMAX], (long)RDB[fib + FINIX_NR_POW], r_cold[0][0]*100, r_cold[0][nnodes[1]-1]*100);

      /* Fast flux tally (only one radial bin) */

      fprintf(fp,"%ld %f %f 1 0 360 1 %f %f 1 15\n", (long)RDB[fib + FINIX_NZ_POW], RDB[fib + FINIX_ZMIN], 
	      RDB[fib + FINIX_ZMAX], r_cold[0][0]*100, r_cold[0][nnodes[1]+nnodes[2]-1]*100);

      /* Get number of timesteps */

      if((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_SRC)
	nt = (double)(RDB[DATA_DYN_NB]*1);
      else
	nt = (double)1;

      /* Write interface data for each timestep */

      for(k=0; k < nt; k++)
	{

	  /* Write number of axial zones */

	  fprintf(fp,"%ld\n",(long)RDB[fib + FINIX_NZ]);

	  axi = (long)RDB[fib + FINIX_PTR_AX];

	  /* Loop over the axial segments */

	  for(i=0; i<nnodes[0]; i++)
	    {
	      /* Get zmin and zmax of this segment */
      
	      zmin = RDB[fib + FINIX_ZMIN]+
		i*(RDB[fib + FINIX_ZMAX]-
		   RDB[fib + FINIX_ZMIN])/RDB[fib + FINIX_NZ];
	      zmax = RDB[fib + FINIX_ZMIN]+
		(i+1)*(RDB[fib + FINIX_ZMAX]-
		       RDB[fib + FINIX_ZMIN])/RDB[fib + FINIX_NZ];      

	      /* Put limits of axial segment (only 1 angular zone) */

	      fprintf(fp,"%f %f 1\n", zmin, zmax);

	      /* Put limits of angular zone and number of radial nodes */

	      if (r_cold[i][0] > 0.0)
		fprintf(fp,"0 360 %ld\n", (long)(nnodes[1]+nnodes[2]+1));
	      else
		fprintf(fp,"0 360 %ld\n", (long)(nnodes[1]+nnodes[2]));

	      /* Print out radial node information */

	      if (r_cold[i][0] > 0.0)
		{
		  /* Central hole */

		  fprintf(fp,"%E %E %E\n", 0.0, 0.0, T[i][0]);
		
		}

	      for(j=0; j < nnodes[1] + nnodes[2]; j++)
		fprintf(fp,"%E %E %E\n", r_cold[i][j]*100, r[i][j]*100, T[i][j]);
      
	      axi = NextItem(axi);
	    }
	}
      fib = NextItem(fib);
    }

  fclose(fp);

}	  

#endif

/*****************************************************************************/
