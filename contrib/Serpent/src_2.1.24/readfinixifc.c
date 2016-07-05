/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readfinixifc.c                                 */
/*                                                                           */
/* Created:       2015/03/11 (VVa)                                           */
/* Last modified: 2015/06/25 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Reads initial guess from interface file for FINIX            */
/*                                                                           */
/* Comments:   -If the interface file does not correspond to the FINIX pin   */
/*              definitions, it won't be read                                */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#ifdef FINIX

#include "./FINIX/defaults.h"
#include "./FINIX/FINIX.h"

#define FUNCTION_NAME "ReadFinixIFC:"

/*****************************************************************************/

void ReadFinixIFC()
{
  long fib, i, j, k;
  long dummy, pins;
  long nz, na, nr;
  double zmin, zmax, amin, amax, rmin, rmax;
  double T, rh, rc;
  double dbl0, dbl1;
  char tmpstr[MAX_STR];
  int *nnodes;
  double **Tptr;
  double **rhptr;
  double **rcptr;
  FILE *fp;

  /* Print interface file name */

  sprintf(tmpstr,"./Finix0.ifc");

  /* Try to open the interface file for writing */
  /* If file cannot be read, just return        */
  /* HZP initial conditions will be used        */

  if ((fp = fopen(tmpstr, "r")) == NULL)
    return;

  fprintf(out, "Reading initial conditions for FINIX pins\n");

  /* Read interface type */
  
  if (fscanf(fp, "%ld", &dummy) == EOF)
    return;

  /* Check that interface type is 6 */

  if (dummy != 6)
    return;

  /* Get output file name (not checked) */

  if (fscanf(fp, "%s", tmpstr) == EOF)
    return;

  /* Get number of pins in interface */

  if (fscanf(fp, "%ld", &dummy) == EOF)
    return;

  /* Count number of finix pins */
  
  fib = (long)RDB[DATA_PTR_FIN0];  
  pins = 0;

  while (fib > VALID_PTR)
    {
      /* Increment pin counter */

      pins++;

      /* Get next finix block*/

      fib = NextItem(fib);
    }

  /* Check that the number of pins in interface file corresponds */
  /* to the number of FINIX-pins in memory                       */

  if (dummy != pins)
    return;

  /* Read fuel rod information */

  for (i = 0; i < pins; i++)
    {
      /* Read pin universe */
	      
      if (fscanf(fp, "%s", tmpstr) == EOF)
	Die(FUNCTION_NAME, "fscanf error");

      /* Find finix block from memory */

      fib = (long)RDB[DATA_PTR_FIN0];
      CheckPointer(FUNCTION_NAME, "(fib)", DATA_ARRAY, fib);

      /* Find correct pin */

      while (fib > VALID_PTR)
	{
	  /* Compare ifc pin name to FINIX pin name */

	  if(!strcmp(tmpstr, GetText(fib + FINIX_PTR_UNI_NAME)))
	    break;

	  /* Next FINIX pin */

	  fib = NextItem(fib);	 

	}

      /* Check if found */

      if (fib < VALID_PTR)
	{
	  fprintf(out, "Interface file contains unknown pin %s\n", tmpstr);
	  return;
	}

      /* Get FINIX pointers */

      Tptr = (double **)((long)RDB[fib + FINIX_PTR_T]);
      rhptr = (double**)((long)RDB[fib + FINIX_PTR_R]);
      rcptr = (double**)((long)RDB[fib + FINIX_PTR_R_COLD]);
      nnodes = (int*)((long)RDB[fib + FINIX_PTR_NNODES]);

      /* Read power output limits */

      if (fscanf(fp, "%ld %lf %lf %ld %lf %lf %ld %lf %lf",
		&nz, &zmin, &zmax, &na, &amin, &amax,
		&nr, &rmin, &rmax) == EOF)
	return;

      /* Check axial power tally size and limits */

      if (nz != (long)RDB[fib + FINIX_NZ_POW])
	return;

      if (zmin != RDB[fib + FINIX_ZMIN])
	return;

      if (zmax != RDB[fib + FINIX_ZMAX])
	return;

      /* Check angular power tally size and limits */

      if (na != 1)
	return;

      if (amin != 0.0)
	return;

      if (amax != 360.0)
	return;

      /* Check radial power tally size and limits */

      if (nr != (long)RDB[fib + FINIX_NR_POW])
	return;

      if (rmin != RDB[fib + FINIX_RMIN])
	return;

      if (rmax != RDB[fib + FINIX_RMAX])
	return;

      /* Read flux output limits (not checked) */

      if (fscanf(fp, "%ld %lf %lf %ld %lf %lf %ld %lf %lf",
		&nz, &zmin, &zmax, &na, &amin, &amax,
		&nr, &rmin, &rmax) == EOF)
	return;

      /* Read flux energy limits (not checked) */

      if (fscanf(fp, "%lf %lf", &zmin, &zmax) == EOF)
	return;

      /* Read number of axial zones */

      if (fscanf(fp, "%ld", &nz) == EOF)
	return;

      /* Check number of axial zones */

      if (nz != nnodes[0])
	return;

      /* Loop over axial zones */

      for (j = 0; j < nnodes[0]; j++)
	{
	  /* Calculate zmin and zmax of this segment */
      
	  zmin = RDB[fib + FINIX_ZMIN]+
	    j*(RDB[fib + FINIX_ZMAX]-
	       RDB[fib + FINIX_ZMIN])/RDB[fib + FINIX_NZ];

	  zmax = RDB[fib + FINIX_ZMIN]+
	    (j+1)*(RDB[fib + FINIX_ZMAX]-
		   RDB[fib + FINIX_ZMIN])/RDB[fib + FINIX_NZ];      

	  /* Read axial limits (not checked) */
	  
	  if (fscanf(fp, "%lf %lf", &dbl0, &dbl1) == EOF)
	    return;

	  /* Check that axial limits for this segment are correct */

	  if (fabs(zmin - dbl0) > 1E-5)
	    return;

	  if (fabs(zmax - dbl1) > 1E-5)
	    return;
	 
	  /* Read number of angular zones */

	  if (fscanf(fp, "%ld", &dummy) == EOF)
	    return;

	  /* Check number of angular zones */

	  if (dummy != 1)
	    return;
	  
	  /* Get angular limits of the zone */

	  if (fscanf(fp, "%lf %lf", &dbl0, &dbl1) == EOF)
	    return;

	  /* Check that angular limits for this segment are correct */

	  if (fabs(0 - dbl0) > 1E-5)
	    return;

	  if (fabs(360.0 - dbl1) > 1E-5)
	    return;

	  /* Read number of radial points */

	  if (fscanf(fp, "%ld", &dummy) == EOF)
	    return;

	  /* Calculate number of supposed radial points */

	  if (rcptr[j][0] > 0.0)
	    {
	      /* Rod with central hole */

	      nr = nnodes[1] + nnodes[2] + 1;

	    }
	  else
	    {

	      /* Rod without central hole */

	      nr = nnodes[1] + nnodes[2];

	    }

	  /* Check number of radial points in file */

	  if (nr != dummy)
	    {
	      fprintf(out, 
		      "Number of radial points doesn't match %ld vs %ld\n", 
		      nr, dummy);
	      return;
	    }

	  /* Read and discard central hole data */

	  if (rcptr[j][0] > 0.0)
	    if (fscanf(fp, "%lf %lf %lf", &rc, &rh, &T) == EOF)
	      return;
	    

	  /* Loop over radial points and get hot radius and temperature */

	  for (k = 0; k < nnodes[1] + nnodes[2]; k++)
	    {

	      /* Read cold radius, hot radius and temperature */

	      if (fscanf(fp, "%lf %lf %lf", &rc, &rh, &T) == EOF)
		return;

	      /* Check cold radius */
	      /* TODO: Ehkä joku suhteellinen vertailu */

	      if(fabs(rcptr[j][k]*100 - rc) > 1E-3)
		return;

	      /* Store hotradius */

	      rhptr[j][k] = rh/100.0;

	      /* Store temperature */

	      Tptr[j][k] = T;
	      
	    }
	  
	}
    }
      
  fprintf(out, "All read\n");

  return;

}	  

#endif

/*****************************************************************************/
