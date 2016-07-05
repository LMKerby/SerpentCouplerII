/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : writesourcefile.c                              */
/*                                                                           */
/* Created:       2012/10/19 (JLe)                                           */
/* Last modified: 2015/01/20 (VVa)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Writes source distribution to file                           */
/*                                                                           */
/* Comments: - The distribution is actually stored into a buffer, and        */
/*             written only after the buffer is full.                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "WriteSourceFile:"

/*****************************************************************************/

void WriteSourceFile(long det, double x, double y, double z, double u, 
		     double v, double w, double E, double wgt, double t, 
		     double flx, long id)
{
  long loc0, ptr, idx, sz;
  double spd, P;
  FILE *fp;
  
  /* Check pointer */

  if (det > VALID_PTR)
    {
      /***********************************************************************/

      /***** Write point in buffer *******************************************/

      /* Get pointer to buffer */
	
      if ((loc0 = (long)RDB[det + DET_WRITE_PTR_BUF]) < VALID_PTR)
	return;
      
      /* Calculate probability */

      if(flx < 0.0)
	{
	  /* Store all points with probability DET_WRITE_PROB  */
	  /* Surface crossing or Fission points                */

	  P = RDB[det + DET_WRITE_PROB];
	}
      else
	{
	  /* storing porbability weighted with neutron density */
	  /* weighted with 1/v*flx    (flx = 1/totxs)          */

	  spd = Speed(PARTICLE_TYPE_NEUTRON, E);

	  P = RDB[det + DET_WRITE_PROB]*wgt*
	    Speed(PARTICLE_TYPE_NEUTRON, RDB[DATA_NEUTRON_EMIN])/spd*
	    MinXS(PARTICLE_TYPE_NEUTRON, spd, id)*flx;
	}

      if(P > 1)
	Warn(FUNCTION_NAME, "P larger than 1 (%E)", P);

      /* Rejection */

      if (RandF(id) > P)
	return;

      /* Get buffer size */

      sz = (long)RDB[det + DET_WRITE_BUF_SZ];

      /* Set OpenMP barrier */

#ifdef OPEN_MP
#pragma omp critical
#endif
      {
	/* Get index */

	idx = (long)RDB[det + DET_WRITE_BUF_IDX];
	CheckValue(FUNCTION_NAME, "idx", "", idx, 0, sz - 1);

	/* Get pointer */

	ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;	    
	CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	/* Put data */
		
	WDB[ptr + SRC_BUF_X] = x;
	WDB[ptr + SRC_BUF_Y] = y;
	WDB[ptr + SRC_BUF_Z] = z;
	WDB[ptr + SRC_BUF_U] = u;
	WDB[ptr + SRC_BUF_V] = v;
	WDB[ptr + SRC_BUF_W] = w;
	WDB[ptr + SRC_BUF_E] = E;
	WDB[ptr + SRC_BUF_WGT] = wgt;
	WDB[ptr + SRC_BUF_T] = t;

	/* Update index */

	idx++;

	/* Check if buffer is full */

	if (idx == sz)
	  {
	    /* Open file */

	    if ((fp = fopen(GetText(det + DET_WRITE_PTR_FILE), "a")) == NULL)
	      Error(det, "Unable to open file \"%s\" for writing", 
		    GetText(det + DET_WRITE_PTR_FILE));

	    /* Loop over buffer */

	    for (idx = 0; idx < sz; idx++)
	      {
		/* Get pointer */

		ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;	    
		CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

		/* Write data */

		fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
		fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
		fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
		fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
		fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
		fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
		fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
		fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
		fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
		fprintf(fp, "\n");
	      }

	    /* Close file */

	    fclose(fp);
	    
	    /* Reset index */

	    idx = 0;
	  }

	/* Put new index */

	WDB[det + DET_WRITE_BUF_IDX] = (double)idx;
      }

      /***********************************************************************/
    }
  else
    {
      /***********************************************************************/

      /***** Dump all buffers ************************************************/

      /* Check OpenMP thread number */

      if (OMP_THREAD_NUM != 0)
	Die(FUNCTION_NAME, "Called from an OpenMP parallel loop");

      /* Loop over detectors */

      det = (long)RDB[DATA_PTR_DET0];
      while (det > VALID_PTR)
	{
	  /* Check pointer to write buffer */
	
	  if ((loc0 = (long)RDB[det + DET_WRITE_PTR_BUF]) > VALID_PTR)
	    {
	      /* Get buffer size */

	      sz = (long)RDB[det + DET_WRITE_BUF_IDX];

	      /* Open file */

	      if ((fp = fopen(GetText(det + DET_WRITE_PTR_FILE), "a")) == NULL)
		Error(det, "Unable to open file \"%s\" for writing", 
		      GetText(det + DET_WRITE_PTR_FILE));	      
	      
	      /* Loop over buffer */

	      for (idx = 0; idx < sz; idx++)
		{
		  /* Get pointer */
		  
		  ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;	    
		  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
		  
		  /* Write data */
		  
		  fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
		  fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
		  fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
		  fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
		  fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
		  fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
		  fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
		  fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
		  fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
		  fprintf(fp, "\n");
		}

	      /* Close file */

	      fclose(fp);

	      /* Reset index */
	      
	      WDB[det + DET_WRITE_BUF_IDX] = 0.0;
	    }

	  /* Next source */

	 det = NextItem(det);
	}

      /* Criticality source detector */

      if ((det = (long)RDB[DATA_PTR_CRIT_SRC_DET]) > VALID_PTR)
	{
	  /* Get pointer to write buffer */

	  loc0 = (long)RDB[det + DET_WRITE_PTR_BUF];
	  CheckPointer(FUNCTION_NAME, "(loc0)", DATA_ARRAY, loc0);

	  /* Get buffer size */
	  
	  sz = (long)RDB[det + DET_WRITE_BUF_IDX];
	  
	  /* Open file */
	  
	  if ((fp = fopen(GetText(det + DET_WRITE_PTR_FILE), "a")) == NULL)
	    Error(det, "Unable to open file \"%s\" for writing", 
		  GetText(det + DET_WRITE_PTR_FILE));	      
	  
	  /* Loop over buffer */

	  for (idx = 0; idx < sz; idx++)
	    {
	      /* Get pointer */
	      
	      ptr = loc0 + idx*SRC_BUF_BLOCK_SIZE;	    
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	      /* Write data */
	      
	      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_X]);
	      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Y]);
	      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_Z]);
	      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_U]);
	      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_V]);
	      fprintf(fp, "%12.5E ", RDB[ptr + SRC_BUF_W]);
	      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_E]);
	      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_WGT]);
	      fprintf(fp, "%11.5E ", RDB[ptr + SRC_BUF_T]);
	      fprintf(fp, "\n");
	    }

	  /* Close file */
	  
	  fclose(fp);
	  
	  /* Reset index */
	  
	  WDB[det + DET_WRITE_BUF_IDX] = 0.0;
	}
      
      /***********************************************************************/
    }
}

/*****************************************************************************/
