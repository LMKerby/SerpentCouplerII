/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readsourcefile.c                               */
/*                                                                           */
/* Created:       2012/03/30 (JLe)                                           */
/* Last modified: 2012/11/02 (JLe)                                           */
/* Version:       2.1.10                                                     */
/*                                                                           */
/* Description: Reads source distribution from file                          */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadSourceFile:"

/*****************************************************************************/

void ReadSourceFile(long src, double *x, double *y, double *z, double *u, 
		    double *v, double *w, double *E, double *wgt, double *t)
{
  long idx, sz, ptr, pos, eof;
  FILE *fp;

  /* Check source pointer */

  CheckPointer(FUNCTION_NAME, "(src)", DATA_ARRAY, src);

  /* Get buffer size */

  sz = (long)RDB[src + SRC_READ_BUF_SZ];

  /* Set OpenMP barrier */

#ifdef OPEN_MP
#pragma omp critical
#endif
  {
    /* Get index to next neutron in buffer */
    
    if ((idx = (long)RDB[src + SRC_READ_BUF_IDX]) > sz - 1)
      {
	/*********************************************************************/

	/***** Fill buffer with source points ********************************/

	/* Open file for reading */
      
	if ((fp = fopen(GetText(src + SRC_READ_PTR_FILE), "r")) == NULL)
	  Error(src, "Unable to open source file \"%s\"", 
		GetText(src + SRC_READ_PTR_FILE));

	/* Get pointer to buffer */
	
	ptr = (long)RDB[src + SRC_READ_PTR_BUF];
	CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	
	/* Reset buffer index */
	
	idx = -1;
	
	/* Loop until buffer is full */
	
	do
	  {
	    /* Get position */
	    
	    pos = (long)RDB[src + SRC_READ_FILE_POS];
	    
	    /* Seek to position */
	    
	    fseek(fp, pos, SEEK_SET);
	    
	    /* Loop over data */
	    
	    while ((eof = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
				 x, y, z, u, v, w, E, wgt, t)) != EOF)
	      {
		/* Update index */
		
		idx++;
		
		/* Get position in file */
		
		pos = ftell(fp);

		/* Renormalization for type 2 */

		if ((long)RDB[src + SRC_READ_FILE_TYPE] == 
		    SRC_FILE_TYPE_S1_RENORM)
		  {
		    *wgt = 1.0;
		    *t = 0.0;
		  }
	
		/* Put data */
		
		WDB[ptr + SRC_BUF_X] = *x;
		WDB[ptr + SRC_BUF_Y] = *y;
		WDB[ptr + SRC_BUF_Z] = *z;
		WDB[ptr + SRC_BUF_U] = *u;
		WDB[ptr + SRC_BUF_V] = *v;
		WDB[ptr + SRC_BUF_W] = *w;
		WDB[ptr + SRC_BUF_E] = *E;
		WDB[ptr + SRC_BUF_WGT] = *wgt;
		WDB[ptr + SRC_BUF_T] = *t;
		
		/* Check if buffer is full */
		
		if (idx == sz - 1)
		  break;
		
		/* Update pointer */
		
		ptr = ptr + SRC_BUF_BLOCK_SIZE;	    
	      }
	    
	    /* Set or reset file position */
	    
	    if (eof == EOF)
	      WDB[src + SRC_READ_FILE_POS] = 0.0;
	    else
	      WDB[src + SRC_READ_FILE_POS] = (double)pos;
	  }
	while (idx < sz - 1);
	
	/* Reset buffer index */
	
	idx = 0;
	WDB[src + SRC_READ_BUF_IDX] = (double)idx;
	
	/* Close file */
	
	fclose(fp);
     
	/*********************************************************************/
      }
  }

  /***************************************************************************/

  /***** Get data from buffer ************************************************/

  /* Check index */

  if ((idx < 0) || (idx > sz - 1))
    Die(FUNCTION_NAME, "Error in buffer index");
  
  /* Get pointer to buffer */

  ptr = (long)RDB[src + SRC_READ_PTR_BUF];
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get pointer to data */
  
  ptr = ptr + idx*SRC_BUF_BLOCK_SIZE;
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Get data */
  
  *x = RDB[ptr + SRC_BUF_X];
  *y = RDB[ptr + SRC_BUF_Y];
  *z = RDB[ptr + SRC_BUF_Z];
  *u = RDB[ptr + SRC_BUF_U];
  *v = RDB[ptr + SRC_BUF_V];
  *w = RDB[ptr + SRC_BUF_W];
  *E = RDB[ptr + SRC_BUF_E];
  *wgt = RDB[ptr + SRC_BUF_WGT];
  *t = RDB[ptr + SRC_BUF_T];

  /* Update buffer index */

#ifdef OPEN_MP
#pragma omp atomic
#endif

  WDB[src + SRC_READ_BUF_IDX] += 1.0;
  
  /***************************************************************************/

  /***** Check values and find position in geometry **************************/

  /* Check energy */

  if ((*E < ZERO) || (*E > 1E+6))
    Error(src, "Invalid particle energy read from source file");

  /* Check weight */

  if ((*wgt < ZERO) || (*wgt > 1E+6))
    Error(src, "Invalid particle weight read from source file");

  /* Check direction cosines */

  if (fabs((*u)*(*u) + (*v)*(*v) + (*w)*(*w) - 1.0) > 1E-5)
    Error(src, "Direction cosines read from source file are not normalized");

  /***************************************************************************/
}

/*****************************************************************************/
