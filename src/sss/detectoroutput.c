#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : detectoroutput.c                               */
/*                                                                           */
/* Created:       2011/03/03 (JLe)                                           */
/* Last modified: 2015/07/03 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Prints detector output                                       */
/*                                                                           */
/* Comments: - Time binejä ei pidä sallia classic lookin kanssa              */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "DetectorOutput:"

/*****************************************************************************/

void DetectorOutput()
{
  long det0, det1, erg, ptr, n0, n1, n2, n, m, i, j;
  long ebins0, ubins0, cbins0, mbins0, lbins0, rbins0, zbins0, ybins0, xbins0;
  long tbins0, ebins1, ubins1, cbins1, mbins1, lbins1, rbins1, zbins1, ybins1;
  long xbins1, tbins1, eb0, ub0, cb0, mb0, lb0, rb0, zb0, yb0, xb0, tb0, eb1;
  long ub1, cb1, mb1, lb1, rb1, zb1, yb1, xb1, tb1, idx0, idx1;
  double min0, max0, min1, max1, min2, max2, x, y, x0, y0, pitch;
  FILE *fp;
  char tmpstr[MAX_STR];

  /* Check detector definitions */

  if ((long)RDB[DATA_PTR_DET0] < VALID_PTR)
    return;

  /* Check mpi task */

  if (mpiid > 0)
    return;

  /* Check if in active cycles */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return;

#ifdef STAB_BURN

  /* Open file for writing. Modified to save pred. and corr. separately (AIs)*/

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP)
    sprintf(tmpstr, "%s_det%ldc%d.m", GetText(DATA_PTR_INPUT_FNAME),
	    (long)RDB[DATA_BURN_STEP], (int)RDB[DATA_BURN_CI_I]);
  else
    sprintf(tmpstr, "%s_det%ldp.m", GetText(DATA_PTR_INPUT_FNAME),
	    (long)RDB[DATA_BURN_STEP]);

#else

  /* Check corrector step */

  if ((long)RDB[DATA_BURN_STEP_PC] == CORRECTOR_STEP) 
    return;

  /* Check branch index */

  if ((long)RDB[DATA_COEF_CALC_IDX] < 0)
    sprintf(tmpstr, "%s_det%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	    (long)RDB[DATA_BURN_STEP]);
  else 
    sprintf(tmpstr, "%s_det%ldb%ld.m", GetText(DATA_PTR_INPUT_FNAME),
	    (long)RDB[DATA_COEF_CALC_BU_IDX], 
	    (long)RDB[DATA_COEF_CALC_IDX]);

#endif
  
  if ((fp = fopen(tmpstr, "w")) == NULL) 
    Die(FUNCTION_NAME, "Unable to open file for writing");

  /* Loop over detectors */
  
  det0 = (long)RDB[DATA_PTR_DET0];
  while(det0 > VALID_PTR)
    {
      /* Get number of bins */
	  
      ebins0 = (long)RDB[det0 + DET_N_EBINS];
      ubins0 = (long)RDB[det0 + DET_N_UBINS];
      cbins0 = (long)RDB[det0 + DET_N_CBINS];
      mbins0 = (long)RDB[det0 + DET_N_MBINS];
      lbins0 = (long)RDB[det0 + DET_N_LBINS];
      rbins0 = (long)RDB[det0 + DET_N_RBINS];
      tbins0 = (long)RDB[det0 + DET_N_TBINS];

      /* Mesh bins */

      if ((ptr = (long)RDB[det0 + DET_PTR_MESH]) > VALID_PTR)
	{
	  xbins0 = (long)RDB[ptr + MESH_N0];
	  ybins0 = (long)RDB[ptr + MESH_N1];
	  zbins0 = (long)RDB[ptr + MESH_N2];
	}
      else
	{
	  xbins0 = 1;
	  ybins0 = 1;
	  zbins0 = 1;
	}
      
      /* Check mode */

      if (tbins0 == 1)
	{
	  /*******************************************************************/

	  /***** Serpent 1 type output (no time bins) ************************/

	  fprintf(fp, "\n");

	  fprintf(fp, "DET%s = [\n", GetText(det0 + DET_PTR_NAME));

	  /* Time bins are not used */

	  tb0 = 0;
	  tb1 = 0;

	  /* Pointer to adjoint detector */
	  
	  if ((det1 = (long)RDB[det0 + DET_PTR_ADJOINT]) < VALID_PTR)
	    {
	      /* Pointer to statistics */
	  
	      ptr = (long)RDB[det0 + DET_PTR_STAT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	      
	      /* Loop over bins */

	      n0 = 0;

	      for (eb0 = 0; eb0 < ebins0; eb0++)
	      for (ub0 = 0; ub0 < ubins0; ub0++)
	      for (cb0 = 0; cb0 < cbins0; cb0++)
	      for (mb0 = 0; mb0 < mbins0; mb0++)
	      for (lb0 = 0; lb0 < lbins0; lb0++)
	      for (rb0 = 0; rb0 < rbins0; rb0++)
	      for (zb0 = 0; zb0 < zbins0; zb0++)
	      for (yb0 = 0; yb0 < ybins0; yb0++)
	      for (xb0 = 0; xb0 < xbins0; xb0++)
		{
		  /* Print indexes */
		  
		  fprintf(fp, "%5ld ", n0 + 1);
		  fprintf(fp, "%4ld ", eb0 + 1);
		  fprintf(fp, "%4ld ", ub0 + 1);
		  fprintf(fp, "%4ld ", cb0 + 1);
		  fprintf(fp, "%4ld ", mb0 + 1);
		  fprintf(fp, "%4ld ", lb0 + 1);
		  fprintf(fp, "%4ld ", rb0 + 1);
		  fprintf(fp, "%4ld ", zb0 + 1);
		  fprintf(fp, "%4ld ", yb0 + 1);
		  fprintf(fp, "%4ld ", xb0 + 1);

		  /* Get index */

		  idx0 = DetIdx(det0, eb0, ub0, cb0, mb0, lb0, zb0, yb0, xb0, 
				tb0);

		  /* Print mean */

		  fprintf(fp, "%12.5E ", Mean(ptr, idx0, rb0));
		  
		  /* Print relative statistical error */
		  
		  fprintf(fp, "%7.5f ", RelErr(ptr, idx0, rb0));
		  
		  /* Print newline */

		  fprintf(fp, "\n");

		  /* Update index */
		  
		  n0++;
		}
	    }
	  else
	    {
	      /* Get number of bins for second detector */
	      
	      ebins1 = (long)RDB[det1 + DET_N_EBINS];
	      ubins1 = (long)RDB[det1 + DET_N_UBINS];
	      cbins1 = (long)RDB[det1 + DET_N_CBINS];
	      mbins1 = (long)RDB[det1 + DET_N_MBINS];
	      lbins1 = (long)RDB[det1 + DET_N_LBINS];
	      rbins1 = (long)RDB[det1 + DET_N_RBINS];
	      tbins1 = (long)RDB[det1 + DET_N_TBINS];

	      /* Mesh bins */
	      
	      if ((ptr = (long)RDB[det1 + DET_PTR_MESH]) > VALID_PTR)
		{
		  xbins1 = (long)RDB[ptr + MESH_N0];
		  ybins1 = (long)RDB[ptr + MESH_N1];
		  zbins1 = (long)RDB[ptr + MESH_N2];
		}
	      else
		{
		  xbins1 = 1;
		  ybins1 = 1;
		  zbins1 = 1;
		}

	      /* Pointer to statistics */
	  
	      ptr = (long)RDB[det0 + DET_PTR_STAT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
	  
	      /* Loop over bins */

	      n0 = 0;

	      for (eb0 = 0; eb0 < ebins0; eb0++)
	      for (ub0 = 0; ub0 < ubins0; ub0++)
	      for (cb0 = 0; cb0 < cbins0; cb0++)
	      for (mb0 = 0; mb0 < mbins0; mb0++)
	      for (lb0 = 0; lb0 < lbins0; lb0++)
	      for (rb0 = 0; rb0 < rbins0; rb0++)
	      for (zb0 = 0; zb0 < zbins0; zb0++)
	      for (yb0 = 0; yb0 < ybins0; yb0++)
	      for (xb0 = 0; xb0 < xbins0; xb0++)
		{
		  /* Loop over second bins */

		  n1 = 0;

		  for (eb1 = 0; eb1 < ebins1; eb1++)
		  for (ub1 = 0; ub1 < ubins1; ub1++)
		  for (cb1 = 0; cb1 < cbins1; cb1++)
		  for (mb1 = 0; mb1 < mbins1; mb1++)
		  for (lb1 = 0; lb1 < lbins1; lb1++)
		  for (rb1 = 0; rb1 < rbins1; rb1++)
		  for (zb1 = 0; zb1 < zbins1; zb1++)
		  for (yb1 = 0; yb1 < ybins1; yb1++)
		  for (xb1 = 0; xb1 < xbins1; xb1++)
		      {
			/* Print indexes */
			
			fprintf(fp, "%5ld ", n0 + 1);
			fprintf(fp, "%4ld ", eb0 + 1);
			fprintf(fp, "%4ld ", ub0 + 1);
			fprintf(fp, "%4ld ", cb0 + 1);
			fprintf(fp, "%4ld ", mb0 + 1);
			fprintf(fp, "%4ld ", lb0 + 1);
			fprintf(fp, "%4ld ", rb0 + 1);
			fprintf(fp, "%4ld ", zb0 + 1);
			fprintf(fp, "%4ld ", yb0 + 1);
			fprintf(fp, "%4ld ", xb0 + 1);

			fprintf(fp, "%5ld ", n1 + 1);
			fprintf(fp, "%4ld ", eb1 + 1);
			fprintf(fp, "%4ld ", ub1 + 1);
			fprintf(fp, "%4ld ", cb1 + 1);
			fprintf(fp, "%4ld ", mb1 + 1);
			fprintf(fp, "%4ld ", lb1 + 1);
			fprintf(fp, "%4ld ", rb1 + 1);
			fprintf(fp, "%4ld ", zb1 + 1);
			fprintf(fp, "%4ld ", yb1 + 1);
			fprintf(fp, "%4ld ", xb1 + 1);
			
			/* Get indexes */

			idx0 = DetIdx(det0, eb0, ub0, cb0, mb0, lb0, zb0, yb0, 
				      xb0, tb0);
			
			idx1 = DetIdx(det1, eb1, ub1, cb1, mb1, lb1, zb1, yb1, 
				      xb1, tb1);

			/* Print mean */

			fprintf(fp, "%12.5E ", Mean(ptr, idx0, idx1, 
						    rb0, rb1));

			/* Print relative statistical error */
			
			fprintf(fp, "%7.5f ", RelErr(ptr, idx0, idx1, 
						     rb0, rb1));
			
			/* Print newline */
			
			fprintf(fp, "\n");
			
			/* Update index */
		  
			n1++;
		      }
		  
		  /* Update index */
		  
		  n0++;		  
		}
	    }

	  fprintf(fp, "];\n\n");

	  /*******************************************************************/
	}
      else
	{
	  /*******************************************************************/

	  /***** Serpent 1 type output (with time bins) **********************/

	  fprintf(fp, "\n");

	  fprintf(fp, "DET%s = [\n", GetText(det0 + DET_PTR_NAME));

	  /* Pointer to adjoint detector */
	  
	  if ((det1 = (long)RDB[det0 + DET_PTR_ADJOINT]) < VALID_PTR)
	    {
	      /* Pointer to statistics */
	  
	      ptr = (long)RDB[det0 + DET_PTR_STAT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Loop over bins */

	      n0 = 0;

	      for (tb0 = 0; tb0 < tbins0; tb0++)
	      for (eb0 = 0; eb0 < ebins0; eb0++)
	      for (ub0 = 0; ub0 < ubins0; ub0++)
	      for (cb0 = 0; cb0 < cbins0; cb0++)
	      for (mb0 = 0; mb0 < mbins0; mb0++)
	      for (lb0 = 0; lb0 < lbins0; lb0++)
	      for (rb0 = 0; rb0 < rbins0; rb0++)
	      for (zb0 = 0; zb0 < zbins0; zb0++)
	      for (yb0 = 0; yb0 < ybins0; yb0++)
	      for (xb0 = 0; xb0 < xbins0; xb0++)
		{
		  /* Print indexes */
		  
		  fprintf(fp, "%5ld ", n0 + 1);
		  fprintf(fp, "%4ld ", tb0 + 1);
		  fprintf(fp, "%4ld ", eb0 + 1);
		  fprintf(fp, "%4ld ", ub0 + 1);
		  fprintf(fp, "%4ld ", cb0 + 1);
		  fprintf(fp, "%4ld ", mb0 + 1);
		  fprintf(fp, "%4ld ", lb0 + 1);
		  fprintf(fp, "%4ld ", rb0 + 1);
		  fprintf(fp, "%4ld ", zb0 + 1);
		  fprintf(fp, "%4ld ", yb0 + 1);
		  fprintf(fp, "%4ld ", xb0 + 1);

		  /* Get index */

		  idx0 = DetIdx(det0, eb0, ub0, cb0, mb0, lb0, rb0, yb0, xb0, 
				tb0);

		  /* Print mean */

		  fprintf(fp, "%12.5E ", Mean(ptr, idx0, rb0));
		  
		  /* Print relative statistical error */
		  
		  fprintf(fp, "%7.5f ", RelErr(ptr, idx0, rb0));
		  
		  /* Print newline */

		  fprintf(fp, "\n");

		  /* Update index */
		  
		  n0++;
		}
	    }
	  else
	    {
	      /* Get number of bins for second detector */
	      
	      ebins1 = (long)RDB[det1 + DET_N_EBINS];
	      ubins1 = (long)RDB[det1 + DET_N_UBINS];
	      cbins1 = (long)RDB[det1 + DET_N_CBINS];
	      mbins1 = (long)RDB[det1 + DET_N_MBINS];
	      lbins1 = (long)RDB[det1 + DET_N_LBINS];
	      rbins1 = (long)RDB[det1 + DET_N_RBINS];
	      tbins1 = (long)RDB[det1 + DET_N_TBINS];

	      /* Mesh bins */
	      
	      if ((ptr = (long)RDB[det1 + DET_PTR_MESH]) > VALID_PTR)
		{
		  xbins1 = (long)RDB[ptr + MESH_N0];
		  ybins1 = (long)RDB[ptr + MESH_N1];
		  zbins1 = (long)RDB[ptr + MESH_N2];
		}
	      else
		{
		  xbins1 = 1;
		  ybins1 = 1;
		  zbins1 = 1;
		}

	      /* Pointer to statistics */
	  
	      ptr = (long)RDB[det0 + DET_PTR_STAT];
	      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	      /* Loop over bins */

	      n0 = 0;

	      for (tb0 = 0; tb0 < tbins0; tb0++)
	      for (eb0 = 0; eb0 < ebins0; eb0++)
	      for (ub0 = 0; ub0 < ubins0; ub0++)
	      for (cb0 = 0; cb0 < cbins0; cb0++)
	      for (mb0 = 0; mb0 < mbins0; mb0++)
	      for (lb0 = 0; lb0 < lbins0; lb0++)
	      for (rb0 = 0; rb0 < rbins0; rb0++)
	      for (zb0 = 0; zb0 < zbins0; zb0++)
	      for (yb0 = 0; yb0 < ybins0; yb0++)
	      for (xb0 = 0; xb0 < xbins0; xb0++)
		{
		  /* Loop over second bins */

		  n1 = 0;

		  for (tb1 = 0; tb1 < ebins1; tb1++)
		  for (eb1 = 0; eb1 < tbins1; eb1++)
		  for (ub1 = 0; ub1 < ubins1; ub1++)
		  for (cb1 = 0; cb1 < cbins1; cb1++)
		  for (mb1 = 0; mb1 < mbins1; mb1++)
		  for (lb1 = 0; lb1 < lbins1; lb1++)
		  for (rb1 = 0; rb1 < rbins1; rb1++)
		  for (zb1 = 0; zb1 < zbins1; zb1++)
		  for (yb1 = 0; yb1 < ybins1; yb1++)
		  for (xb1 = 0; xb1 < xbins1; xb1++)
		      {
			/* Print indexes */
			
			fprintf(fp, "%5ld ", n0 + 1);
			fprintf(fp, "%4ld ", tb0 + 1);
			fprintf(fp, "%4ld ", eb0 + 1);
			fprintf(fp, "%4ld ", ub0 + 1);
			fprintf(fp, "%4ld ", cb0 + 1);
			fprintf(fp, "%4ld ", mb0 + 1);
			fprintf(fp, "%4ld ", lb0 + 1);
			fprintf(fp, "%4ld ", rb0 + 1);
			fprintf(fp, "%4ld ", zb0 + 1);
			fprintf(fp, "%4ld ", yb0 + 1);
			fprintf(fp, "%4ld ", xb0 + 1);

			fprintf(fp, "%5ld ", n1 + 1);
			fprintf(fp, "%4ld ", tb1 + 1);
			fprintf(fp, "%4ld ", eb1 + 1);
			fprintf(fp, "%4ld ", ub1 + 1);
			fprintf(fp, "%4ld ", cb1 + 1);
			fprintf(fp, "%4ld ", mb1 + 1);
			fprintf(fp, "%4ld ", lb1 + 1);
			fprintf(fp, "%4ld ", rb1 + 1);
			fprintf(fp, "%4ld ", zb1 + 1);
			fprintf(fp, "%4ld ", yb1 + 1);
			fprintf(fp, "%4ld ", xb1 + 1);
			
			/* Get indexes */

			idx0 = DetIdx(det0, eb0, ub0, cb0, mb0, lb0, zb0, yb0, 
				      xb0, tb0);

			idx1 = DetIdx(det1, eb1, ub1, cb1, mb1, lb1, zb1, yb1, 
				      xb1, tb1);

			/* Print mean */

			fprintf(fp, "%12.5E ", Mean(ptr, idx0, idx1, 
						    rb0, rb1));
			
			/* Print relative statistical error */
			
			fprintf(fp, "%7.5f ", RelErr(ptr, idx0, idx1,
						     rb0, rb1));
			
			/* Print newline */
			
			fprintf(fp, "\n");
			
			/* Update index */
		  
			n1++;
		      }
		  
		  /* Update index */
		  
		  n0++;		  
		}
	    }

	  fprintf(fp, "];\n\n");

	  /*******************************************************************/
	}
    
      /***********************************************************************/
      
      /***** Print energy intervals ******************************************/

      if ((erg = (long)RDB[det0 + DET_PTR_EGRID]) > VALID_PTR)
	{
	  /* Pointer to values */
	  
	  ptr = (long)RDB[erg + ENERGY_GRID_PTR_DATA];
	  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);

	  /* Loop over energy bins and print values */
	  
	  fprintf(fp, "\nDET%sE = [\n", GetText(det0 + DET_PTR_NAME));
	  
	  for (n0 = 0; n0 < (long)RDB[erg + ENERGY_GRID_NE] - 1; n0++)
	    fprintf(fp, "%12.5E %12.5E %12.5E\n", RDB[ptr + n0],
		    RDB[ptr + n0 + 1], 
		    (RDB[ptr + n0] + RDB[ptr + n0 + 1])/2.0);
	  
	  fprintf(fp, "];\n");
	}

      /***********************************************************************/

      /***** Print time intervals ********************************************/

      if ((ptr = (long)RDB[det0 + DET_PTR_TME]) > VALID_PTR)
	{
	  /* Loop over time bins and print values */
	  
	  fprintf(fp, "\nDET%sT = [\n", GetText(det0 + DET_PTR_NAME));

	  for (n0 = 0; n0 < (long)RDB[det0 + DET_N_TBINS]; n0++)
	    fprintf(fp, "%12.5E %12.5E %12.5E\n", RDB[ptr + n0],
		    RDB[ptr + n0 + 1],
		    (RDB[ptr + n0] + RDB[ptr + n0 + 1])/2.0);
	  
	  fprintf(fp, "];\n");
	}
      
      /***********************************************************************/
      
      /***** Print mesh intervals ********************************************/

      /* Pointer to mesh */

      if ((ptr = (long)RDB[det0 + DET_PTR_MESH]) > VALID_PTR)
	{
	  /* Get sizes */

	  n0 = (long)RDB[ptr + MESH_N0];
	  n1 = (long)RDB[ptr + MESH_N1];
	  n2 = (long)RDB[ptr + MESH_N2];
	  
	  /* Get dimensions */
	  
	  min0 = RDB[ptr + MESH_MIN0];
	  max0 = RDB[ptr + MESH_MAX0];
	  min1 = RDB[ptr + MESH_MIN1];
	  max1 = RDB[ptr + MESH_MAX1];
	  min2 = RDB[ptr + MESH_MIN2];
	  max2 = RDB[ptr + MESH_MAX2];

	  /* Check type */

	  if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_CARTESIAN)
	    {
	      /* x-direction */

	      if (n0 > 1)
		{
		  fprintf(fp, "\nDET%sX = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  for (n = 0; n < n0; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n0)*(max0 - min0) + min0,
			    ((double)n + 1.0)/((double)n0)*(max0 - min0) + min0,
			    ((double)n + 0.5)/((double)n0)*(max0 - min0) + min0);
	  
		  fprintf(fp, "];\n");
		}

	      /* y-direction */

	      if (n1 > 1)
		{
		  fprintf(fp, "\nDET%sY = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  for (n = 0; n < n1; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n1)*(max1 - min1) + min1,
			    ((double)n + 1.0)/((double)n1)*(max1 - min1) + min1,
			    ((double)n + 0.5)/((double)n1)*(max1 - min1) + min1);
	  
		  fprintf(fp, "];\n");
		}

	      /* z-direction */

	      if (n2 > 1)
		{
		  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  for (n = 0; n < n2; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);
	  
		  fprintf(fp, "];\n");
		}
	    }
	  else if (((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXX) ||
		   ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXY))
	    {
	      /* Print cell center coordinates */

	      if ((n0 > 0) && (n1 > 0))
		{
		  /* Get pitch */
	      
		  pitch = RDB[ptr + MESH_MAX0];
		  
		  /* Get center coordinates */
		  
		  x0 = RDB[ptr + MESH_MIN0] + (1 - (n0 % 2))*0.5*pitch;
		  y0 = RDB[ptr + MESH_MIN1] + (1 - (n1 % 2))*0.5*pitch;

		  /* Print */

		  fprintf(fp, "\nDET%sCOORD = [\n", 
			  GetText(det0 + DET_PTR_NAME));

		  /* Avoid compiler warning */

		  x = 0.0;
		  y = 0.0;

		  /* Loop over lattice */

		  i = -(long)((double)n0/2.0);
		  for (n = 0; n < n0; n++)
		    {
		      j = -(long)((double)n1/2.0);
		      for (m = 0; m < n1; m++)
			{
			  if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXX)
			    {
			      x = x0 + (i + COS60*j)*pitch; 
			      y = y0 + j*SIN60*pitch; 
			    }
			  else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_HEXY)
			    {
			      x = x0 + j*SIN60*pitch; 
			      y = y0 + (i + COS60*j)*pitch; 
			    }
			  
			  j++;
			  
			  fprintf(fp, "%E %E\n", x, y);
			}
		      i++;
		    }
		  
		  fprintf(fp, "];\n");
		}

	      /* z-direction */

	      if (n2 > 1)
		{
		  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  for (n = 0; n < n2; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);
	  
		  fprintf(fp, "];\n");
		}
	    }
	  else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_CYLINDRICAL)
	    {
	      /* r-direction */

	      if (n0 > 1)
		{
		  fprintf(fp, "\nDET%sR = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  for (n = 0; n < n0; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n0)*(max0 - min0) + min0,
			    ((double)n + 1.0)/((double)n0)*(max0 - min0) + min0,
			    ((double)n + 0.5)/((double)n0)*(max0 - min0) + min0);
	  
		  fprintf(fp, "];\n");
		}

	      /* phi-direction */

	      if (n1 > 1)
		{
		  fprintf(fp, "\nDET%sPHI = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  min1 = min1*180.0/PI;
		  max1 = max1*180.0/PI;

		  for (n = 0; n < n1; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n1)*(max1 - min1) + min1,
			    ((double)n + 1.0)/((double)n1)*(max1 - min1) + min1,
			    ((double)n + 0.5)/((double)n1)*(max1 - min1) + min1);
	  
		  fprintf(fp, "];\n");
		}

	      /* z-direction */

	      if (n2 > 1)
		{
		  fprintf(fp, "\nDET%sZ = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  for (n = 0; n < n2; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);
	  
		  fprintf(fp, "];\n");
		}
	    }
	  else if ((long)RDB[ptr + MESH_TYPE] == MESH_TYPE_SPHERICAL)
	    {
	      /* r-direction */

	      if (n0 > 1)
		{
		  fprintf(fp, "\nDET%sR = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  for (n = 0; n < n0; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n0)*(max0 - min0) + min0,
			    ((double)n + 1.0)/((double)n0)*(max0 - min0) + min0,
			    ((double)n + 0.5)/((double)n0)*(max0 - min0) + min0);
	  
		  fprintf(fp, "];\n");
		}

	      /* phi-direction */

	      if (n1 > 1)
		{
		  fprintf(fp, "\nDET%sPHI = [\n", GetText(det0 + DET_PTR_NAME));
	  
		  min1 = min1*180.0/PI;
		  max1 = max1*180.0/PI;

		  for (n = 0; n < n1; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n1)*(max1 - min1) + min1,
			    ((double)n + 1.0)/((double)n1)*(max1 - min1) + min1,
			    ((double)n + 0.5)/((double)n1)*(max1 - min1) + min1);
	  
		  fprintf(fp, "];\n");
		}

	      /* theta-direction */

	      if (n2 > 1)
		{
		  fprintf(fp, "\nDET%sTHETA = [\n", GetText(det0 + DET_PTR_NAME));

		  min2 = min2*180.0/PI;
		  max2 = max2*180.0/PI;
	  
		  for (n = 0; n < n2; n++)
		    fprintf(fp, "%12.5E %12.5E %12.5E\n", 
			    ((double)n)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 1.0)/((double)n2)*(max2 - min2) + min2,
			    ((double)n + 0.5)/((double)n2)*(max2 - min2) + min2);
		  
		  fprintf(fp, "];\n");
		}
	    }
	  else
	    Die(FUNCTION_NAME, "Invalid mesh type");	  
	}
  
      /***********************************************************************/

      /* Next detector */
      
      det0 = NextItem(det0);
    }
      
  /* Close file */

  fclose(fp);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
