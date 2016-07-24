#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processtransformations.c                       */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2016/02/20 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Links transformations to universes                           */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ProcessTransformations:"

/*****************************************************************************/

void ProcessTransformations()
{
  long ptr, uni, lvl, surf, cell, next, dummy, bra, trp;
  double det, a11, a12, a13, a21, a22, a23, a31, a32, a33;

  /* Loop over transformations */

  ptr = RDB[DATA_PTR_TR0];
  while (ptr > VALID_PTR)
    {
      /* Get pointer to next */

      next = NextItem(ptr);

      /* Remove transformation from list */

      RemoveItem(ptr);

      /* Check type */

      if ((long)RDB[ptr + TRANS_TYPE] == TRANS_TYPE_UNI)
	{
	  /*******************************************************************/

	  /***** Universe transformation *************************************/

	  /* Find universe */
	  
	  uni = (long)RDB[DATA_PTR_U0];
	  CheckPointer(FUNCTION_NAME, "(uni)", DATA_ARRAY, uni);
	  
	  uni = SeekListStr(uni, UNIVERSE_PTR_NAME, 
			    GetText(ptr + TRANS_PTR_UNI));

	  /* Check */
	  
	  if (uni > VALID_PTR)
	    {
	      /* Set pointer */
	      
	      WDB[ptr + TRANS_PTR_UNI] = (double)uni;

	      /* Add to transformation list */

	      if ((long)RDB[uni + UNIVERSE_PTR_TRANS] > VALID_PTR)
		AddItem(uni + UNIVERSE_PTR_TRANS, ptr);
	      else
		{
		  /* No previous, must create a dummy first to get the list */

		  NewItem(uni + UNIVERSE_PTR_TRANS, TRANS_BLOCK_SIZE);
		  AddItem(uni + UNIVERSE_PTR_TRANS, ptr);

		  /* Remove first */

		  dummy = FirstItem(ptr);
		  RemoveItem(dummy);
		}

	      /* Set used flag */
	      
	      SetOption(ptr + TRANS_OPTIONS, OPT_USED);
	      
	      /* Find level */
	      
	      lvl = (long)RDB[DATA_PTR_LVL0];
	      CheckPointer(FUNCTION_NAME, "(lvl)", DATA_ARRAY, lvl);
	      
	      lvl = SeekList(lvl, LVL_NUMBER, (long)RDB[ptr + TRANS_PTR_LVL],
			     SORT_MODE_ASCEND);
	      
	      /* Check */
	      
	      if (lvl > VALID_PTR)
		{
		  /* Pointer to private data */
		  
		  lvl = (long)RDB[lvl + LVL_PTR_PRIVATE_DATA];
		  CheckPointer(FUNCTION_NAME, "(lvl)", PRIVA_ARRAY, lvl);

		  /* Set pointers */
		  
		  WDB[ptr + TRANS_PTR_LVL] = (double)lvl;
		} 
	      else if ((long)RDB[ptr + TRANS_PTR_LVL] > -1)
		Note(ptr, 
		     "Level %ld in transformation of universe %s does not exist",
		     (long)RDB[ptr + TRANS_PTR_LVL], 
		     GetText(uni + UNIVERSE_PTR_NAME));
	    } 
	  else
	    {
	      /* Check if transformation is used in branch calculation */
	      /* (still removed but avoid printing warning) */

	      bra = (long)RDB[DATA_PTR_BRA0];
	      while (bra > VALID_PTR)
		{
		  /* Loop over transformations */

		  trp = (long)RDB[bra + DEP_BRA_PTR_TRANS];
		  while (trp > VALID_PTR)
		    {
		      /* Compare */

		      if (CompareStr(ptr + TRANS_PTR_UNI,
				     trp + DEP_BRA_TRANS_PTR_TRANS))
			break;
		      
		      /* Next */

		      trp = NextItem(trp);
		    }	

		  if (trp > VALID_PTR)
		    break;

		  /* Next */

		  bra = NextItem(bra);
		}
	      
	      /* Check used-flag */

	      if (bra < VALID_PTR)
		Note(ptr, "Universe %s in transformation does not exist", 
		     GetText(ptr + TRANS_PTR_UNI));
	    }

	  /*******************************************************************/
	}
      else if ((long)RDB[ptr + TRANS_TYPE] == TRANS_TYPE_SURF)
	{
	  /*******************************************************************/

	  /***** Surface transformation **************************************/

	  /* Find surface */
	  
	  surf = (long)RDB[DATA_PTR_S0];
	  CheckPointer(FUNCTION_NAME, "(surf)", DATA_ARRAY, surf);
	  
	  surf = SeekListStr(surf, SURFACE_PTR_NAME, 
			     GetText(ptr + TRANS_PTR_SURF));
	  
	  /* Check pointer */

	  if (surf > VALID_PTR)
	    {
	      /* Add to transformation list */

	      if ((long)RDB[surf + SURFACE_PTR_TRANS] > VALID_PTR)
		AddItem(surf + SURFACE_PTR_TRANS, ptr);
	      else
		{
		  /* No previous, must create a dummy first to get the list */

		  NewItem(surf + SURFACE_PTR_TRANS, TRANS_BLOCK_SIZE);
		  AddItem(surf + SURFACE_PTR_TRANS, ptr);

		  /* Remove first */

		  dummy = FirstItem(ptr);
		  RemoveItem(dummy);
		}
	      
	      /* Set used flag */
	      
	      SetOption(ptr + TRANS_OPTIONS, OPT_USED);
	    }
	  else
	    {
	      /* Check if transformation is used in branch calculation */
	      /* (still removed but avoid printing warning) */

	      bra = (long)RDB[DATA_PTR_BRA0];
	      while (bra > VALID_PTR)
		{
		  /* Loop over transformations */

		  trp = (long)RDB[bra + DEP_BRA_PTR_TRANS];
		  while (trp > VALID_PTR)
		    {
		      /* Compare */

		      if (CompareStr(ptr + TRANS_PTR_SURF,
				     trp + DEP_BRA_TRANS_PTR_TRANS))
			break;
		      
		      /* Next */

		      trp = NextItem(trp);
		    }	

		  if (trp > VALID_PTR)
		    break;

		  /* Next */

		  bra = NextItem(bra);
		}
	      
	      /* Check used-flag */

	      if (bra < VALID_PTR)
		Note(ptr, "Surface %s in transformation does not exist",
		     GetText(ptr + TRANS_PTR_SURF));
	    }

	  /*******************************************************************/
	}
      else if ((long)RDB[ptr + TRANS_TYPE] == TRANS_TYPE_FILL)
	{
	  /*******************************************************************/

	  /***** Fill transformation *****************************************/

	  /* Find cell */
	  
	  cell = (long)RDB[DATA_PTR_C0];
	  CheckPointer(FUNCTION_NAME, "(cell)", DATA_ARRAY, cell);
	  
	  cell = SeekListStr(cell, CELL_PTR_NAME, 
			     GetText(ptr + TRANS_PTR_CELL));
	  
	  /* Check pointer */

	  if (cell > VALID_PTR)
	    {
	      /* Check fill pointer */

	      if ((long)RDB[cell + CELL_PTR_FILL] < VALID_PTR)
		Error(ptr, "Cell %s is not filled by a universe",
		      GetText(cell + CELL_PTR_NAME));

	      /* Add to transformation list */

	      if ((long)RDB[cell + CELL_PTR_TRANS] > VALID_PTR)
		AddItem(cell + CELL_PTR_TRANS, ptr);
	      else
		{
		  /* No previous, must creat a dummy first to get the list */

		  NewItem(cell + CELL_PTR_TRANS, TRANS_BLOCK_SIZE);
		  AddItem(cell + CELL_PTR_TRANS, ptr);

		  /* Remove first */

		  dummy = FirstItem(ptr);
		  RemoveItem(dummy);
		}
	      
	      /* Set used flag */
	      
	      SetOption(ptr + TRANS_OPTIONS, OPT_USED);
	    }
	  else
	    {
	      /* Check if transformation is used in branch calculation */
	      /* (still removed but avoid printing warning) */

	      bra = (long)RDB[DATA_PTR_BRA0];
	      while (bra > VALID_PTR)
		{
		  /* Loop over transformations */

		  trp = (long)RDB[bra + DEP_BRA_PTR_TRANS];
		  while (trp > VALID_PTR)
		    {
		      /* Compare */

		      if (CompareStr(ptr + TRANS_PTR_CELL,
				     trp + DEP_BRA_TRANS_PTR_TRANS))
			break;
		      
		      /* Next */

		      trp = NextItem(trp);
		    }	

		  if (trp > VALID_PTR)
		    break;

		  /* Next */

		  bra = NextItem(bra);
		}
	      
	      /* Check used-flag */

	      if (bra < VALID_PTR)
		Note(ptr, "Cell %s in transformation does not exist",
		     GetText(ptr + TRANS_PTR_CELL));
	    }

	  /*******************************************************************/
	}
      else 
	Die(FUNCTION_NAME, "Invalid transformation type");

      /**********************************************************************/

      /***** Check matrices *************************************************/
      
      if ((long)RDB[ptr + TRANS_ROT] == YES)
	{
	  /* Get coefficients */
	  
	  a11 = RDB[ptr + TRANS_RX1];
	  a12 = RDB[ptr + TRANS_RX2];
	  a13 = RDB[ptr + TRANS_RX3];
	  a21 = RDB[ptr + TRANS_RX4];
	  a22 = RDB[ptr + TRANS_RX5];
	  a23 = RDB[ptr + TRANS_RX6];
	  a31 = RDB[ptr + TRANS_RX7];
	  a32 = RDB[ptr + TRANS_RX8];
	  a33 = RDB[ptr + TRANS_RX9];
	  
	  /* Calculate determinant */
	  
	  det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + 
	    a13*(a21*a32 - a22*a31);
	  
	  /* Check */ 
	  
	  if (fabs(det - 1.0) > 1E-9)
	    Error(ptr, "Transformation matrix must be isometric");

	  /* Renormalize */
	  
	  WDB[ptr + TRANS_RX1] = a11/det;
	  WDB[ptr + TRANS_RX2] = a12/det;
	  WDB[ptr + TRANS_RX3] = a13/det;
	  WDB[ptr + TRANS_RX4] = a21/det;
	  WDB[ptr + TRANS_RX5] = a22/det;
	  WDB[ptr + TRANS_RX6] = a23/det;
	  WDB[ptr + TRANS_RX7] = a31/det;
	  WDB[ptr + TRANS_RX8] = a32/det;
	  WDB[ptr + TRANS_RX9] = a33/det;
	  
	  /* Get coefficients */
	  
	  a11 = RDB[ptr + TRANS_RY1];
	  a12 = RDB[ptr + TRANS_RY2];
	  a13 = RDB[ptr + TRANS_RY3];
	  a21 = RDB[ptr + TRANS_RY4];
	  a22 = RDB[ptr + TRANS_RY5];
	  a23 = RDB[ptr + TRANS_RY6];
	  a31 = RDB[ptr + TRANS_RY7];
	  a32 = RDB[ptr + TRANS_RY8];
	  a33 = RDB[ptr + TRANS_RY9];
	  
	  /* Calculate determinant */
	  
	  det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + 
	    a13*(a21*a32 - a22*a31);
	  
	  /* Check */ 
	  
	  if (fabs(det - 1.0) > 1E-9)
	    Error(ptr, "Transformation matrix must be isometric");
	  
	  /* Renormalize */
	  
	  WDB[ptr + TRANS_RY1] = a11/det;
	  WDB[ptr + TRANS_RY2] = a12/det;
	  WDB[ptr + TRANS_RY3] = a13/det;
	  WDB[ptr + TRANS_RY4] = a21/det;
	  WDB[ptr + TRANS_RY5] = a22/det;
	  WDB[ptr + TRANS_RY6] = a23/det;
	  WDB[ptr + TRANS_RY7] = a31/det;
	  WDB[ptr + TRANS_RY8] = a32/det;
	  WDB[ptr + TRANS_RY9] = a33/det;
	  
	  /* Get coefficients */
	  
	  a11 = RDB[ptr + TRANS_RZ1];
	  a12 = RDB[ptr + TRANS_RZ2];
	  a13 = RDB[ptr + TRANS_RZ3];
	  a21 = RDB[ptr + TRANS_RZ4];
	  a22 = RDB[ptr + TRANS_RZ5];
	  a23 = RDB[ptr + TRANS_RZ6];
	  a31 = RDB[ptr + TRANS_RZ7];
	  a32 = RDB[ptr + TRANS_RZ8];
	  a33 = RDB[ptr + TRANS_RZ9];
	  
	  /* Calculate determinant */
	  
	  det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + 
	    a13*(a21*a32 - a22*a31);
	  
	  /* Check */ 
	  
	  if (fabs(det - 1.0) > 1E-9)
	    Error(ptr, "Transformation matrix must be isometric");
	  
	  /* Renormalize */
	  
	  WDB[ptr + TRANS_RZ1] = a11/det;
	  WDB[ptr + TRANS_RZ2] = a12/det;
	  WDB[ptr + TRANS_RZ3] = a13/det;
	  WDB[ptr + TRANS_RZ4] = a21/det;
	  WDB[ptr + TRANS_RZ5] = a22/det;
	  WDB[ptr + TRANS_RZ6] = a23/det;
	  WDB[ptr + TRANS_RZ7] = a31/det;
	  WDB[ptr + TRANS_RZ8] = a32/det;
	  WDB[ptr + TRANS_RZ9] = a33/det;
	}
      
      /**********************************************************************/

      /* Next transformation */
      
      ptr = next;
    }
}

/****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
