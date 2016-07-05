/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : processtransformations.c                       */
/*                                                                           */
/* Created:       2010/10/07 (JLe)                                           */
/* Last modified: 2014/08/14 (JLe)                                           */
/* Version:       2.1.22                                                     */
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
  long ptr, uni, lvl, surf, next, dummy;

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
	  /********************************************************************/

	  /***** Universe transformation **************************************/

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
		  /* No previous, must creat a dummy first to get the list */

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
	    Note(ptr, "Universe %s in transformation does not exist", 
		 GetText(ptr + TRANS_PTR_UNI));

	  /********************************************************************/
	}
      else if ((long)RDB[ptr + TRANS_TYPE] == TRANS_TYPE_SURF)
	{
	  /********************************************************************/

	  /***** Surface transformation ***************************************/

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
		  /* No previous, must creat a dummy first to get the list */

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
	    Note(ptr, "Surface %s in transformation does not exist",
		 GetText(ptr + TRANS_PTR_SURF));

	  /********************************************************************/
	}
      else 
	Die(FUNCTION_NAME, "Invalid transformation type");

      /* Next transformation */
      
      ptr = next;
    }
}

/*****************************************************************************/
