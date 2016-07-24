#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : readofbatches.c                                */
/*                                                                           */
/* Created:       2014/03/30 (JLe)                                           */
/* Last modified: 2014/03/30 (JLe)                                           */
/* Version:       2.1.19                                                     */
/*                                                                           */
/* Description: Reads batch data from OpenFOAM format files                  */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "ReadOFBatches:"

/*****************************************************************************/

void ReadOFBatches(long loc0)
{
  long n, loc1;
  char tmpstr[MAX_STR], c;
  FILE *fp;
  
  /* Loop over files */

  for (n = 0; n < 6; n++)
    {
      /* Get file name */

      if (n == 0) 
	sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_PFILE));
      else if (n == 1) 
	sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_FFILE));
      else if (n == 2) 
	sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_OFILE));
      else if (n == 3) 
	sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_NFILE));
      else if (n == 4) 
	sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_RFILE));
      else if (n == 5) 
	sprintf(tmpstr, "%s", GetText(loc0 +  IFC_PTR_OF_TFILE));
      else
	Die(FUNCTION_NAME, "Overflow");

      /* Open file for reading */

      if ((fp = fopen(tmpstr, "r")) == NULL)
	continue;

      /* Look over file and look for key word "boundaryField" */

      while (fscanf(fp, "%s", tmpstr) != EOF)
	{
	  /* Check string */
	  
	  if (!strcmp(tmpstr, "boundaryField"))
	    {
	      /* Loop until begin marker */
	      
	      while ((c = fgetc(fp)) != EOF)
		if (c == '{')
		  break;
	      
	      /* Loop */
	      
	      while (1 != 2)
		{
		  /* Read next word */
		  
		  if (fscanf(fp, "%s", tmpstr) == EOF)
		    break;
		  
		  /* Check end marker */
		  
		  if (tmpstr[0] == '}')
		    break;
		  else 
		    {
		      /* Find match */

		      loc1 = (long)RDB[loc0 + IFC_PTR_OF_BATCHES];
		      while (loc1 > VALID_PTR)
			{
			  /* Compare names */

			  if (!strcmp(GetText(loc1 + IFC_OF_BATCH_PTR_NAME),
				      tmpstr))
			    break;

			  /* Next */

			  loc1 = NextItem(loc1);
			}
		      
		      /* Check pointer */
		      
		      if (loc1 < VALID_PTR)
			{
			  /* Add new batch */

			  loc1 = NewItem(loc0 + IFC_PTR_OF_BATCHES,
					 IFC_OF_BATCH_BLOCK_SIZE);
			  
			  /* Put name */

			  WDB[loc1 + IFC_OF_BATCH_PTR_NAME] = PutText(tmpstr);
			}
		      
		      /* Loop until end marker */
		      
		      while ((c = fgetc(fp)) != EOF)
			if (c == '}')
			  break;
		    }
		}
	    }
	}
    }
}

/*****************************************************************************/



#ifdef __cplusplus 
} 
#endif 
