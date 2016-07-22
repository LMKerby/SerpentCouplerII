/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : element.c                                      */
/*                                                                           */
/* Created:       2016/02/12 (JLe)                                           */
/* Last modified: 2016/02/12 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Decomposes natural elements into isotopic compositions       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "natural_elements.h"
#include "element_data.h"

#define FUNCTION_NAME "Element:"

/*****************************************************************************/

void Element(char *name, char *dens, char *id)
{
  double d, tot;
  long Z, A, n;

  /* Get Z */

  if ((Z = atol(name)) < 1)
    for (n = 0; n < NUMBER_OF_ELEMENTS; n++)
      if (!strcasecmp(name, element_symbols[n]))
	{
	  /* Set Z */

	  Z = n;

	  /* Break loop */

	  break;
	}

  /* Check */

  if ((Z < 1) || (Z > NUMBER_OF_ELEMENTS))
    Error(0, "Element %s does not exist", name);

  /* Check if natural composition exists */

  n = 0;
  while(nat_frac[n][0] > 0)
    {
      /* Check Z */
      
      if (nat_frac[n][0] == Z)
	break;
      
      /* Next */
      
      n++;
    }
  
  /* Check */

  if (nat_frac[n][0] == Z)
    fprintf(out, "\nIsotopic composition for natural %s:\n\n", 
	    element_names[Z]);
  else
    {
      /* Not found */
      
      fprintf(out, "\nElement %s has no natural isotopes.\n\n", 
	     element_names[Z]);

      /* Exit */

      return;
    }

  /* Check density */

  if ((d = atof(dens)) < 0.0)
    {
      /* Calculate total for mass fraction */
      
      tot = 0.0;

      n = 0;
      while(nat_frac[n][0] > 0)
	{
	  /* Check Z and add to total */

	  if (nat_frac[n][0] == Z)
	    tot = tot + nat_frac[n][2]*nat_frac[n][3];
	  
	  /* Next */
	  
	  n++;
	}

      /* Loop over data */

      n = 0;
      while(nat_frac[n][0] > 0)
	{
	  /* Check Z and print */

	  if (nat_frac[n][0] == Z)
	    {
	      /* Get A */
	      
	      A = nat_frac[n][1];
	      
	      /* Check if id is given */
	      
	      if (id == NULL)
		fprintf(out, "%6ld  %1.5E\n", 1000*Z + A, 
			nat_frac[n][2]*nat_frac[n][3]*d/tot);
	      else
		fprintf(out, "%6ld.%s  %1.5E\n", 1000*Z + A, id,
			nat_frac[n][2]*nat_frac[n][3]*d/tot); 
	    }
	  
	  /* Next */
	  
	  n++;
	}

    }
  else if (d > 0.0)
    {
      /* Loop over data */

      n = 0;
      while(nat_frac[n][0] > 0)
	{
	  /* Check Z and print */

	  if (nat_frac[n][0] == Z)
	    {
	      /* Get A */
	      
	      A = nat_frac[n][1];
	      
	      /* Check if id is given */
	      
	      if (id == NULL)
		fprintf(out, "%6ld  %1.5E\n", 1000*Z + A, nat_frac[n][3]*d);
	      else
		fprintf(out, "%6ld.%s  %1.5E\n", 1000*Z + A, id, 
			nat_frac[n][3]*d);
	    }
	  
	  /* Next */
	  
	  n++;
	}
    }
  else
    Error(0, "Zero density");

  /* Exit OK */

  fprintf(out, "\n");
}

/*****************************************************************************/
