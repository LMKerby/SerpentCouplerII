#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getlatticeindexes.c                            */
/*                                                                           */
/* Created:       2010/10/08 (JLe)                                           */
/* Last modified: 2011/10/28 (JLe)                                           */
/* Version:       2.1.0                                                      */
/*                                                                           */
/* Description: Returns indexes in squre and hexagonal lattices              */
/*                                                                           */
/* Comments: - Taken from Serpent 1.1.0                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetLatticeIndexes:"

/*****************************************************************************/

long GetLatticeIndexes(double px, double py, double pz, double x0, double y0, 
		       double z0, long *i, long *j, long *k, long type)
{
  double x, y, n1, n2, n3;
  double mid1, mid2, mid3;

  if ((type == LAT_TYPE_S) || (type == LAT_TYPE_INFS))
    {
      /* Neliöhila, triviaalitapaus */
      
      *i = (long)rint(x0/px);
      *j = (long)rint(y0/py);
      *k = (long)rint(z0/pz);
      
      return 0;
    }

  else
    {
      /* Kolmiohila, käytetään R. Mattilan johtamaa kaavaa */

      if ((type == LAT_TYPE_HY) || (type == LAT_TYPE_INFHY))
	{
	  /* Y-tyyppi -> vaihdetaan koordinaatit */
	  
	  y = x0/px;
	  x = y0/py;
	}
      else
	{
	  /* X-tyyppi -> käytetään suoraan */

	  x = x0/px;
	  y = y0/py;
	}
      
      n1 =        2.0*x - 0.5;
      n2 = -x + SQRT3*y - 0.5;
      n3 = -x - SQRT3*y - 0.5;
      
      mid1 = rint(n1);
      mid2 = rint(n2);
      mid3 = rint(n3);

      *i = (long)floor(0.5 + (mid1 - mid2)/3.0);
      *j = (long)floor(0.5 + (mid2 - mid3)/3.0);

      /* Z-indeksi */

      *k = (long)rint(z0/pz);

      return 0;
    }
  
  fprintf(err, "%s Invalid lattice type %ld.\n", FUNCTION_NAME, type);

  exit(-1);
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
