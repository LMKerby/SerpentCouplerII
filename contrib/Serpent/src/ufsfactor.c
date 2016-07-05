/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : ufsfactor.c                                    */
/*                                                                           */
/* Created:       2012/02/11 (JLe)                                           */
/* Last modified: 2014/05/25 (JLe)                                           */
/* Version:       2.1.21                                                     */
/*                                                                           */
/* Description: Returns multiplier for uniform fission source method         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "UFSFactor:"

/*****************************************************************************/

double UFSFactor(double x, double y, double z, long id)
{
  long msh, ptr, idx, n0, n1, n2, lat, nx, ny, nz, type, i, j, k;
  double f, zmin, zmax, pr, pz;
  
  /* Check if method is in use */

  if ((long)RDB[DATA_UFS_MODE] == UFS_MODE_NONE)
    return 1.0;

  /* Check simulation mode */

  if ((long)RDB[DATA_SIMULATION_MODE] != SIMULATION_MODE_CRIT)
    return 1.0;

  /* Get pointer to mesh */

  msh = (long)RDB[DATA_UFS_PTR_SRC_MESH];
  CheckPointer(FUNCTION_NAME, "(msh)", DATA_ARRAY, msh);

  /* Get mesh size */
  
  n0 = (long)RDB[msh + MESH_N0];
  n1 = (long)RDB[msh + MESH_N1];
  n2 = (long)RDB[msh + MESH_N2];
  
  /* Check values */
  
  CheckValue(FUNCTION_NAME, "n0", "", n0, 1, 10000);
  CheckValue(FUNCTION_NAME, "n1", "", n1, 1, 10000);
  CheckValue(FUNCTION_NAME, "n2", "", n2, 1, 10000);
      
  /* Check if all data is collected */

  if (RDB[DATA_CYCLE_IDX] < RDB[DATA_CRIT_SKIP])
    return 1.0;

  /* Check lattice pointer */
      
  if ((lat = (long)RDB[DATA_UFS_PTR_LAT]) > VALID_PTR)
    {
      /* Get lattice type */
      
      type = (long)RDB[lat + LAT_TYPE];

      /* Check */

      if ((type != LAT_TYPE_S) && (type != LAT_TYPE_HX) && 
	  (type != LAT_TYPE_HY))
	Die(FUNCTION_NAME, "Lattice type %ld not supported with UFS", type);

      /* Get lattice size */
      
      nx = (long)RDB[lat + LAT_NX];
      CheckValue(FUNCTION_NAME, "nx", "", nx, 1, 10000);

      ny = (long)RDB[lat + LAT_NY];
      CheckValue(FUNCTION_NAME, "ny", "", ny, 1, 10000);

      /* Get axial dimension */
      
      nz = (long)RDB[DATA_UFS_NZ];
      CheckValue(FUNCTION_NAME, "nz", "", nz, 1, 10000);

      zmin = RDB[DATA_UFS_ZMIN];
      CheckValue(FUNCTION_NAME, "zmin", "", zmin, -INFTY, INFTY);

      zmax = RDB[DATA_UFS_ZMAX];
      CheckValue(FUNCTION_NAME, "zmax", "", zmax, zmin, INFTY);

      /* Get radial and axial pitch */

      pr = RDB[lat + LAT_PITCH];
      pz = (zmax - zmin)/((double)nz);

      CheckValue(FUNCTION_NAME, "pr", "", pr, ZERO, INFTY);
      CheckValue(FUNCTION_NAME, "pz", "", pz, ZERO, INFTY);

      /* Transfer co-ordinates */
      
      x = x - RDB[lat + LAT_ORIG_X0];
      y = y - RDB[lat + LAT_ORIG_Y0];

      /* Get indexes */
      
      GetLatticeIndexes(pr, pr, pz, x, y, z, &i, &j, &k, type);

      /* Adjust */

      i = i + (long)((double)nx/2.0);
      j = j + (long)((double)ny/2.0);

      /* Check */

      if ((i < 0) || (i > nx - 1) || (j < 0) || (j > ny - 1) ||
	  (k < 0) || (k > nz - 1))
	return 1.0;

      /* Calculate index */

      idx = i + j*n0 + k*n0*n1;
    }
  else
    {
      /* Get index from coordinates */

      idx = MeshIndex(msh, x, y, z);
    }

  /* Check */

  if (idx > -1)      
    {
      /* Get pointer to direct data */
      
      ptr = (long)RDB[DATA_UFS_PTR_FACTORS];
      CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr + idx);
      
      /* Read value */
      
      f = RDB[ptr + idx];
    }
  else
    f = 1.0;

  /* Check factor */

  CheckValue(FUNCTION_NAME, "f", "", f, RDB[DATA_UFS_MIN], RDB[DATA_UFS_MAX]);

  /* Return value */

  return f;
}

/*****************************************************************************/
