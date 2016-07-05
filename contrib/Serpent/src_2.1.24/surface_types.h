/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : surface_types.h                                */
/*                                                                           */
/* Created:       2010/10/06 (JLe)                                           */
/* Last modified: 2015/03/05 (JLe)                                           */
/* Version:       2.1.23                                                     */
/*                                                                           */
/* Description: Surface types, names and number of parameters                */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

/* Surface names */

char *surf_types[] = {
  "cyl",
  "px",
  "py",
  "pz",
  "inf",
  "sqc",
  "hexyc",
  "hexxc",
  "sph",
  "cross",
  "pad",
  "cube",
  "cone",
  "svc",
  "cuboid",
  "hexyprism",
  "hexxprism",
  "dode",
  "octa",
  "astra",
  "plane",
  "quadratic",
  "cylx",
  "cyly",
  "cylz",
  "gcross",
  "ppd",
  "rect",
  "usr",
  "\0"
};

/* Number of parameters */
  
int surf_params[SURFACE_TYPES][2] = {
  {3,5},
  {1,1},
  {1,1},
  {1,1},
  {0,0},
  {3,4},
  {3,4},
  {3,4},
  {4,4},
  {4,5},
  {6,6},
  {4,4},
  {5,5},
  {3,3},
  {6,6},
  {5,5},
  {5,5},
  {3,4},
  {3,4},
  {5,5},
  {1,9},
  {1,10},
  {3,5},
  {3,5},
  {3,5},
  {3, MAX_SURFACE_PARAMS},
  {9,9},
  {4,4},
  {0, MAX_SURFACE_PARAMS},
};

/*****************************************************************************/
