#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : fixhexmesh.c                                   */
/*                                                                           */
/* Created:       2015/03/27 (VVa)                                           */
/* Last modified: 2015/03/27 (VVa)                                           */
/* Version:       2.1.24                                                     */
/*                                                                           */
/* Description: Stores bounding box of tetrahedron based on its points       */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "TetPutBoundingBox:"

/*****************************************************************************/

void TetPutBoundingBox(long ncgns, long pointlist, long pts[4])
{

  /* Reset bounding box for cell */

  WDB[ncgns + IFC_TET_MSH_XMIN] = INFTY;
  WDB[ncgns + IFC_TET_MSH_XMAX] = -INFTY;
  WDB[ncgns + IFC_TET_MSH_YMIN] = INFTY;
  WDB[ncgns + IFC_TET_MSH_YMAX] = -INFTY;
  WDB[ncgns + IFC_TET_MSH_ZMIN] = INFTY;
  WDB[ncgns + IFC_TET_MSH_ZMAX] = -INFTY;

  /* Put bounding box */

  if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[pts[0] + 0])
    WDB[ncgns + IFC_TET_MSH_XMIN]   = RDB[pts[0] + 0];
  if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[pts[0] + 0])
    WDB[ncgns + IFC_TET_MSH_XMAX]   = RDB[pts[0] + 0];      
  if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[pts[0] + 1])
    WDB[ncgns + IFC_TET_MSH_YMIN]   = RDB[pts[0] + 1];
  if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[pts[0] + 1])
    WDB[ncgns + IFC_TET_MSH_YMAX]   = RDB[pts[0] + 1];      
  if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[pts[0] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMIN]   = RDB[pts[0] + 2];
  if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[pts[0] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMAX]   = RDB[pts[0] + 2];      

  if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[pts[1] + 0])
    WDB[ncgns + IFC_TET_MSH_XMIN]   = RDB[pts[1] + 0];
  if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[pts[1] + 0])
    WDB[ncgns + IFC_TET_MSH_XMAX]   = RDB[pts[1] + 0];      
  if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[pts[1] + 1])
    WDB[ncgns + IFC_TET_MSH_YMIN]   = RDB[pts[1] + 1];
  if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[pts[1] + 1])
    WDB[ncgns + IFC_TET_MSH_YMAX]   = RDB[pts[1] + 1];      
  if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[pts[1] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMIN]   = RDB[pts[1] + 2];
  if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[pts[1] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMAX]   = RDB[pts[1] + 2];      

  if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[pts[2] + 0])
    WDB[ncgns + IFC_TET_MSH_XMIN]   = RDB[pts[2] + 0];
  if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[pts[2] + 0])
    WDB[ncgns + IFC_TET_MSH_XMAX]   = RDB[pts[2] + 0];      
  if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[pts[2] + 1])
    WDB[ncgns + IFC_TET_MSH_YMIN]   = RDB[pts[2] + 1];
  if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[pts[2] + 1])
    WDB[ncgns + IFC_TET_MSH_YMAX]   = RDB[pts[2] + 1];      
  if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[pts[2] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMIN]   = RDB[pts[2] + 2];
  if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[pts[2] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMAX]   = RDB[pts[2] + 2];      

  if (RDB[ncgns + IFC_TET_MSH_XMIN] > RDB[pts[3] + 0])
    WDB[ncgns + IFC_TET_MSH_XMIN]   = RDB[pts[3] + 0];
  if (RDB[ncgns + IFC_TET_MSH_XMAX] < RDB[pts[3] + 0])
    WDB[ncgns + IFC_TET_MSH_XMAX]   = RDB[pts[3] + 0];      
  if (RDB[ncgns + IFC_TET_MSH_YMIN] > RDB[pts[3] + 1])
    WDB[ncgns + IFC_TET_MSH_YMIN]   = RDB[pts[3] + 1];
  if (RDB[ncgns + IFC_TET_MSH_YMAX] < RDB[pts[3] + 1])
    WDB[ncgns + IFC_TET_MSH_YMAX]   = RDB[pts[3] + 1];      
  if (RDB[ncgns + IFC_TET_MSH_ZMIN] > RDB[pts[3] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMIN]   = RDB[pts[3] + 2];
  if (RDB[ncgns + IFC_TET_MSH_ZMAX] < RDB[pts[3] + 2])
    WDB[ncgns + IFC_TET_MSH_ZMAX]   = RDB[pts[3] + 2];      

}
#ifdef __cplusplus 
} 
#endif 
