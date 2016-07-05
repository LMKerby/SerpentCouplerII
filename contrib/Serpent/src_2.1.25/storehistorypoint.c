/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : storehistorypoint.c                            */
/*                                                                           */
/* Created:       2011/05/13 (JLe)                                           */
/* Last modified: 2015/10/02 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Adds point in particle track in history array                */
/*                                                                           */
/* Comments: - Fixed source -laskussa pitää jotenkin huomioida se että       */
/*             lähdeneutroneilla / protoneilla ei oo historiaa               */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "StoreHistoryPoint:"

/*****************************************************************************/

void StoreHistoryPoint(long part, long mat, long rea, double x, double y, 
		       double z, double u, double v, double w, double E, 
		       double t, double wgt, double flx, long trk) 
{  
#ifdef OLD_HIST
 
  long hst;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Pointer to history data */

  if ((hst = (long)RDB[part + PARTICLE_PTR_HIST]) < VALID_PTR)
    return;

  /* Store data */

  WDB[hst + HIST_X] = x;
  WDB[hst + HIST_Y] = y;
  WDB[hst + HIST_Z] = z;
  WDB[hst + HIST_U] = u;
  WDB[hst + HIST_V] = v;
  WDB[hst + HIST_W] = w;
  WDB[hst + HIST_E] = E;
  WDB[hst + HIST_T] = t;
  WDB[hst + HIST_WGT] = wgt;
  WDB[hst + HIST_PTR_REA] = (double)rea;
  WDB[hst + HIST_PTR_MAT] = (double)mat;
  WDB[hst + HIST_TRK] = (double)trk;
  WDB[hst + HIST_FLX] = (double)flx;
  WDB[hst + HIST_IDX] = RDB[part + PARTICLE_HISTORY_IDX];

  /* Get pointer to next */

  hst = NextItem(hst);
  CheckPointer(FUNCTION_NAME, "(hst)", DATA_ARRAY, hst);

  /* Update pointer */

  WDB[part + PARTICLE_PTR_HIST] = (double)hst;

#ifdef DEBUG

  /* Pointer to previous */

  hst = PrevItem(hst);
  CheckPointer(FUNCTION_NAME, "(hst)", DATA_ARRAY, hst);

  hst = PrevItem(hst);
  CheckPointer(FUNCTION_NAME, "(hst)", DATA_ARRAY, hst);

  /* Check time */

  if ((RDB[hst + HIST_WGT]) > 0.0)
    if (t < RDB[hst + HIST_T])
      Die(FUNCTION_NAME, "Error in time");
    
#endif

#else

  long ptr, flags, ok;

  /* Check if events are recorded */

  if ((flags = (long)RDB[DATA_EVENT_RECORD_FLAGS]) == 0)
    return;

  /* Check pointer */

  CheckPointer(FUNCTION_NAME, "(part)", DATA_ARRAY, part);

  /* Reset OK flag */

  ok = NO;

  /* Types for plotter */

  if ((trk == TRACK_END_STRT) || (trk == TRACK_END_COLL) || 
      (trk == TRACK_END_LEAK) || (trk == TRACK_END_SURF) ||
      (trk == TRACK_END_TCUT) || (trk == TRACK_END_ECUT) ||
      (trk == TRACK_END_WCUT) || (trk == TRACK_END_BC) ||
      (trk == TRACK_END_WWIN))
    if (flags & RECORD_EVENT_PLOTTER)
      ok = YES;

  /* Types for importance */

  if (flx > 0.0) 
    if ((trk == TRACK_END_COLL) || (trk == TRACK_END_VIRT))
      if (flags & RECORD_EVENT_IMPORTANCE)
	ok = YES;

  /* Types for IFP */

  if (trk == EVENT_TYPE_FISS)
    if (flags & RECORD_EVENT_IFP)
      ok = YES;

  /* Check ok flag */

  if (ok == NO)
    return;

  /* New event from bank */
  
  ptr = EventFromBank(part);
  CheckPointer(FUNCTION_NAME, "(ptr)", DATA_ARRAY, ptr);
  
  /* Put type */
  
  WDB[ptr + EVENT_TYPE] = (double)trk;
  
  /* Put coordinates */
  
  CheckValue(FUNCTION_NAME, "x", "", x, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "y", "", y, -INFTY, INFTY);
  CheckValue(FUNCTION_NAME, "z", "", z, -INFTY, INFTY);
  
  WDB[ptr + EVENT_X] = x;
  WDB[ptr + EVENT_Y] = y;
  WDB[ptr + EVENT_Z] = z;
  
  /* Put energy and time */
  
  CheckValue(FUNCTION_NAME, "E", "", E, ZERO, INFTY);
  CheckValue(FUNCTION_NAME, "t", "", t, 0.0, INFTY);
  
  WDB[ptr + EVENT_E] = E;
  WDB[ptr + EVENT_T] = t;
  
  /* Weight and flux */
  
  WDB[ptr + EVENT_WGT] = wgt;
  WDB[ptr + EVENT_FLX] = flx;

  /* Put material */

  WDB[ptr + EVENT_PTR_MAT] = mat;

#endif

  /***************************************************************************/
}

/*****************************************************************************/
