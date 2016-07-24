#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : inithistories.                                 */
/*                                                                           */
/* Created:       2011/04/01 (JLe)                                           */
/* Last modified: 2016/02/01 (VVa)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Initializes particle stacks, ques, source and bank for       */
/*              transport simulation                                         */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "InitHistories:"

/*****************************************************************************/

void InitHistories()
{
  long ptr, loc0, id, np, n;

  fprintf(out, "Allocating memory for particle structures...\n");

  /***************************************************************************/

  /***** Create lists ********************************************************/

  /* Neutron stacks */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_NSTACK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Gamma stacks */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_GSTACK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Precursor stacks */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_PSTACK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Particle ques */

  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_QUE] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Banks */
     
  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_BANK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Stores for dynamic simulation mode */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      /* BOI Store */
     
      loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
      WDB[DATA_PART_PTR_BOI_STORE] = (double)loc0;

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	{
	  ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
	  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
	  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
	}

      /* EOI Store */
     
      loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
      WDB[DATA_PART_PTR_EOI_STORE] = (double)loc0;

      for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
	{
	  ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
	  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
	  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
	}
    }

  /* Separate banks for track plotting */
     
  loc0 = ReallocMem(DATA_ARRAY, (long)RDB[DATA_OMP_MAX_THREADS]);
  WDB[DATA_PART_PTR_TRK_BANK] = (double)loc0;

  for (id = 0; id < (long)RDB[DATA_OMP_MAX_THREADS]; id++)
    {
      ptr = NewItem(loc0++, PARTICLE_BLOCK_SIZE);
      WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
      WDB[ptr + PARTICLE_RNG_IDX] = -1.0;
    }

  /* Source (no division to threads) */
      
  ptr = NewItem(DATA_PART_PTR_SOURCE, PARTICLE_BLOCK_SIZE);
  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;

  /* Precursor source (no division to threads) */
      
  ptr = NewItem(DATA_PART_PTR_PSOURCE, PARTICLE_BLOCK_SIZE);
  WDB[ptr + PARTICLE_TYPE] = (double)PARTICLE_TYPE_DUMMY;
  WDB[ptr + PARTICLE_RNG_IDX] = -1.0;

  /* Set fraction of weight to emit on next interval to 2.0 */
  /* Precursor source will be sorted based on this in a     */
  /* descending order (precursors should have this value    */
  /* in interval (0,1) */

  WDB[ptr + PARTICLE_U] = 2.0;  

  /***************************************************************************/

  /***** Allocate memory for neutrons and photons ****************************/

  if ((long)RDB[DATA_NEUTRON_TRANSPORT_MODE] == YES)
    {
      /* Number of particles */

      if (RDB[DATA_CRIT_POP] > RDB[DATA_SRC_POP])
	np = (long)(RDB[DATA_PART_NBUF_FACTOR]*RDB[DATA_CRIT_POP]);
      else
	np = (long)(RDB[DATA_PART_NBUF_FACTOR]*RDB[DATA_SRC_POP]*
		    (1.0 + 0.25*(RDB[DATA_OMP_MAX_THREADS] - 1.0)));

      /* Allocate memory for neutrons */

      AllocParticleStack(PARTICLE_TYPE_NEUTRON, np);
    }

  if ((long)RDB[DATA_PHOTON_TRANSPORT_MODE] == YES)
    {
      /* Number of particles */
      
      np = (long)(RDB[DATA_PART_GBUF_FACTOR]*RDB[DATA_SRC_POP]*
		  RDB[DATA_OMP_MAX_THREADS]);

      /* Allocate memory for neutrons */

      AllocParticleStack(PARTICLE_TYPE_GAMMA, np);
    }
  

  if ((long)RDB[DATA_PRECURSOR_TRANSPORT_MODE] == PREC_MODE_POINT)
    {
      /* Number of particles */

      if (RDB[DATA_CRIT_POP] > RDB[DATA_SRC_POP])
	np = (long)(RDB[DATA_PART_PBUF_FACTOR]*RDB[DATA_CRIT_POP]);
      else
	np = (long)(RDB[DATA_PART_PBUF_FACTOR]*RDB[DATA_SRC_POP]*
		    (1.0 + 0.25*(RDB[DATA_OMP_MAX_THREADS] - 1.0)));

      np = np*2.0;

      /* Allocate memory for neutrons */

      AllocParticleStack(PARTICLE_TYPE_PRECURSOR, np);
    }
  
  /***************************************************************************/

  /***** Allocate memory for events ******************************************/

  /* Record if any flags are set */

  if ((long)RDB[DATA_EVENT_RECORD_FLAGS] > 0)
    {
      /* Multiply number of events by population size */

      if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_CRIT)
	np = (long)(RDB[DATA_EVENT_BANK_SZ]*RDB[DATA_CRIT_POP]);
      else
	np = (long)RDB[DATA_EVENT_BANK_SZ]*RDB[DATA_SRC_POP];

      /* Allocate memory for events */

      for (n = 0; n < np; n++)
	NewLIFOItem(DATA_PTR_EVENT_BANK, EVENT_BLOCK_SIZE);

      /* Put new value */

      WDB[DATA_EVENT_BANK_SZ] = (double)np;
    }

  /* Allocate memory for particle counts */

  if ((long)RDB[DATA_SIMULATION_MODE] == SIMULATION_MODE_DYN)
    {
      n = (long)RDB[DATA_OMP_MAX_THREADS];
      np = (long)RDB[DATA_SRC_BATCHES];

      ptr = ReallocMem(DATA_ARRAY, n*np*2);
      WDB[DATA_PTR_DYN_PARTCOUNT] = (double)ptr;
    }
  
  /* History index for debugging */
  
  ptr = AllocPrivateData(1, PRIVA_ARRAY);
  WDB[DATA_PTR_PRIVA_HIS_IDX] = (double)ptr;

  /***************************************************************************/

  /* Everything is OK */

  fprintf(out, "OK.\n\n");
  
  /***************************************************************************/
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
