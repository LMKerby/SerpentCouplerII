/* Added to compile on Mac, LMK 6/2016 */

#ifdef __cplusplus
extern "C" {
#endif

#include <header.h>

/* Arrays */

double *ACE;
double *WDB;
double *PRIVA;
double *BUF;
double *RES1;
double *RES2;

const double *RDB;

char *ASCII;

unsigned long *SEED;
unsigned long *SEED0;

/*****************************************************************************/

/* Output pointers */

FILE *err;
FILE *out;

/* Number of mpi tasks and id */

int mpitasks;
int mpiid;

/* Random number seed */

unsigned long parent_seed;

/* Collision counter */

/* Timers */
/* Created tag to compile on Mac, LMK 6/2016 */
struct Timer timer[TOT_TIMERS + 1];


#ifdef __cplusplus
}
#endif
