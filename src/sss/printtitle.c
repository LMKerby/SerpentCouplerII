#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : printtitle.c                                   */
/*                                                                           */
/* Created:       2010/09/10 (JLe)                                           */
/* Last modified: 2016/01/18 (JLe)                                           */
/* Version:       2.1.25                                                     */
/*                                                                           */
/* Description: Prints logo and title                                        */
/*                                                                           */
/* Comments:                                                                 */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "PrintTitle:"

/*****************************************************************************/

void PrintTitle()
{
  char *path;

  fprintf(out, "\n");
  fprintf(out, "  _                   .-=-.           .-=-.          .-==-.       \n");
  fprintf(out, " { }      __        .' O o '.       .' O o '.       /  -<' )--<   \n");
  fprintf(out, " { }    .' O'.     / o .-. O \\     / o .-. O \\     /  .---`       \n");
  fprintf(out, " { }   / .-. o\\   /O  /   \\  o\\   /O  /   \\  o\\   /O /            \n");
  fprintf(out, "  \\ `-` /   \\ O`-'o  /     \\  O`-'o  /     \\  O`-`o /             \n");
  fprintf(out, "   `-.-`     '.____.'       `._____.'       `.____.'              \n");
  
  fprintf(out, "\nSerpent 2 beta\n\n");
  fprintf(out, "A Continuous-energy Monte Carlo Reactor Physics Burnup ");
  fprintf(out, "Calculation Code\n\n");

  fprintf(out, " - Version %s (%s) -- Contact: %s\n\n", CODE_VERSION, CODE_DATE, 
	 CODE_AUTHOR);

  fprintf(out, " - Reference: J. Leppanen, et al. \"The Serpent Monte Carlo code: Status,\n              development and applications in 2013.\" Ann. Nucl. Energy,\n              82 (2015) 142-150.\n\n");

#if defined(__DATE__) && defined(__TIME__)

  fprintf(out, " - Compiled %s %s\n\n", __DATE__, __TIME__);

#endif

#ifdef MPI
  
  fprintf(out, " - MPI Parallel calculation mode available\n\n");

#else

  fprintf(out, " - MPI Parallel calculation mode not available\n\n");

#endif

#ifdef OPEN_MP
  
  fprintf(out, " - OpenMP Parallel calculation mode available\n\n");

#else

  fprintf(out, " - OpenMP Parallel calculation mode not available\n\n");

#endif

#ifdef NOGFX

  fprintf(out, " - Geometry and mesh plotting not available\n\n");

#else
  
  fprintf(out, " - Geometry and mesh plotting available\n\n");

#endif

#ifdef DEBUG

  fprintf(out, " - Source code compiled in debugger mode\n\n");

#endif 

  if ((path = getenv("SERPENT_DATA")) != NULL)
    fprintf(out, " - Default data path set to: \"%s\"\n\n", path);
  else
    fprintf(out, " - Default data path not set\n\n");
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
