#ifdef __cplusplus 
extern "C" { 
#endif 
/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : warn.c                                         */
/*                                                                           */
/* Created:       2010/09/14 (JLe)                                           */
/* Last modified: 2014/08/15 (JLe)                                           */
/* Version:       2.1.22                                                     */
/*                                                                           */
/* Description: Prints note or warning message for user                      */
/*                                                                           */
/* Comments: - Tähän yhdeksi argumentiksi pointteri laskuriin?               */
/*           - Noi messaget vois kerätä myös erilliseen fileen               */
/*           - Printataan stdouttiin?                                        */
/*                                                                           */
/*****************************************************************************/

#include "header.h"

#define FUNCTION_NAME "Note:"

/*****************************************************************************/

void Note(long dummy, ...)
{
  va_list argp;
  va_start (argp, dummy);

  fprintf(out, "\n ***** Warning: ");
  vfprintf(out, va_arg(argp, char *), argp);
  fprintf(out, "\n\n");
}

/*****************************************************************************/
#ifdef __cplusplus 
} 
#endif 
