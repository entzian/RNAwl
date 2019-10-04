/*
 * globals.h : global definitions for Wang-Landau sampling
 * Last changed Time-stamp: <2014-07-22 15:52:54 mtw>
 */

#ifndef GLOBALS_H
#define GLOBALS_H

#include "config.h"
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>


gsl_histogram *h;     /* histogram of energies seen in current iteration */
gsl_histogram *s;     /* true DOS of the lowest energy range (if
                       * available); required for normalization, which
                       * is computed based on the lowest-energy bins */

/* function pointers */
void          (*pre_process_model)(void);


/* functions */
void
process_commandline(int   argc,
                    char  *argv[]);


void
wanglandau(void);


void
wanglandau_free_memory(void);


void
sighandler(int);


void
dealloc_gengetopt(void);


#endif
