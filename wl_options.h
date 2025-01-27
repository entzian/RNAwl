/* Last changed Time-stamp: <2014-07-22 16:04:13 mtw> */

#ifndef WL_OPTIONS_H
#define WL_OPTIONS_H

#include <stdio.h>

typedef struct _options {
  FILE      *INFILE;            /* input file */
  char      *basename;          /* base name of processed file */
  char      *structures_file;    /* output file for sampled structures */
  long int  checksteps;         /* wl steps before histogram is checked for
                                 * flatness */
  char      *sequence;          /* sequence */
  char      *structure;         /* start structure */
  int       len;                /* sequence length */
  int       bins;               /* # of equidistant bins in histogram */
  double    ffinal;             /* modification parameter f */
  float     flat;               /* flatness criterion */
  long int  seed;               /* seed for random number generator */
  int       seed_given;         /* whether seed was given at the command line */
  long int  steplimit;          /* maximum number of MC steps to perform */
  float     T;                  /* fold temperature */
  float     erange;             /* energy range for subopt */
  int       norm;               /* # of normalization-bins */
  double    max;                /* upper energy bound of sampling range */
  int       max_given;          /* whether max was given at the command line */
  double    res;                /* histogram bin width */
  int       res_given;          /* whether res was given at the command line */
  int       truedosbins;        /* # of bins that get overwritten by true DOS */
  int       truedosbins_given;  /* whether truedosbins was given */
  int       verbose;            /* be verbose */
  int       debug;              /* debug mode */
} options;

options wanglandau_opt;

#endif
