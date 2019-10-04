/*
 * wl_rna.c RNA-related routines for Wang-Landau sampling
 * Last changed Time-stamp: <2014-07-22 15:42:10 mtw>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include "config.h"
#include "wl_options.h"
#include "globals.h"
#include "wl_rna.h"

static void
subopt_of_lowest_bins_RNA(float);


/* ==== */
void
pre_process_RNA(void)
{
  subopt_of_lowest_bins_RNA(wanglandau_opt.erange);
}


/* ==== */
static void
subopt_of_lowest_bins_RNA(float e)
{
  int       strucs, have_lowest_bin, status;
  size_t    i;
  SOLUTION  *sol = NULL;

  have_lowest_bin = 0;
  if (wanglandau_opt.verbose) {
    fprintf(stderr, "[[subopt_of_lowest_bins_RNA()]]\nq");
    fprintf(stderr, "computing subopt -e %g\n", e);
  }

  /* compute suboptimal structures within energy range mfe+e */
  sol = subopt(wanglandau_opt.sequence, NULL, e * 100, NULL);
  for (strucs = 0; sol[strucs].structure != NULL; strucs++) {
    status = gsl_histogram_find(s, sol[strucs].energy, &i);
    if (status) {
      if (status == GSL_EDOM)
        printf("error: %s\n", gsl_strerror(status));
      else
        fprintf(stderr, "GSL error: gsl_errno=%d\n", status);

      exit(EXIT_FAILURE);
    }

    if (i == 0)
      have_lowest_bin = 1;

    if (wanglandau_opt.verbose)
      printf("%s %6.2f %ld\n", sol[strucs].structure, sol[strucs].energy, i);

    gsl_histogram_increment(s, sol[strucs].energy);
    free(sol[strucs].structure);
  }
  free(sol);
  if (have_lowest_bin != 1) {
    fprintf(stderr,
            "Lowest bin has not been populated in subopt_first_bin()\n");
    fprintf(stderr,
            "Please decrease the number of bins for this simulation\n");
    fprintf(stderr, "exiting ...\n");
    exit(EXIT_FAILURE);
  }

  /* be verbose about the lower energy bins */
  if (wanglandau_opt.verbose) {
    fprintf(stderr,
            "histogram s (first bin required for normalization)\n");
    for (i = 0; i < wanglandau_opt.truedosbins; i++) {
      double value = gsl_histogram_get(s, i);
      fprintf(stderr, "s[%ld]: %7g\n", i, value);
    }
  }
}
