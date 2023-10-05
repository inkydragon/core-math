/* Check correctness of univariate binary64 function on worst cases.

Copyright (c) 2022-2023 Stéphane Glondu and Paul Zimmermann, Inria.

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#define _POSIX_C_SOURCE 200809L  /* for getline */

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#ifndef CORE_MATH_NO_OPENMP
#include <omp.h>
#endif

#include "function_under_test.h"

double cr_function_under_test (double);
double ref_function_under_test (double);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;

static void
readstdin(double **result, int *count)
{
  char *buf = NULL;
  size_t buflength = 0;
  ssize_t n;
  int allocated = 512;

  *count = 0;
  if (NULL == (*result = malloc(allocated * sizeof(double)))) {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }

  while ((n = getline(&buf, &buflength, stdin)) >= 0) {
    if (n > 0 && buf[0] == '#') continue;
    if (*count >= allocated) {
      int newsize = 2 * allocated;
      double *newresult = realloc(*result, newsize * sizeof(double));
      if (NULL == newresult) {
        fprintf(stderr, "realloc(%d) failed\n", newsize);
        exit(1);
      }
      allocated = newsize;
      *result = newresult;
    }
    double *item = *result + *count;
    if (sscanf(buf, "%la", item) == 1) {
      (*count)++;
    }
  }
}

static inline uint64_t
asuint64 (double f)
{
  union
  {
    double f;
    uint64_t i;
  } u = {f};
  return u.i;
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
is_nan (double x)
{
  uint64_t u = asuint64 (x);
  int e = u >> 52;
  return (e == 0x7ff || e == 0xfff) && (u << 12) != 0;
}

static inline int
is_equal (double x, double y)
{
  if (is_nan (x))
    return is_nan (y);
  if (is_nan (y))
    return is_nan (x);
  return asuint64 (x) == asuint64 (y);
}

void
doloop(void)
{
  double *items;
  int count, tests = 0, failures = 0, skipped = 0;

  readstdin(&items, &count);

#ifndef CORE_MATH_NO_OPENMP
#pragma omp parallel for reduction(+: failures,tests,skipped)
#endif
  for (int i = 0; i < count; i++) {
    ref_init();
    ref_fesetround(rnd);
    fesetround(rnd1[rnd]);
    double x = items[i];
    double z1 = ref_function_under_test(x);
    double z2 = cr_function_under_test(x);
    tests ++;
    /* Note: the test z1 != z2 would not distinguish +0 and -0. */
    if (z2 == 0) skipped++;
    if (z2 != 0 && is_equal (z1, z2) == 0) {
      printf("FAIL x=%la ref=%la z=%la\n", x, z1, z2);
      fflush(stdout);
#ifdef DO_NOT_ABORT
      failures ++;
#else
      exit(1);
#endif
    }
#ifdef WORST_SYMMETRIC
    x = -x;
    z1 = ref_function_under_test(x);
    z2 = cr_function_under_test(x);
    tests ++;
    if (z2 == 0) skipped++;
    if (z2 != 0 && is_equal (z1, z2) == 0) {
      printf("FAIL x=%la ref=%la z=%la\n", x, z1, z2);
      fflush(stdout);
#ifdef DO_NOT_ABORT
      failures ++;
#else
      exit(1);
#endif
    }
#endif /* WORST_SYMMETRIC */
  }

  free(items);
  printf("%d tests passed, %d failure(s), %d skipped\n", tests, failures, skipped);
}

int
main (int argc, char *argv[])
{
  while (argc >= 2)
    {
      if (strcmp (argv[1], "--rndn") == 0)
        {
          rnd = 0;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "--rndz") == 0)
        {
          rnd = 1;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "--rndu") == 0)
        {
          rnd = 2;
          argc --;
          argv ++;
        }
      else if (strcmp (argv[1], "--rndd") == 0)
        {
          rnd = 3;
          argc --;
          argv ++;
        }
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  doloop();
}
