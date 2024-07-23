/* Generate exact cases for cbrt testing.

Copyright (c) 2022 Stéphane Glondu and Paul Zimmermann, Inria.

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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#include <unistd.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int ref_init (void);
int ref_fesetround (int);

double ref_cbrt (double);
double cr_cbrt (double);

int rnd = 0;
int verbose = 0;

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

static void
check (double x)
{
  double z1, z2;
  z1 = ref_cbrt (x);
  fesetround(rnd1[rnd]);
  z2 = cr_cbrt (x);
  if (asuint64 (z1) != asuint64 (z2)) {
    printf("FAIL x=%la ref=%la z=%la\n", x, z1, z2);
    fflush(stdout);
    exit(1);
  }
}

/* check exact cases in binade 2^(i-1) <= x < 2^i */
static void
check_exact (int i)
{
  uint64_t t0, t1;
  if (verbose)
    printf ("Checking 2^%d <= x < 2^%d\n", i-1, i);
  if (i == 0)
  {
    t0 = 208064;
    t1 = 262144;
  }
  if (i == 1)
  {
    t0 = 262144;
    t1 = 330282;
  }
  if (i == 2)
  {
    t0 = 330282;
    t1 = 416128;
  }
  for (uint64_t t = t0; t < t1; t += 2)
  {
    double x = ldexp (t * t * t, -54);
    check (x);
  }
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
      else if (strcmp (argv[1], "--verbose") == 0)
        {
          verbose = 1;
          argc --;
          argv ++;
        }
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  printf ("Checking exact cube roots\n");
  check_exact (0);
  check_exact (1);
  check_exact (2);

  ref_init ();
  ref_fesetround (rnd);

  printf ("Checking random values\n");
#define K 1000000000UL /* total number of tests */
#define BUF_SIZE 1000

  long seed = getpid ();
  srand48 (seed);

  double buf[BUF_SIZE];
  uint64_t N = K / BUF_SIZE;
  for (uint64_t n = 0; n < N; n++)
  {
    /* warning: drand48 is not thread-safe, thus we put it outside
       the parallel loop */
    for (int i = 0; i < BUF_SIZE; i++)
    {
      buf[i] = drand48 ();
      if (i & 1)
        buf[i] = -buf[i];
    }
#pragma omp parallel for
    for (int i = 0; i < BUF_SIZE; i++)
      check (buf[i]);
  }

  return 0;
}
