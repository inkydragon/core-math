/* Generate special cases for exp testing.

Copyright (c) 2022-2024 Stéphane Glondu and Paul Zimmermann, Inria.

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

int ref_fesetround (int);
void ref_init (void);

double cr_exp (double);
double ref_exp (double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

#define MAX_THREADS 192

static unsigned int Seed[MAX_THREADS];

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

static void
check (double x)
{
  double y1 = ref_exp (x);
  fesetround (rnd1[rnd]);
  double y2 = cr_exp (x);
  if (! is_equal (y1, y2))
  {
    printf ("FAIL x=%la ref=%la z=%la\n", x, y1, y2);
    fflush (stdout);
    exit (1);
  }
}

typedef union {double f; uint64_t u;} b64u64_u;

static double
get_random (int tid)
{
  b64u64_u v;
  v.u = rand_r (Seed + tid);
  v.u |= (uint64_t) rand_r (Seed + tid) << 31;
  v.u |= (uint64_t) rand_r (Seed + tid) << 62;
  return v.f;
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

  ref_init();
  ref_fesetround (rnd);

#ifndef CORE_MATH_TESTS
#define CORE_MATH_TESTS 1000000000UL /* total number of tests */
#endif

  /* check subnormal results */
  /* x0 is the smallest x such that 2^-1075 <= exp(x) */
  double x0 = -0x1.74910d52d3051p+9;
  /* x1 is the smallest x such that 2^-1022 <= exp(x) */
  double x1 = -0x1.6232bdd7abcd2p+9;
  /* x2 is the smallest x such that 2^-1074 <= exp(x) */
  double x2 = -0x1.74385446d71c3p+9;
  int64_t n0 = ldexp (x0, 43); /* n0 = -6554261109157969 */
  int64_t n1 = ldexp (x1, 43); /* n1 = -6231120794008786 */
  int64_t n2 = ldexp (x2, 43);
  // first check values around x0
#define K 1000000
  printf ("Checking results around log(2^-1075)\n");
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int64_t n = n0 - K; n < n0 + K; n++)
    check (ldexp ((double) n, -43));
  // then check values around x2
  printf ("Checking results around log(2^-1074)\n");
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int64_t n = n2 - K; n < n2 + K; n++)
    check (ldexp ((double) n, -43));
  // then check values around x1
  printf ("Checking results around log(2^-1022)\n");
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int64_t n = n1 - K; n < n1 + K; n++)
    check (ldexp ((double) n, -43));
  printf ("Checking results in subnormal range\n");
  int64_t skip = (n1 - n0) / CORE_MATH_TESTS + 1;
  n0 += getpid () % skip;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (int64_t n = n0; n < n1; n += skip)
    check (ldexp ((double) n, -43));

  unsigned int seed = getpid ();
  for (int i = 0; i < MAX_THREADS; i++)
    Seed[i] = seed + i;

  printf ("Checking random values\n");
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint64_t n = 0; n < CORE_MATH_TESTS; n++)
  {
    ref_init ();
    ref_fesetround (rnd);
    int tid;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
    tid = omp_get_thread_num ();
#else
    tid = 0;
#endif
    double x = get_random (tid);
    check (x);
  }

  return 0;
}
