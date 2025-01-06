/* Generate special cases for tgamma testing.

Copyright (c) 2022-2024 Paul Zimmermann, Inria.

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
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

int ref_init (void);
int ref_fesetround (int);

double cr_tgamma (double);
double ref_tgamma (double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

#define MAX_THREADS 192

static unsigned int Seed[MAX_THREADS];

typedef union {double f; uint64_t u;} b64u64_u;

static inline uint64_t
asuint64 (double f)
{
  b64u64_u u = {.f = f};
  return u.u;
}

static inline double
asfloat64 (uint64_t n)
{
  b64u64_u u = {.u = n};
  return u.f;
}

static double
get_random (int tid)
{
  b64u64_u v;
  v.u = rand_r (Seed + tid);
  v.u |= (uint64_t) rand_r (Seed + tid) << 31;
  v.u |= (uint64_t) rand_r (Seed + tid) << 62;
  return v.f;
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
  ref_init ();
  ref_fesetround (rnd);
  double y1 = ref_tgamma (x);
  fesetround (rnd1[rnd]);
  double y2 = cr_tgamma (x);
  if (!is_equal (y1,y2))
  {
    printf ("FAIL x=%la ref=%la z=%la\n", x, y1, y2);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
}

typedef union { double f; uint64_t i; } d64u64;

static void
check_negative (void)
{
  for (int n = -1000000; n < 0; n++)
  {
    check (nextafter ((double) n, 0.0));
    check (nextafter ((double) (n+1), (double) n));
    check ((double) n + 0.5);
  }
}

// return x0 such that |gamma(x0)| is minimal on (n,n+1), for n < 0
static double
find_min (int n)
{
  double x0, x1, x2, x3, v0, v1, v2, v3;
  x0 = nextafter (n, n+1);
  x3 = nextafter (n+1, n);
  // use trichotomy
  while (1)
  {
    x1 = (2.0 * x0 + x3) / 3.0;
    x2 = (x0 + 2.0 * x3) / 3.0;
    if (x0 == x1 || x1 == x2 || x2 == x3)
      break;
    v1 = fabs (cr_tgamma (x1));
    v2 = fabs (cr_tgamma (x2));
    if (v1 < v2)
      x3 = x2;
    else
      x0 = x1;
  }
  v0 = fabs (cr_tgamma (x0));
  v1 = fabs (cr_tgamma (x1));
  v2 = fabs (cr_tgamma (x2));
  v3 = fabs (cr_tgamma (x3));
  if (v0 <= v1 && v0 <= v2 && v0 <= v3)
    return x0;
  if (v1 <= v2 && v1 <= v3)
    return x1;
  if (v2 <= v3)
    return x2;
  return x3;
}

#ifndef CORE_MATH_TESTS
#define CORE_MATH_TESTS 1000000000UL /* total number of tests */
#endif

static void
check_subnormal_aux (double x1, double x2)
{
  if (!(x1 <= x2))
    return;
  uint64_t n1 = asuint64 (x1);
  uint64_t n2 = asuint64 (x2);
  // for negative numbers, x1 < x2 means n2 < n1
  int64_t d = n1 - n2;
  // with s=d/40000, we perform 929699 tests
  uint64_t s = d / (CORE_MATH_TESTS / 23);
  if (s == 0)
    s = 1;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint64_t n = n2; n <= n1; n += s)
    check (asfloat64 (n));
}

static void
check_subnormal (void)
{
  double x0, y0, a, b, c, x1, x2, x3, x4;
  for (int k = -184; k <= -171; k++)
  {
    x0 = find_min (k);
    // |gamma(x)| is decreasing on (k,x0), increasing on (x0,k+1)
    y0 = fabs (cr_tgamma (x0));
    if (y0 >= 0x1p-1022)
      continue; // no subnormal value
    a = k;
    b = x0;
    while (nextafter (a, x0) != b)
    {
      c = (a + b) / 2.0;
      if (fabs (cr_tgamma (c)) >= 0x1p-1022)
        a = c;
      else
        b = c;
    }
    x1 = b; // smallest value in (k,x0) such that |gamma(x)| < 2^-1022
    a = k;
    b = x0;
    while (nextafter (a, x0) != b)
    {
      c = (a + b) / 2.0;
      if (fabs (cr_tgamma (c)) >= 0x1p-1074)
        a = c;
      else
        b = c;
    }
    x2 = a; // largest value in (k,x0) such that |gamma(x)| >= 2^-1074
    check_subnormal_aux (x1, x2);
    a = x0;
    b = k+1;
    while (nextafter (a, k+1) != b)
    {
      c = (a + b) / 2.0;
      if (fabs (cr_tgamma (c)) >= 0x1p-1022)
        b = c;
      else
        a = c;
    }
    x4 = a; // smallest value in (x0,k+1) such that |gamma(x)| < 2^-1022
    a = x0;
    b = k+1;
    while (nextafter (a, k+1) != b)
    {
      c = (a + b) / 2.0;
      if (fabs (cr_tgamma (c)) >= 0x1p-1074)
        b = c;
      else
        a = c;
    }
    x3 = b; // largest value in (k,x0) such that |gamma(x)| >= 2^-1074
    check_subnormal_aux (x3, x4);
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
  ref_init ();
  ref_fesetround (rnd);

  printf ("Check subnormal output\n");
  check_subnormal ();

  printf ("Check negative inputs\n");
  check_negative ();

  long seed = getpid ();
  for (int i = 0; i < MAX_THREADS; i++)
    Seed[i] = seed + i;
  
  printf ("Checking random numbers...\n");
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
