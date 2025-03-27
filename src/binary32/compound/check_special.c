/* Generate special cases for compoundf testing.

Copyright (c) 2022-2025 Stéphane Glondu and Paul Zimmermann, Inria.

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
#include <unistd.h>
#include <fenv.h>
#include <mpfr.h>
#include <math.h>
#include <assert.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

float cr_compoundf (float, float);
float ref_compound (float, float);
int ref_fesetround (int);
void ref_init (void);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

int mid = 1; // if mid=1, also check midpoint cases

typedef union {float f; uint32_t u;} b32u32_u;

static inline uint32_t
asuint (float f)
{
  b32u32_u v = {.f = f};
  return v.u;
}

static inline float
asfloat (uint32_t u)
{
  b32u32_u v = {.u = u};
  return v.f;
}

static float
get_random (void)
{
  b32u32_u v;
  int64_t l;
  l = rand ();
  v.u = l;
  // rand generates only 31 bits
  l = rand ();
  v.u |= (uint32_t) l << 31;
  return v.f;
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
is_nan (float x)
{
  uint32_t u = asuint (x);
  int e = u >> 23;
  return (e == 0xff || e == 0x1ff) && (u << 9) != 0;
}

static inline int
is_qnan (float x)
{
  uint32_t u = asuint (x);
  int e = u >> 23;
  return (e == 0xff || e == 0x1ff) && (u & 0x400000) != 0;
}

static inline int
is_equal (float x, float y)
{
  if (is_nan (x))
    return is_nan (y);
  if (is_nan (y))
    return is_nan (x);
  return asuint (x) == asuint (y);
}

static void
check (float x, float y)
{
  float z1, z2;
  ref_init();
  ref_fesetround(rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  z1 = ref_compound (x, y);
#ifdef CORE_MATH_CHECK_INEXACT
  mpfr_flags_t inex1 = mpfr_flags_test (MPFR_FLAGS_INEXACT);
#endif
  fesetround(rnd1[rnd]);
  feclearexcept (FE_INEXACT);
  z2 = cr_compoundf(x, y);
#ifdef CORE_MATH_CHECK_INEXACT
  int inex2 = fetestexcept (FE_INEXACT);
#endif
  if (!is_equal (z1, z2)) {
    printf("FAIL x,y=%a,%a ref=%a z=%a\n", x, y, z1, z2);
    fflush(stdout);
    exit(1);
  }
#ifdef CORE_MATH_CHECK_INEXACT
  if ((inex1 == 0) && (inex2 != 0))
  {
    printf ("Spurious inexact exception for x=%a y=%a (z=%a)\n", x, y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
  if ((inex1 != 0) && (inex2 == 0))
  {
    printf ("Missing inexact exception for x=%a y=%a (z=%a)\n", x, y, z1);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit(1);
#endif
  }
#endif
}

#ifndef CORE_MATH_TESTS
#define CORE_MATH_TESTS 10000000ul // total number of tests
#endif

static void
check_random (int seed, int nthreads)
{
  ref_init ();
  ref_fesetround (rnd);
  fesetround(rnd1[rnd]);
  float x, y;
  srand (seed);
  for (uint64_t n = 0; n < CORE_MATH_TESTS; n += nthreads)
  {
    x = get_random ();
    y = get_random ();
    check (x, y);
  }
}

// return the largest y such that (1+x)^y >= threshold for -1 < x < 0
static float
underflow_threshold_neg (float x, float threshold)
{
  int rnd = fegetround ();
  fesetround (FE_TOWARDZERO);
  float a, b, c;
  a = 0x1p0f;
  b = 0x1.fffffep+127f;
  while (1) {
    c = (a + b) /2;
    if (c == a || c == b) break;
    float z = ref_compound (x, c);
    if (z < threshold) // (1+x)^y is decreasing with y for -1 < x < 0
      b = c;
    else
      a = c;
  }
  fesetround (rnd);
  return a;
}

// return the smallest y such that (1+x)^y >= threshold for x > 0
static float
underflow_threshold_pos (float x, float threshold)
{
  int rnd = fegetround ();
  fesetround (FE_TOWARDZERO);
  float a, b, c;
  a = -0x1.fffffep+127f;
  b = -1.0f;
  while (1) {
    c = (a + b) /2;
    if (c == a || c == b) break;
    float z = ref_compound (x, c);
    if (z < threshold) // (1+x)^y is increasing with y for x > 0
      a = c;
    else
      b = c;
  }
  fesetround (rnd);
  return b;
}

// return the largest y such that (1+x)^y <= threshold for -1 < x < 0
static float
overflow_threshold_neg (float x, float threshold)
{
  int rnd = fegetround ();
  fesetround (FE_TOWARDZERO);
  float a, b, c;
  a = -0x1.fffffep+127f;
  b = -1.0f;
  while (1) {
    c = (a + b) /2;
    if (c == a || c == b) break;
    float z = ref_compound (x, c);
    if (z <= threshold)
      b = c;
    else
      a = c;
  }
  fesetround (rnd);
  return b;
}

// return the largest y such that (1+x)^y <= threshold for x > 0
static float
overflow_threshold_pos (float x, float threshold)
{
  int rnd = fegetround ();
  fesetround (FE_TOWARDZERO);
  float a, b, c;
  a = 1.0f;
  b = 0x1.fffffep+127f;
  while (1) {
    c = (a + b) /2;
    if (c == a || c == b) break;
    float z = ref_compound (x, c);
    if (z <= threshold)
      a = c;
    else
      b = c;
  }
  fesetround (rnd);
  return a;
}

#define SKIP 1000

/* check corner values: for each x > -1, check the extremal values such that
   (1+x)^y is near the underflow or overflow thresholds */
static void
check_corner (void)
{
  // check near underflow
  uint32_t u0, u1;
  u0 = asuint (-0x1.fffffep-1f) - (rand () % SKIP);
  u1 = asuint (-0x1.9d1da2p-122f);
  assert (u0 > u1);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint32_t u = u0; u >= u1; u -= SKIP)
  {
    float x = asfloat (u);
    float y = underflow_threshold_neg (x, 0x1p-149f);
    check (x, y);
    check (x, nextafterf (y, 0x1.fffffep+127f));
  }
  u0 = asuint (0x1.9d1da2p-122f) + (rand () % SKIP);
  u1 = asuint (0x1.fffffep+127f);
  assert (u0 < u1);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint32_t u = u0; u <= u1; u += SKIP)
  {
    float x = asfloat (u);
    float y = underflow_threshold_pos (x, 0x1p-149f);
    check (x, y);
    check (x, nextafterf (y, -0x1.fffffep+127f));
    x = nextafterf (x, 0x1.fffffep+127f);
  }

  // check near overflow
  u0 = asuint (-0x1.fffffep-1f) - (rand () % SKIP);
  u1 = asuint (-0x1.62e432p-122f);
  assert (u0 > u1);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint32_t u = u0; u >= u1; u -= SKIP)
  {
    float x = asfloat (u);
    float y = overflow_threshold_neg (x, 0x1.fffffep+127f);
    check (x, y);
    check (x, nextafterf (y, -0x1.fffffep+127f));
    x = nextafterf (x, 0);
  }
  u0 = asuint (0x1.62e432p-122f) + (rand () % SKIP);
  u1 = asuint (0x1.fffffep+127f);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for (uint32_t u = u0; u <= u1; u += SKIP)
  {
    float x = asfloat (u);
    float y = overflow_threshold_pos (x, 0x1.fffffep+127f);
    check (x, y);
    check (x, nextafterf (y, 0x1.fffffep+127f));
    x = nextafterf (x, 0x1.fffffep+127f);
  }
}

static void
check_random_all (void)
{
  int nthreads = 1;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel
  nthreads = omp_get_num_threads ();
#endif
#pragma omp parallel for
  for (int i = 0; i < nthreads; i++)
    check_random (getpid () + i, nthreads);
}

static void
check_exact_midpoint (void)
{
  for (float x = 1.0f; x <= 100.0f; x += 1.0f)
    for (float y = -100.0f; y <= 100.0f; y += 1.0f)
      check (x, y);
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

  printf ("Checking exact/midpoint values\n");
  check_exact_midpoint ();

  printf ("Checking near the underflow threshold\n");
  check_corner ();

  printf ("Checking random values\n");
  check_random_all ();

  return 0;
}
