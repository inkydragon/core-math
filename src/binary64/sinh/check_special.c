/* Check sinh on random inputs.

Copyright (c) 2022-2025 Paul Zimmermann, Inria.

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
#include <inttypes.h>
#include <assert.h>
#include <mpfr.h>
#include "function_under_test.h"
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

int ref_init (void);
int ref_fesetround (int);

double cr_sinh (double);
double ref_sinh (double);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;
unsigned long tested = 0;

#define MAX_THREADS 192

static unsigned int Seed[MAX_THREADS];

static inline double tfun(double x){
  return cr_function_under_test(x);
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

static void
check (double x)
{
  int bug;
#pragma omp atomic update
  tested ++;
  double y1 = ref_sinh (x);
  fesetround (rnd1[rnd]);
  double y2 = cr_sinh (x);
  if (isnan (y1))
    bug = !isnan (y2);
  else if (isnan (y2))
    bug = !isnan (y1);
  else
    bug = asuint64 (y1) != asuint64 (y2);
  if (bug)
  {
    printf ("FAIL x=%la ref=%la z=%la\n", x, y1, y2);
    fflush (stdout);
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
}

static inline double
asfloat64 (uint64_t i)
{
  b64u64_u u = {.u = i};
  return u.f;
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
is_nan (double x)
{
  uint64_t u = asuint64 (x);
  uint64_t e = u >> 52;
  return (e == 0x7ff || e == 0xfff) && (u << 12) != 0;
}

// When x is a NaN, returns 1 if x is an sNaN and 0 if it is a qNaN
static inline int is_signaling(double x) {
  b64u64_u _x = {.f = x};

  return !(_x.u & (1ull << 51));
}

static inline int
is_inf (double x)
{
  uint64_t u = asuint64 (x);
  uint64_t e = u >> 52;
  return (e == 0x7ff || e == 0xfff) && (u << 12) == 0;
}

static void
check_invalid (void)
{
  double plusInf = asfloat64(0x7ff0000000000000ull);
  double minInf = asfloat64(0xfff0000000000000ull);
  double qnan = asfloat64(0x7ff8000000000000ull);

  // Check +Inf
  feclearexcept (FE_INVALID);
  double y = cr_sinh(plusInf);

  /* In accordance with IEEE-754 sinh(+Inf) is +Inf and sinh(-Inf) is -Inf  */
  if (is_nan (y))
  {
    fprintf (stderr, "Error, foo(+Inf) should be +Inf, got %la=%"PRIx64"\n",
             y, asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // check the invalid exception was not set
  int flag = fetestexcept (FE_INVALID);
  if (flag)
  {
    printf ("Spurious invalid exception for x=+Inf\n");
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // Check that y is Inf
  if (!is_inf(y) || y < 0)
  {
    fprintf (stderr, "Error, foo(+Inf) should be +Inf, got %la=%"PRIx64"\n",
      y, asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }

  // Check -Inf
  feclearexcept (FE_INVALID);
  y = cr_sinh(minInf);
  if (is_nan (y))
  {
    fprintf (stderr, "Error, foo(-Inf) should be -Inf, got %la=%"PRIx64"\n",
             y, asuint64 (y));
#ifndef DO_NOT_ABORT
      exit (1);
#endif
  }
  // check the invalid exception was set
  flag = fetestexcept (FE_INVALID);
  if (flag)
  {
    printf ("Spurious invalid exception for x=-Inf\n");
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }
  // Check that y is Inf
  if (!is_inf(y) || y > 0)
  {
    fprintf (stderr, "Error, foo(-Inf) should be -Inf, got %la=%"PRIx64"\n",
      y, asuint64 (y));
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }

  // Check qnan
  feclearexcept (FE_INVALID);
  y = cr_sinh(qnan);
  // Check that sinh(qNaN) = qNaN
  if (!is_nan (y) || is_signaling(y))
  {
  fprintf (stderr, "Error, foo(qNaN) should be qNaN, got %la=%"PRIx64"\n",
               y, asuint64 (y));
#ifndef DO_NOT_ABORT
        exit (1);
#endif
  }
  // check the invalid exception was set
  flag = fetestexcept (FE_INVALID);
  if (flag)
  {
    printf ("Spurious invalid exception for x=-Inf\n");
#ifndef DO_NOT_ABORT
    exit (1);
#endif
  }

  // Check 
  double T[] = {asfloat64(0x7ff0000000000000ull), // +Inf
    0x1.fffffffffffffp+1023, // DBL_MAX
    0x1.633ce8fb9f87cp+9,
    0x1.633ce8fb9f87dp+9,
    0x1.633ce8fb9f87ep+9,
    0x1.633ce8fb9f87fp+9};

  for (unsigned long i = 0; i < sizeof(T)/sizeof(T[0]); i++) {
    feclearexcept (FE_INVALID);
    double x = T[i];
     y = cr_sinh (x);
    if (x >= 0x1.633ce8fb9f87ep+9 && !is_inf (y))
    {
      fprintf (stderr, "Error, foo(%la) should be +Inf, got %la=%"PRIx64"\n",
               x, y, asuint64 (y));
#ifndef DO_NOT_ABORT
      exit (1);
#endif
    }
    // check the invalid exception is not set
    flag = fetestexcept (FE_INVALID);
    if (flag)
    {
      printf ("Spurious invalid exception for x=%la\n", x);
#ifndef DO_NOT_ABORT
      exit (1);
#endif
    }

    // Check -x
    feclearexcept (FE_INVALID);
    y = cr_sinh (-x);
    if (x >= 0x1.633ce8fb9f87ep+9 && !is_inf (y))
    {
      fprintf (stderr, "Error, foo(%la) should be -Inf, got %la=%"PRIx64"\n",
               x, y, asuint64 (y));
#ifndef DO_NOT_ABORT
      exit (1);
#endif
    }
    // check the invalid exception is not set
    flag = fetestexcept (FE_INVALID);
    if (flag)
    {
      printf ("Spurious invalid exception for x=%la\n", x);
#ifndef DO_NOT_ABORT
      exit (1);
#endif
    }  
  }
}

// put in h+l a double-double approximation of sinh(x)
static void
dd_sinh (double *h, double *l, double x)
{
  mpfr_t t;
  mpfr_init2 (t, 107);
  mpfr_set_d (t, x, MPFR_RNDN);
  mpfr_sinh (t, t, MPFR_RNDN);
  *h = mpfr_get_d (t, MPFR_RNDN);
  mpfr_sub_d (t, t, *h, MPFR_RNDN);
  *l = mpfr_get_d (t, MPFR_RNDN);
  mpfr_clear (t);
}

// test |n| values starting from x (downwards if n < 0)
static void scan_consecutive(int64_t n, double x){
  ref_init();
  ref_fesetround(rnd);
  fesetround(rnd1[rnd]);
  if (n < 0) {
    n = -n;
    x = asfloat64 (asuint64 (x) - n);
  }
  int64_t n0 = n;
  while (n) {
    double h, l, d, dd;
    dd_sinh (&h, &l, x);
    d = sqrt (1.0 + h * h); // derivative cosh(x) = sqrt(1+sinh(x)^2)
    dd = fabs (h); // absolute value of 2nd derivative
    int e;
    frexp (x, &e);
    /* 2^(e-1) <= |x| < 2^e thus ulp(x) = 2^(e-53) */
    d = ldexp (d, e - 53); // multiply d by ulp(x)
    dd = ldexp (dd, 2 * (e - 53)); // multiply dd by ulp(x)^2
    /* we want j^2*dd < 2^-11 ulp(h) so that the 2nd-order term
       produces an error bounded by 2^-11 ulp(h), to that MPFR
       will be called with probability about 2^-11.
       Thus approximately j^2*dd < 2^-64 h,
       or j < 2^-32 sqrt(h/dd) */
    int64_t jmax = 0x1p-32 * sqrt (fabs (h) / dd);
    if (jmax > n) jmax = n; // cap to n
    assert (jmax > 0);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
    for(int64_t j=0;j<jmax;j++){
      b64u64_u v = {.f = x};
      v.u += j;
      double t = tfun (v.f);
      // asinh(x+j*u) is approximated by h + l + j*d
      double w = h + __builtin_fma (j, d, l);
      if (t != w) // expensive test
        check(v.f);
    }
    n -= jmax;
    x += jmax * ldexp (1.0, e - 53);
  }
  printf ("checked %lu values, expensive checks %lu\n", n0, tested);
}

int
main (int argc, char *argv[])
{
  int conseq = 0;  // scan consecutive values
  double a;        // starting value for scan_consecutive
  unsigned long C; // length for scan_consecutive
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
      else if (strcmp (argv[1], "-C") == 0)
        {
          conseq = 1;
          C = strtoul (argv[2], NULL, 0);
          argc -= 2;
          argv += 2;
        }
      else if (strcmp (argv[1], "-a") == 0)
        {
          conseq = 1;
          a = strtod (argv[2], NULL);
          argc -= 2;
          argv += 2;
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

  if (conseq) {
    scan_consecutive (C, a);
    return 0;
  }

  check_invalid ();

#ifndef CORE_MATH_TESTS
#define CORE_MATH_TESTS 1000000000UL /* total number of tests */
#endif

  unsigned int seed = getpid ();
  for (int i = 0; i < MAX_THREADS; i++)
    Seed[i] = seed + i;

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
    double x;
    do x = get_random (tid); while (fabs (x) >= 0x1.633ce8fb9f87ep+9);
    check (x);
  }

  return 0;
}
