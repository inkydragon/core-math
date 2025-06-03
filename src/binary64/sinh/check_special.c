/* Check sinh on random inputs.

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
#include <inttypes.h>
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
    exit (1);
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
static inline int issignaling(double x) {
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
  if (!is_inf(y) && y > 0)
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
  if (!is_inf(y) && y < 0)
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
  if (!is_nan (y) || issignaling(y))
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
