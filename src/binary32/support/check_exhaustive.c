/* Check correctness of binary32 function by exhaustive search.

Copyright (c) 2022 Alexei Sibidanov.
Copyright (c) 2022 Paul Zimmermann, INRIA.

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
#include <mpfr.h>
#include <errno.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

#include "function_under_test.h"

float cr_function_under_test (float);
float ref_function_under_test (float);
int ref_fesetround (int);
void ref_init (void);

/* the code below is to check correctness by exhaustive search */

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int keep = 0;

typedef union { uint32_t n; float x; } union_t;

float
asfloat (uint32_t n)
{
  union_t u;
  u.n = n;
  return u.x;
}

static inline uint32_t
asuint (float f)
{
  union
  {
    float f;
    uint32_t i;
  } u = {f};
  return u.i;
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
is_nan (float x)
{
  uint32_t u = asuint (x);
  int e = u >> 23;
  return (e == 0xff || e == 0x1ff) && (u << 9) != 0;
}

/* define our own is_inf function to avoid depending from math.h */
static inline int
is_inf (float x)
{
  uint32_t u = asuint (x);
  int e = u >> 23;
  return (e == 0xff || e == 0x1ff) && (u << 9) == 0;
}

static int
is_equal (float y1, float y2)
{
  if (is_nan (y1))
    return is_nan (y2);
  if (is_nan (y2))
    return is_nan (y1);
  return asuint (y1) == asuint (y2);
}

void
doit (uint32_t n)
{
  float x, y, z;
  x = asfloat (n);
  ref_init ();
  ref_fesetround (rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  y = ref_function_under_test (x);
#ifdef CORE_MATH_CHECK_INEXACT
  mpfr_flags_t inex_y = mpfr_flags_test (MPFR_FLAGS_INEXACT);
#endif
  fesetround (rnd1[rnd]);
  feclearexcept (FE_INEXACT | FE_UNDERFLOW);
  errno = 0;
  z = cr_function_under_test (x);
  int inex_z = fetestexcept (FE_INEXACT);
  /* Note: the test y != z would not distinguish +0 and -0, instead we compare
     the 32-bit encodings. */
  if (!is_equal (y, z))
  {
    printf ("FAIL x=%a ref=%a y=%a\n", x, y, z);
    fflush (stdout);
    if (!keep) exit (1);
  }
  // check spurious underflow
  if ((y < -0x1p-126f || 0x1p-126f < y) && fetestexcept (FE_UNDERFLOW))
  {
    printf ("Spurious underflow exception for x=%a (y=%a)\n", x, y);
    fflush (stdout);
    if (!keep) exit (1);
  }
#ifdef CORE_MATH_CHECK_INEXACT
  if ((inex_y == 0) && (inex_z != 0))
  {
    printf ("Spurious inexact exception for x=%a (y=%a)\n", x, y);
    fflush (stdout);
    if (!keep) exit (1);
  }
  if ((inex_y != 0) && (inex_z == 0))
  {
    printf ("Missing inexact exception for x=%a (y=%a)\n", x, y);
    fflush (stdout);
    if (!keep) exit (1);
  }
#endif
#ifdef CORE_MATH_SUPPORT_ERRNO
  /* If x is a normal number and y is NaN, we should have errno = EDOM.
     If x is a normal number and y is +/-Inf, we should have errno = ERANGE.
  */
  if (!is_nan (x) && !is_inf (x))
  {
    if (is_nan (y) && errno != EDOM)
    {
      printf ("Missing errno=EDOM for x=%a (y=%a)\n", x, y);
      fflush (stdout);
      if (!keep) exit (1);
    }
    if (is_inf (y) && errno != ERANGE)
    {
      printf ("Missing errno=ERANGE for x=%a (y=%a)\n", x, y);
      fflush (stdout);
      if (!keep) exit (1);
    }
  }
#endif
}

// When x is a NaN, returns 1 if x is an sNaN and 0 if it is a qNaN
static inline int issignaling(float x) {
  union_t _x = {.x = x};

  return !(_x.n & (1ull << 22));
}

/* check for signaling NaN input */
static void
check_signaling_nan (void)
{
  float snan = asfloat (0x7f800001ul);
  float y = cr_function_under_test (snan);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(sNaN) should be NaN, got %la=%x\n",
             y, asuint (y));
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(sNaN) should be qNaN, got sNaN=%x\n",
             asuint (y));
    exit (1);
  }
  // also test sNaN with sign bit set
  snan = asfloat (0xff800001ul);
  y = cr_function_under_test (snan);
  // check that foo(NaN) = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, foo(sNaN) should be NaN, got %la=%x\n",
             y, asuint (y));
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, foo(sNaN) should be qNaN, got sNaN=%x\n",
             asuint (y));
    exit (1);
  }
}

static void
check_exceptions_aux (uint32_t n)
{
  float x = asfloat (n);
  feclearexcept (FE_INEXACT);
  float y = cr_function_under_test (x);
  int inex = fetestexcept (FE_INEXACT);
  // there should be no inexact exception if the result is NaN, +/-Inf or +/-0
  if (inex && (is_nan (y) || is_inf (y) || y == 0))
  {
    fprintf (stderr, "Error, for x=%a=%x, inexact exception set (y=%a=%x)\n",
             x, asuint (x), y, asuint (y));
    exit (1);
  }
  feclearexcept (FE_OVERFLOW);
  y = cr_function_under_test (x);
  inex = fetestexcept (FE_OVERFLOW);
  if (inex)
  {
    fprintf (stderr, "Error, for x=%a, overflow exception set (y=%a)\n", x, y);
    exit (1);
  }
  feclearexcept (FE_UNDERFLOW);
  y = cr_function_under_test (x);
  inex = fetestexcept (FE_UNDERFLOW);
  if (inex)
  {
    fprintf (stderr, "Error, for x=%a, underflow exception set (y=%a)\n", x, y);
    exit (1);
  }
}

// check that no overflow/underflor/inexact is set for input NaN, Inf, 0
// when the output is also NaN, Inf, 0
static void
check_exceptions (void)
{
  // check +sNaN and -sNaN
  check_exceptions_aux (0x7f800001);
  check_exceptions_aux (0xff800001);
  // check +qNaN and -qNaN
  check_exceptions_aux (0x7fc00000);
  check_exceptions_aux (0xffc00000);
  // check +Inf and -Inf
  check_exceptions_aux (0x7f800000);
  check_exceptions_aux (0xff800000);
  // check +0 and -0
  check_exceptions_aux (0x0);
  check_exceptions_aux (0x80000000);
}

static int doloop (void)
{
  // check sNaN
  doit (0x7f800001);
  doit (0xff800001);
  // check qNaN
  doit (0x7fc00000);
  doit (0xffc00000);
  // check +Inf and -Inf
  doit (0x7f800000);
  doit (0xff800000);

  check_signaling_nan ();

  check_exceptions ();

  // check regular numbers
  uint32_t nmin = asuint (0x0p0f), nmax = asuint (0x1.fffffep+127f);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
  /* Use a static schedule with small chunks, since the function might be
     very easy to evaluate in some ranges, for example log of x < 0. */
#pragma omp parallel for schedule(static,1024)
#endif
  for (uint32_t n = nmin; n <= nmax; n++)
  {
    doit (n);
    doit (n | 0x80000000);
  }
  printf ("all ok\n");
  return 0;
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
      else if (strcmp (argv[1], "--keep") == 0)
        {
          keep = 1;
          argc --;
          argv ++;
        }
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[1]);
          exit (1);
        }
    }

  return doloop();
}
