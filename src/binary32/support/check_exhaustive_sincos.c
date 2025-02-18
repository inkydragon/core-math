/* Check correctness of binary32 function like sincos by exhaustive search.

Copyright (c) 2022 Alexei Sibidanov.
Copyright (c) 2022-2024 Paul Zimmermann, INRIA.

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
#include <mpfr.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

#include "function_under_test.h"

void cr_function_under_test (float, float*, float*);
void ref_function_under_test (float, float*, float*);
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

static int
is_equal (float y1, float y2)
{
  if (isnan (y1))
    return isnan (y2);
  if (isnan (y2))
    return isnan (y1);
  return asuint (y1) == asuint (y2);
}

void
doit (uint32_t n)
{
  float x, y1, y2, z1, z2;
  x = asfloat (n);
  ref_init ();
  ref_fesetround (rnd);
  mpfr_flags_clear (MPFR_FLAGS_INEXACT);
  ref_function_under_test (x, &y1, &y2);
#ifdef CORE_MATH_CHECK_INEXACT
  mpfr_flags_t inex_y = mpfr_flags_test (MPFR_FLAGS_INEXACT);
#endif
  fesetround (rnd1[rnd]);
  feclearexcept (FE_INEXACT);
  cr_function_under_test (x, &z1, &z2);
  int inex_z = fetestexcept (FE_INEXACT);
  if (!is_equal (y1, z1) || !is_equal (y2, z2))
  {
    printf ("FAIL x=%a ref=(%a,%a) z=(%a,%a)\n", x, y1, y2, z1, z2);
    fflush (stdout);
    if (!keep) exit (1);
  }
#ifdef CORE_MATH_CHECK_INEXACT
  if ((inex_y == 0) && (inex_z != 0))
  {
    printf ("Spurious inexact exception for x=%a\n", x);
    fflush (stdout);
    if (!keep) exit (1);
  }
  if ((inex_y != 0) && (inex_z == 0))
  {
    printf ("Missing inexact exception for x=%a\n", x);
    fflush (stdout);
    if (!keep) exit (1);
  }
#endif
}

// When x is a NaN, returns 1 if x is an sNaN and 0 if it is a qNaN
static inline int issignaling(float x) {
  union_t _x = {.x = x};

  return !(_x.n & (1ull << 22));
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
is_nan (float x)
{
  uint32_t u = asuint (x);
  int e = u >> 23;
  return (e == 0xff || e == 0x1ff) && (u << 9) != 0;
}

/* check for signaling NaN input */
static void
check_signaling_nan (void)
{
  float snan = asfloat (0x7f800001);
  float y, z;
  cr_function_under_test (snan, &y, &z);
  // check that y = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, 1st return value should be NaN, got %la=%x\n",
             y, asuint (y));
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, 1st return value should be qNaN, got sNaN=%x\n",
             asuint (y));
    exit (1);
  }
  // check that z = NaN
  if (!is_nan (z))
  {
    fprintf (stderr, "Error, 2nd return value should be NaN, got %la=%x\n",
             z, asuint (z));
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (z))
  {
    fprintf (stderr, "Error, 2nd return value should be qNaN, got sNaN=%x\n",
             asuint (z));
    exit (1);
  }
  // also test with the sign bit set
  snan = asfloat (0xff800001);
  cr_function_under_test (snan, &y, &z);
  // check that y = NaN
  if (!is_nan (y))
  {
    fprintf (stderr, "Error, 1st return value should be NaN, got %la=%x\n",
             y, asuint (y));
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (y))
  {
    fprintf (stderr, "Error, 1st return value should be qNaN, got sNaN=%x\n",
             asuint (y));
    exit (1);
  }
  // check that z = NaN
  if (!is_nan (z))
  {
    fprintf (stderr, "Error, 2nd return value should be NaN, got %la=%x\n",
             z, asuint (z));
    exit (1);
  }
  // check that the signaling bit disappeared
  if (issignaling (z))
  {
    fprintf (stderr, "Error, 2nd return value should be qNaN, got sNaN=%x\n",
             asuint (z));
    exit (1);
  }
}

static inline int doloop (void)
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

  // check regular numbers
  uint32_t nmin = asuint (0x0p0f), nmax = asuint (0x1.fffffep+127);
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
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
