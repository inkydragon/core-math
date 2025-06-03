/* Check correctness of univariate binary32 function by exhaustive search.

Copyright (c) 2022-2025 Alexei Sibidanov and Paul Zimmermann

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

_Float16 cr_function_under_test (_Float16);
_Float16 ref_function_under_test (_Float16);
int mpfr_function_under_test (mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
int ref_fesetround (int);
void ref_init (void);

/* the code below is to check correctness by exhaustive search */

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };
int rnd2[] = { MPFR_RNDN,    MPFR_RNDZ,     MPFR_RNDU, MPFR_RNDD   };

int rnd = 0;
int keep = 0;

typedef union { uint16_t n; _Float16 x; } union_t;

float
asfloat (uint16_t n)
{
  union_t u;
  u.n = n;
  return u.x;
}

static inline uint16_t
asuint (_Float16 f)
{
  union
  {
    _Float16 f;
    uint16_t i;
  } u = {f};
  return u.i;
}

/* define our own is_nan function to avoid depending from math.h */
static inline int
is_nan (_Float16 x)
{
  uint16_t u = asuint (x);
  int e = u >> 10;
  return (e == 0x1f || e == 0x3f) && (u << 6) != 0;
}

/* define our own is_inf function to avoid depending from math.h */
static inline int
is_inf (_Float16 x)
{
  uint16_t u = asuint (x);
  int e =  u >> 10;
  return (e == 0x1f || e == 0x3f) && (u << 6) == 0;
}

static int
is_equal (_Float16 y1, _Float16 y2)
{
  if (is_nan (y1))
    return is_nan (y2);
  if (is_nan (y2))
    return is_nan (y1);
  return asuint (y1) == asuint (y2);
}

// IL MANQUE PLEIN DE CHOSES :
// 	underflow_before
//	overflow
//	inexact
//	errno
//	is_signaling

void
doit (uint16_t n)
{
  _Float16 x, y, z;
  x = asfloat (n);
  ref_init ();
  ref_fesetround (rnd);
  		// mpfr_flags_clear (MPFR_FLAGS_INEXACT | MPFR_FLAGS_UNDERFLOW | MPFR_FLAGS_OVERFLOW);
  y = ref_function_under_test (x);
		// #if defined(CORE_MATH_CHECK_INEXACT) || defined(CORE_MATH_SUPPORT_ERRNO)
  		// mpfr_flags_t inex_y = mpfr_flags_test (MPFR_FLAGS_INEXACT);
		// #endif
  fesetround (rnd1[rnd]);
  		// feclearexcept (FE_INEXACT | FE_UNDERFLOW | FE_OVERFLOW);
		// #ifdef CORE_MATH_SUPPORT_ERRNO
		//  errno = 0;
		// #endif
  z = cr_function_under_test (x);
		// #ifdef CORE_MATH_CHECK_INEXACT
  		// int inex_z = fetestexcept (FE_INEXACT);
		// #endif
  		/* Note: the test y != z would not distinguish +0 and -0, instead we compare
    		 the 32-bit encodings. */
  if (!is_equal (y, z))
  {
    printf ("FAIL x=%a ref=%a y=%a\n", (float) x, (float) y, (float) z);
    fflush (stdout);
    if (!keep) exit (1);	
  }
  								// BOUT DE CODE SUPPRIMÉ
}

static int doloop (void)
{
  // check regular numbers
  uint16_t nmin = 0x0001, nmax = 0x7bff;
		// #if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
  		/* Use a static schedule with small chunks, since the function might be
     		very easy to evaluate in some ranges, for example log of x < 0. */
		// #pragma omp parallel for schedule(static,1024)
		// #endif
  for (uint16_t n = nmin; n <= nmax; n++)
  {
    doit (n);
    doit (n | 0x8000);
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

  // check_underflow_before ();

  return doloop();
}
