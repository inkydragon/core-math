/* Correctly-rounded reciprocal square root function for binary32 value.

Copyright (c) 2022 Alexei Sibidanov.

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

#include <mpfr.h>
#include "fenv_mpfr.h"

/* code from MPFR */

float
ref_rsqrt (float x)
{
  /* mpfr_rec_sqrt differs from IEEE 754-2019: IEEE 754-2019 says that
     rsqrt(-0) should give -Inf, whereas mpfr_rec_sqrt(-0) gives +Inf */
  if (x == 0.0f && 1.0f / x < 0.0f)
    return 1.0f / x;
  mpfr_t y;
  mpfr_init2 (y, 24);
  mpfr_set_flt (y, x, MPFR_RNDN);
  mpfr_rec_sqrt (y, y, rnd2[rnd]);
  /* since |x| < 2^128 for non-zero x, we have 2^-64 < 1/sqrt(x),
     thus no underflow can happen, and there is no need to call
     mpfr_subnormalize */
  float ret = mpfr_get_flt (y, MPFR_RNDN);
  mpfr_clear (y);
  return ret;
}
