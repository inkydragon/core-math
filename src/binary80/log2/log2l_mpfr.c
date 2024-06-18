/* Correctly-rounded log2l of binary64 value.

Copyright (c) 2024 Paul Zimmermann (Inria)

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

/* reference function using MPFR */
long double
ref_log2l (long double x)
{
  mpfr_t y;
  mpfr_exp_t emin = mpfr_get_emin ();
  mpfr_set_emin (-16444);
  mpfr_init2 (y, 64);
  mpfr_set_ld (y, x, MPFR_RNDN);
  mpfr_log2 (y, y, rnd2[rnd]);
  // no need to call mpfr_subnormalize
  long double ret = mpfr_get_ld (y, MPFR_RNDN);
  mpfr_clear (y);
  mpfr_set_emin (emin);
  return ret;
}
