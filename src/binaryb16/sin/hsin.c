/* Correctly-rounded sine function for bfloat16 value.

Copyright (c) 2025 Paul Zimmermann

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

#include <stdint.h>
#include <errno.h>
#include <math.h> // only used during performance tests
#include <stdio.h>
#include <fenv.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {__bf16 f; uint16_t u;} b16u16_u;

#include "/tmp/S1.h"
#include "/tmp/C1.h"
#include "/tmp/S2.h"
#include "/tmp/C2.h"

__bf16 cr_hsin (__bf16 x){
  // int bug = x == 0x1.d4p-4;
  int bug = 0;
  b16u16_u v = {.f = x};
  uint16_t u = v.u, au = u & 0x7fff;

  static const float s[] = {1.0, -1.0};

  // for |x| <= 0x1.dp-4, sin(x) rounds to x to nearest
  if (au <= 0x3de8) { // |x| <= 0x1.dp-4
    if (au == 0) return x; // x = +0 or -0
    // avoid spurious underflow for |x|=0x1p-126
    if (au == 0x80) return s[u>>15] * 0x1.fffffffffffffp-127;
    float res = x * 0x1.fffffep-1f;
    return res;
  }

  // now 0x3de8 < au, |x| > 0x1.dp-4
  
  if (au < 0x4580) { // |x| < 4096
    // the exponent part lies in [123,138]
    uint16_t i1 = (au - 0x3d80) >> 3; // 0x3d80 is 123<<7
    uint16_t i2 = (((au - 0x3d80) >> 7) << 3) | (au & 0x7);
    // we use a FMA to fix the evaluation order
    if (bug) printf ("i1=%d i2=%d %a %a %a %a\n", i1, i2,
                     S1[i1], C2[i2], C1[i1], S2[i2]);
    float res = __builtin_fmaf (S1[i1], C2[i2], C1[i1] * S2[i2]);
    if (bug) printf ("res=%a\n", res);
    return s[u>>15] * res;
  }

  if ((au >> 7) == 0xff) { // NaN or Inf
    b16u16_u qnan = {.u = 0x7fc0};
    return qnan.f + x;
  }

  // now |x| >= 4096
  return (__bf16) sinf ((float) x);
}

// dummy function since GNU libc does not provide it
__bf16 hsin (__bf16 x) {
  return (__bf16) sinf ((float) x);
}
