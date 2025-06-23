/* Correctly-rounded cubuc root of binary16 value.

Copyright (c) 2025 Maxence Ponsardin.

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
#include <math.h> // only used during performance tests
#include <fenv.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON


typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_cbrtf16(_Float16 x){
	b16u16_u t = {.f = x};
	if ((t.u & 0x03ff) == 0x0151) { // only wrong case is 0x1.544pk with k = 1 mod 3
		int expo = ((t.u & 0x7fff) >> 10);
		if (expo % 3 == 1 && expo < 31) { // avoid sNaN and k != 1 mod 3
			t.u = (((expo - 16) / 3 + 15) << 10) + 0x018b + ((t.u >> 15) << 15);
			return (float) t.f + 0x1p-16 * ((t.u >> 15) - 0.5f);
		}
	}
	b32u32_u y = {.f = cbrtf ((float) x)};
	// exact cases have always 7 ending 0
	if (!(y.u % (1 << 20)) && y.f * y.f * y.f == x) feclearexcept(FE_INEXACT);
	return y.f;
}

// dummy function since GNU libc does not provide it
_Float16 cbrtf16 (_Float16 x) {
  return (_Float16) cbrtf ((float) x);
}
