/* Correctly-rounded reciprocal square root function for binary16 value.

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
#include <errno.h>
#include <math.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON


typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_rsqrtf16(_Float16 x){
	b16u16_u t = {.f = x};
	// two types of wrong cases : 0xxxx11100010011 and 0xxxx10100011111
	if ((t.u | 0x7800) == 0x7f13) {
		if (t.u == 0x7f13) return 0.0f / 0.0f; // if x=sNaN return NaN
		int expo = -((t.u >> 10) - 15) / 2 + 14;
		t.u = (expo << 10) + 0x204;
		return -0x1p-20 + t.f;
	} else if ((t.u | 0x7800) == 0x7d1f) {
		if (t.u == 0x7d1f) return 0.0f / 0.0f; // if x=sNaN return NaN
		int expo = -((t.u >> 10) - 15) / 2 + 14;
		t.u = (expo << 10) + 0x312;
		return 0x1p-20 + t.f;
	}
#ifdef CORE_MATH_SUPPORT_ERRNO
	if (x == 0.0f)
		errno = ERANGE;
#endif
	return 1.0f / sqrtf ((float) x);
}

// dummy function since GNU libc does not provide it
_Float16 rsqrtf16 (_Float16 x) {
	return (_Float16) (1.0f / sqrtf ((float) x));
}
