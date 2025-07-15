/* Correctly-rounded natural exponential function for binary16 value.

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
#include <math.h> // only used during performance tests

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_expf16(_Float16 x){
	static const b32u32_u x0 = {.f = -0x1.154p+4f}; // smallest _Float16 such that exp(x0-) < MIN_SUBNORMALIZE <= exp(x0)
 	static const b32u32_u x1 = {.f = 0x1.62cp3f}; // largest _Float16 such that exp(x1) <= MAX_FLOAT16 < exp(x1+)
	static const float tb[] = // tabulate value of exp(log(2)*i/32) for i in [0, 31]
		{0x1p+0f, 0x1.059b0cp+0f, 0x1.0b5586p+0f, 0x1.11301ep+0f,  
		 0x1.172b84p+0f, 0x1.1d4874p+0f, 0x1.2387a6p+0f, 0x1.29e9ep+0f,  
		 0x1.306fep+0f, 0x1.371a74p+0f, 0x1.3dea64p+0f, 0x1.44e086p+0f,  
		 0x1.4bfdaep+0f, 0x1.5342b6p+0f, 0x1.5ab07ep+0f, 0x1.6247ecp+0f,  
		 0x1.6a09e6p+0f, 0x1.71f75ep+0f, 0x1.7a1148p+0f, 0x1.82589ap+0f,  
		 0x1.8ace54p+0f, 0x1.93737cp+0f, 0x1.9c4918p+0f, 0x1.a5503cp+0f,  
		 0x1.ae89fap+0f, 0x1.b7f77p+0f, 0x1.c199bep+0f, 0x1.cb720ep+0f,  
		 0x1.d5819p+0f, 0x1.dfc974p+0f, 0x1.ea4afap+0f, 0x1.f50766p+0f};
	b32u32_u v = {.f = x};
#ifdef CORE_MATH_SUPPORT_ERRNO
	if (v.f > x1.f || v.f < -0x1.368p+3f)
		errno = ERANGE;
#endif
	if ((v.u & 0x7fffffff) > 0x41316000) { // in this case, we have x > min(x0, x1) in abs value
		if ((v.u & 0x7f800000) == 0x7f800000) { // if x is nan or x is inf
			if (v.u == 0xff800000) return 0x0p0;
			else return x + x;
		}
                /* With -DCORE_MATH_SUPPORT_ERRNO, gcc 14.2.0 emits a spurious
                   underflow for x=0x1.63p+3 (for example). This is due
                   to https://gcc.gnu.org/bugzilla/show_bug.cgi?id=120910. */
		else if (v.u > x0.u) return 0x1p-25f; // x smaller than x0
		else if (v.f > x1.f) return 0x1.ffcp15f + 0x1p5f; // x greater than x1
	}
	static const float thirtytwo_over_log2 = 0x1.715476p+5f;
	static const float minus_log2_over_thirtytwo = -0x1.62e430p-6f;
	float j = __builtin_roundevenf(thirtytwo_over_log2 * v.f);
	int32_t jint = j;
	int i = jint & 0x1f;
	float xp = __builtin_fmaf(minus_log2_over_thirtytwo, j, v.f);
	// xf = j*log(2)/32 + xp = (j>>5)*log(2) + log(2)*i/32 + xp
	// so exp(xf) = 2^(j>>5) * exp(log(2)*i/32) * exp(xp)
	xp = __builtin_fmaf(__builtin_fmaf(0.5f, xp, 1.0f), xp, 1.0f);
	b32u32_u w = {.f = xp * tb[i]};
	w.u += (jint >> 5) * (1l << 23);
	return w.f; // conversion float -> _Float16 (with rounding)
}

// dummy function since GNU libc does not provide it
_Float16 expf16 (_Float16 x) {
  return (_Float16) expf ((float) x);
}
