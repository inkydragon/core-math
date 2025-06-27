/* Correctly-rounded 2^x function for binary16 value.

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

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_exp2f16(_Float16 x){
	static const b16u16_u x0 = {.f = -0x1.9p+4}; // smallest _Float16 such that 2^x0- < MIN_SUBNORMALIZE <= 2^x0
 	static const b16u16_u x1 = {.f = 0x1.ffcp+3f}; // largest _Float16 such that 2^x1 <= MAX_FLOAT16 < 2^x1+
	static const float tb[] = // tabulate value of 2^(i/32) for i in [0, 31]
		{0x1p+0f, 0x1.059b0ep+0f, 0x1.0b5586p+0f, 0x1.11301ep+0f,  
		 0x1.172b84p+0f, 0x1.1d4874p+0f, 0x1.2387aap+0f, 0x1.29e9ep+0f,  
		 0x1.306fep+0f, 0x1.371a74p+0f, 0x1.3dea64p+0f, 0x1.44e086p+0f,  
		 0x1.4bfdaep+0f, 0x1.5342b6p+0f, 0x1.5ab07ep+0f, 0x1.6247ecp+0f,  
		 0x1.6a09e6p+0f, 0x1.71f75ep+0f, 0x1.7a1148p+0f, 0x1.82589ap+0f,  
		 0x1.8ace54p+0f, 0x1.93737cp+0f, 0x1.9c4918p+0f, 0x1.a5503cp+0f,  
		 0x1.ae89fap+0f, 0x1.b7f77p+0f, 0x1.c199bep+0f, 0x1.cb720ep+0f,  
		 0x1.d5818cp+0f, 0x1.dfc976p+0f, 0x1.ea4afap+0f, 0x1.f5076ap+0f};
	b16u16_u v = {.f = x};
	if ((v.u & 0x7fff) > 0x4bff) { // in this case, we have x > min(x0, x1) in abs value
		if ((v.u & 0x7c00) == 0x7c00) { // if x is nan or x is inf
			if (v.u == 0xfc00) return 0x0p0;
			else return x + x;
		}
		else if (v.u > x0.u) return 0x1p-25f; // x smaller than x0
		else if (v.f > x1.f) return 0x1.ffcp15f + 0x1p5f; // x greater than x1
	}
	float xf = x; // exact conversion from _Float16 to float
	static const float thirtytwo = 0x1p5f;
	static const float minus_one_over_thirtytwo = -0x1p-5f;
	float j = __builtin_roundevenf(thirtytwo * xf);
	uint32_t jint = j;
	int i = jint & 0x1f;
	if (i == 0) { // two only wrong cases
		if (v.u == 0x11c5) return 0x1.004p+0f - 0x1p-12f;
		if (v.u == 0xa39d) return 0x1.facp-1f - 0x1p-13f;
	}
	float xp = __builtin_fmaf(minus_one_over_thirtytwo, j, xf);
	// xf = j/32 + xp = (j>>5) + i/32 + xp
	// so exp(xf) = 2^(j>>5) * 2^(i/32) * 2^xp
	xp = __builtin_fmaf(__builtin_fmaf(0x1.ebfbep-3f, xp, 0x1.62e43p-1f), xp, 1.0f);
	b32u32_u w = {.f = xp * tb[i]};
	w.u += (jint >> 5) * (1l << 23);
	return w.f; // conversion float -> _Float16 (with rounding)
}

// dummy function since GNU libc does not provide it
_Float16 exp2f16 (_Float16 x) {
  return (_Float16) exp2f ((float) x);
}
