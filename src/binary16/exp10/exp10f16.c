/* Correctly-rounded 10^x function for binary16 value.

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

#define _GNU_SOURCE // needed to define expf10
#include <stdint.h>
#include <errno.h>
#include <fenv.h>
#include <math.h> // only used during performance tests

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_exp10f16(_Float16 x){
	static const b16u16_u x0 = {.f = -0x1.e18p+2f}; // smallest _Float16 such that 10^x0- < MIN_SUBNORMALIZE <= 10^x0
 	static const b16u16_u x1 = {.f = 0x1.344p+2f}; // largest _Float16 such that 10^x1 <= MAX_FLOAT16 < 10^x1+
	static const float tb[] = // tabulate value of 10^(log10(2)*i/64) for i in [0, 63]
		{0x1p+0f, 0x1.02c9a4p+0f, 0x1.059b0ep+0f, 0x1.087452p+0f,  
		 0x1.0b5586p+0f, 0x1.0e3ec4p+0f, 0x1.11301ep+0f, 0x1.1429aap+0f,  
		 0x1.172b84p+0f, 0x1.1a35bep+0f, 0x1.1d4874p+0f, 0x1.2063b8p+0f,  
		 0x1.2387a6p+0f, 0x1.26b456p+0f, 0x1.29e9ep+0f, 0x1.2d285ap+0f,  
		 0x1.306fep+0f, 0x1.33c08cp+0f, 0x1.371a74p+0f, 0x1.3a7db4p+0f,
		 0x1.3dea68p+0f, 0x1.4160a2p+0f, 0x1.44e086p+0f, 0x1.486a2cp+0f,  
		 0x1.4bfdaep+0f, 0x1.4f9b28p+0f, 0x1.5342b6p+0f, 0x1.56f474p+0f,  
		 0x1.5ab07ep+0f, 0x1.5e76f2p+0f, 0x1.6247ecp+0f, 0x1.662388p+0f,  
		 0x1.6a09dep+0f, 0x1.6dfb24p+0f, 0x1.71f75ep+0f, 0x1.75feb6p+0f,  
		 0x1.7a1148p+0f, 0x1.7e2f34p+0f, 0x1.82589ap+0f, 0x1.868d9ap+0f,  
		 0x1.8ace54p+0f, 0x1.8f1aeap+0f, 0x1.93737cp+0f, 0x1.97d82ap+0f,  
		 0x1.9c4918p+0f, 0x1.a0c668p+0f, 0x1.a5503cp+0f, 0x1.a9e6b6p+0f,  
		 0x1.ae89fap+0f, 0x1.b33a2cp+0f, 0x1.b7f77p+0f, 0x1.bcc1eap+0f,  
		 0x1.c199bep+0f, 0x1.c67f12p+0f, 0x1.cb720ep+0f, 0x1.d072d4p+0f,  
		 0x1.d5818ep+0f, 0x1.da9e6p+0f, 0x1.dfc974p+0f, 0x1.e502eep+0f,  
		 0x1.ea4af6p+0f, 0x1.efa1cp+0f, 0x1.f50766p+0f, 0x1.fa7c1ap+0f};
	b16u16_u v = {.f = x};
#ifdef CORE_MATH_SUPPORT_ERRNO
	if (v.f > x1.f || v.f <= -0x1.0dcp+2)
		errno = ERANGE;
	if (v.f == x1.f && 0x1.fffffep127f + 0x1p104f == 0x1.fffffep127f)
		errno = ERANGE;
#endif
	if ((v.u & 0x7fff) > 0x44d1) { // in this case, we have x > min(x0, x1) in abs value
		if ((v.u & 0x7c00) == 0x7c00) { // if x is nan or x is inf
			if (v.u == 0xfc00) return 0x0p0;
			else return x + x;
		}
		else if (v.u > x0.u) return 0x1p-25f; // x smaller than x0
		else if (v.f > x1.f) return 0x1.ffcp15f + 0x1p5f; // x greater than x1
	}
	// 7 wrong cases
	if (v.u == 0xbf0a) return 0x1.1ccp-6f + 0x1p-18f;
	if (v.u == 0x9766) return 0x1.fdcp-1f + 0x1p-13f;
	if (v.u == 0xb2d4) return 0x1.394p-1f - 0x1p-13f;
	if (!(v.u & 0x81ff)) { // exact cases (x = 1, 2, 3, 4)
		if (v.u == 0x4000) return 0x1.9p+6f;
		if (v.u == 0x4200) return 0x1.f4p+9f;
		if (v.u == 0x4400) return 0x1.388p+13f;
		if (v.u == 0x3c00) return 0x1.4p+3f;	
	}
	float xf = x; // exact conversion from _Float16 to float
	static const float sixtyfour_over_log10_2 = 0x1.a934fp+7f; 
	static const float minus_log10_2_over_sixtyfour = -0x1.344136p-8f;
	float j = __builtin_roundevenf(sixtyfour_over_log10_2 * xf);
	int32_t jint = j;
	int i = jint & 0x3f;
	float xp = __builtin_fmaf(minus_log10_2_over_sixtyfour, j, xf);
	// xf = j*log10(2)/64 + xp = (j>>6)*log10(2) + log10(2)*i/64 + xp
	// so 10^xf = 2^(j>>6) * 10^(log10(2)*i/64) * 10^xp
	xp = __builtin_fmaf(__builtin_fmaf(0x1.53524ep1f, xp, 0x1.26bb1cp1f), xp, 1.0f);; 
	b32u32_u w = {.f = xp * tb[i]};
	w.u += (jint >> 6) * (1l << 23);
	return w.f; // conversion float -> _Float16 (with rounding)
}

// dummy function since GNU libc does not provide it
_Float16 exp10f16 (_Float16 x) {
  return (_Float16) exp10f ((float) x);
}
