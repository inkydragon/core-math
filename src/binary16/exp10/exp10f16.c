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
#include <math.h> // only used during performance tests

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_exp10f16(_Float16 x){
	uint16_t x0 = 0xc786; // binary representation of x0 in order to compare uint16_t rather than flt16
 	_Float16 x1 = 0x1.344p2; // largest _Float16 such that 10^x1 <= MAX_FLOAT16 < 10^x1+
	static const float tb[] = // tabulate value of 10^(i/2^6) for i in [-2^5*log(2)/log(10), 2^5*log(2)/log(10)], size(tb) = 39
		{0x8.13b01p-4, 0x8.5f6eep-4, 0x8.adf41p-4, 0x8.ff59ap-4,
		 0x9.53ba8p-4, 0x9.ab32bp-4, 0xa.05df3p-4, 0xa.63dep-4,  
		 0xa.c54e5p-4, 0xb.2a506p-4, 0xb.9305cp-4, 0xb.ff911p-4,  
		 0xc.70165p-4, 0xc.e4bacp-4, 0xd.5da52p-4, 0xd.dafd7p-4,  
		 0xe.5ced3p-4, 0xe.e39f8p-4, 0xf.6f41p-4, 0x1p+0,  
		 0x1.0960c6p+0, 0x1.13198p+0, 0x1.1d2d64p+0, 0x1.279fcap+0,  
		 0x1.32742ap+0, 0x1.3dae1ep+0, 0x1.49515p+0, 0x1.5561aap+0,  
		 0x1.61e31ep+0, 0x1.6ed9eap+0, 0x1.7c4a4p+0, 0x1.8a38ap+0,  
		 0x1.98a9a4p+0, 0x1.a7a218p+0, 0x1.b726fp+0, 0x1.c73d52p+0,  
		 0x1.d7ea92p+0, 0x1.e93436p+0, 0x1.fb1ffcp+0};
	b16u16_u v = {.f = x};
	if ((v.u & 0x7c00) == 0x7c00) { // if x is nan or x is inf
		if (v.u == 0xfc00) return 0x0p0;
		else return x + x;
	}
	else if (v.u > x0) return (_Float16) 0x1p-25f;
	else if (x > x1) return (_Float16) 0x1.ffcp15f + 0x1p5f; 
	else {
		if (x == -0x1.c28p+0) return 0x1.1ccp-6 + 0x1p-18;
		if (x == 0x1p+2) return 0x1.388p+13;
		if (x == 0x1p+1) return 0x1.9p+6;
		if (x == 0x1.8p+1) return 0x1.f4p+9;
		if (x == 0x1p+0) return 0x1.4p3;
		if (x == -0x1.2ap-6) return 0x1.ebp-1 - 0x1p-13;
		if (x == -0x1.b5p-3) return 0x1.394p-1 - 0x1p-13;
		float log10_on_log2 = 0x1.a934fp1f;
		float minus_log2_on_log10 = -0x1.344136p-2f;
		float k = __builtin_roundevenf((float) x * log10_on_log2); 
		float xp = __builtin_fmaf(minus_log2_on_log10, k, (float) x);
		int i = 0x1p6 * xp;
		float xpp = __builtin_fmaf((float) i, -0x1p-6, xp); // x = klog(2)/log(10) + i/2^6 + xpp
																		   									// So, 10^x = 2^k * 10^(i/2^6) * 10^xpp

		// result
		xpp = __builtin_fmaf(__builtin_fmaf(__builtin_fmaf(__builtin_fmaf(0x1.2b9e52p0, xpp, 0x1.046efap1), xpp, 0x1.53524ep1), xpp, 0x1.26bb1cp1), xpp, 1.0);
		b32u32_u w = {.f = xpp * tb[i + 19]};
    w.u += (int) k * (1 << 23);
  	return w.f;
	}
}

// dummy function since GNU libc does not provide it
_Float16 exp10f16 (_Float16 x) {
  return (_Float16) exp10f ((float) x);
}
