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
#include <math.h> // only used during performance tests

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_expf16(_Float16 x){
	uint16_t x0 = 0xcc55; // binary representation of x0 in order to compare uint16_t rather than flt16
 	_Float16 x1 = 0x1.62cp3; // largest _Float16 such that exp(x1) <= MAX_FLOAT16 < exp(x1+)
	static const float tb[] = // tabulate value of exp(i/2^6) for i in [-2^6*log(2)/2, 2^6*log(2)/2], size(tb) = 45
		{0xb.587fcp-4, 0xb.863cfp-4, 0xb.b4b29p-4, 0xb.e3e38p-4,
		 0xc.13d2bp-4, 0xc.44832p-4, 0xc.75f7dp-4, 0xc.a8340p-4,
		 0xc.db3a8p-4, 0xd.0f0edp-4, 0xd.43b41p-4, 0xd.792d8p-4,
		 0xd.af7e9p-4, 0xd.e6aaap-4, 0xe.1eb51p-4, 0xe.57a17p-4,
		 0xe.91735p-4, 0xe.cc2e4p-4, 0xf.07d6p-4, 0xf.446e3p-4,
		 0xf.81fabp-4, 0xf.c07f5p-4, 0x1p0, 0x1.04080ap+0,
		 0x1.082056p+0, 0x1.0c4924p+0, 0x1.1082b6p+0, 0x1.14cd5p+0,
		 0x1.192938p+0, 0x1.1d96bp+0, 0x1.221604p+0, 0x1.26a77ap+0,
		 0x1.2b4b58p+0, 0x1.3001ecp+0, 0x1.34cb82p+0, 0x1.39a862p+0,
		 0x1.3e98dep+0, 0x1.439d44p+0, 0x1.48b5e4p+0, 0x1.4de30ep+0,
		 0x1.532518p+0, 0x1.587c54p+0, 0x1.5de918p+0, 0x1.636bbap+0,
		 0x1.690492p+0};

	b16u16_u v = {.f = x};
	if ((v.u & 0x7c00) == 0x7c00) { // if x is nan or x is inf
		if (v.u == 0xfc00) return 0x0p0;
		else return x + x;
	}
	else if (v.u > x0) return (_Float16) 0x1p-25f;
	else if (x > x1) return (_Float16) 0x1.ffcp15f + 0x1p5f; 
	else {
                float xf = x; // exact conversion from _Float16 to float
		float minus_log2 = -0x1.62e430p-1;
		float inv_log2 = 0x1.715476p0f;
		float k = __builtin_roundevenf(inv_log2 * xf);
		float xp = __builtin_fmaf(k, minus_log2, xf); // xp is a float such that |xp| is minimal and x = klog(2) + xp
		int i = 0x1p6 * xp;
		if ((uint16_t) (i & 0x80000001) <= 1) { // some wrong cases
			if (xf == 0x1.de4p-8) return (0x1.01cp+0 + 0x1p-12);
			if (xf == 0x1.73cp-6) return (0x1.05cp+0 + 0x1p-12);
			if (xf == 0x1.62cp+3) return (0x1.fdcp+15 - 1);
		}
		float xpp = __builtin_fmaf((float) i , -0x1p-6f, xp); // x = klog(2) + i/2^6 + xpp
																													// So, exp(x) = 2^k * exp(i/2^6) * exp(xpp)
		// result
		xpp = __builtin_fmaf(__builtin_fmaf(__builtin_fmaf(xpp, 0x1.555644p-3, 0.5), xpp, 1.0), xpp, 1.0);
		b32u32_u w = {.f = xpp * tb[i + 22]};
		w.u += (int) k << 23;
		return w.f; // conversion float -> _Float16 (with rounding)
	}
}

// dummy function since GNU libc does not provide it
_Float16 expf16 (_Float16 x) {
  return (_Float16) expf ((float) x);
}
