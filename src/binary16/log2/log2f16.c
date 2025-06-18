/* Correctly-rounded binary logarithm function for binary16 value.

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
b16u16_u neginf = {.u = 0xfc00};

_Float16 cr_log2f16(_Float16 x){
	b16u16_u t = {.f = x};
	if (t.u == 0) return neginf.f;
	else if (t.u >> 10 >= 0x1f) {
		if (t.u == 0x8000) return neginf.f;
		else if (t.u >> 15) return 0.0 / 0.0;
		else return x + x;
	}
	static const float tb[] =
		{0x0p+0, 0xb.5d69cp-8, 0x1.663f7p-4, 0x2.118b1p-4,
		 0x2.b80348p-4, 0x3.59ebc4p-4, 0x3.f782d8p-4, 0x4.9101e8p-4,
		 0x5.269e1p-4, 0x5.b8887p-4, 0x6.46eeap-4, 0x6.d1fbp-4,
		 0x7.59d4f8p-4, 0x7.dea158p-4, 0x8.60828p-4, 0x8.df989p-4,
		 0x9.5c01ap-4, 0x9.d5dap-4, 0xa.4d3c2p-4, 0xa.c2411p-4,
		 0xb.35004p-4, 0xb.a58ffp-4, 0xc.1404fp-4, 0xc.80731p-4,
		 0xc.eaedp-4, 0xd.53848p-4, 0xd.ba4a4p-4, 0xe.1f4e5p-4,
		 0xe.829fbp-4, 0xe.e44cdp-4, 0xf.44636p-4, 0xf.a2f04p-4};
	static const float tl[] =
		{0x1p+0, 0xf.83e1p-4, 0xf.0f0f1p-4, 0xe.a0ea1p-4,
		 0xe.38e39p-4, 0xd.d67c9p-4, 0xd.79436p-4, 0xd.20d21p-4,
		 0xc.ccccdp-4, 0xc.7ce0cp-4, 0xc.30c31p-4, 0xb.e82fap-4,
		 0xb.a2e8cp-4, 0xb.60b61p-4, 0xb.21643p-4, 0xa.e4c41p-4,
		 0xa.aaaabp-4, 0xa.72f05p-4, 0xa.3d70ap-4, 0xa.0a0a1p-4,
		 0x9.d89d9p-4, 0x9.a90e8p-4, 0x9.7b426p-4, 0x9.4f209p-4,
		 0x9.24925p-4, 0x8.fb824p-4, 0x8.d3dcbp-4, 0x8.ad8f3p-4,
		 0x8.88889p-4, 0x8.64ba18p-4, 0x8.42108p-4, 0x8.20821p-4};
	b32u32_u xf = {.f = x};
	int expo = (xf.u >> 23) - 127; // used float instead of flaot16 to avoid working with subnormalized
	int i = (xf.u & 0x007c0000) >> 18;
	xf.f = 0x1p-23 * (xf.u & 0x0003ffff) * tl[i];
	// We have, x = 2^expo * (1 + i2^-5 + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2^-5) + log(1 + xf.f / (1 + i2^-5))
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.555554p-2, xf.f, -0.5f), xf.f, 1.0f) * 0x1.715476p+0;
	return (float) expo + tb[i] + xf.f;
}

// dummy function since GNU libc does not provide it
_Float16 log2f16 (_Float16 x) {
  return (_Float16) log2f ((float) x);
}
