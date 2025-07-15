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
#include <errno.h>
#include <math.h> // only used during performance tests

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;
static const b16u16_u neginf = {.u = 0xfc00};

/* This code is largely inspired by TANG, P. T. P. Table-driven
implementation of the logarithm function in IEEE floating-point
arithmetic. ACM Trans. Math. Softw. 16, 4 (1990), 378–400.
https://dl.acm.org/doi/10.1145/98267.98294 */

_Float16 cr_log2f16(_Float16 x){
	b32u32_u xf = {.f = x};
	if (xf.f == 0) {
#ifdef CORE_MATH_SUPPORT_ERRNO
		errno = ERANGE;
#endif
		return neginf.f;
	}
	else if (xf.u >> 23 >= 0xff) {
#ifdef CORE_MATH_SUPPORT_ERRNO
		errno = EDOM;
#endif
		if (xf.u >> 31) return 0.0f / 0.0f;
		else return x + x;
	}
	static const float tb[] = // tabulate value of log2(1 + i2^-5) for i in [0, 31]
		{0x0p+0f, 0xb.5d69cp-8f, 0x1.663f7p-4f, 0x2.118b1p-4f,
		 0x2.b80348p-4f, 0x3.59ebc4p-4f, 0x3.f782d8p-4f, 0x4.9101e8p-4f,
		 0x5.269e1p-4f, 0x5.b8887p-4f, 0x6.46eeap-4f, 0x6.d1fbp-4f,
		 0x7.59d4f8p-4f, 0x7.dea158p-4f, 0x8.60828p-4f, 0x8.df989p-4f,
		 0x9.5c01ap-4f, 0x9.d5dap-4f, 0xa.4d3c2p-4f, 0xa.c2411p-4f,
		 0xb.35004p-4f, 0xb.a58ffp-4f, 0xc.1404fp-4f, 0xc.80731p-4f,
		 0xc.eaedp-4f, 0xd.53848p-4f, 0xd.ba4a4p-4f, 0xe.1f4e5p-4f,
		 0xe.829fbp-4f, 0xe.e44cdp-4f, 0xf.44636p-4f, 0xf.a2f04p-4f};
	static const float tl[] = // tabulate value of 1 / (1 + i2^-5) for i in [0, 31]
		{0x1p-23f, 0xf.83e1p-27f, 0xf.0f0f1p-27f, 0xe.a0ea1p-27f,
		 0xe.38e39p-27f, 0xd.d67c9p-27f, 0xd.79436p-27f, 0xd.20d21p-27f,
		 0xc.ccccdp-27f, 0xc.7ce0cp-27f, 0xc.30c31p-27f, 0xb.e82fap-27f,
		 0xb.a2e8cp-27f, 0xb.60b61p-27f, 0xb.21643p-27f, 0xa.e4c41p-27f,
		 0xa.aaaabp-27f, 0xa.72f05p-27f, 0xa.3d70ap-27f, 0xa.0a0a1p-27f,
		 0x9.d89d9p-27f, 0x9.a90e8p-27f, 0x9.7b426p-27f, 0x9.4f209p-27f,
		 0x9.24925p-27f, 0x8.fb824p-27f, 0x8.d3dcbp-27f, 0x8.ad8f3p-27f,
		 0x8.88889p-27f, 0x8.64ba18p-27f, 0x8.42108p-27f, 0x8.20821p-27f};
	int expo = (xf.u >> 23) - 127; // used float instead of flaot16 to avoid working with subnormalized
	int i = (xf.u & 0x007c0000) >> 18;
	xf.f = (xf.u & 0x0003ffff) * tl[i];
	// We have, x = 2^expo * (1 + i2^-5 + xf.f)
	// Thus, log2(x) = expo + log2(1 + i2^-5) + log2(1 + xf.f / (1 + i2^-5))
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.ec709ep-2f, xf.f, -0x1.715476p-1f), xf.f, 0x1.715476p+0f);
	return (float) expo + tb[i] + xf.f;
}

// dummy function since GNU libc does not provide it
_Float16 log2f16 (_Float16 x) {
  return (_Float16) log2f ((float) x);
}
