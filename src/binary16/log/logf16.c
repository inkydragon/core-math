/* Correctly-rounded logarithm function for binary16 value.

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
static const b16u16_u neginf = {.u = 0xfc00};

/* This code is largely inspired by TANG, P. T. P. Table-driven
implementation of the logarithm function in IEEE floating-point
arithmetic. ACM Trans. Math. Softw. 16, 4 (1990), 378–400.
https://dl.acm.org/doi/10.1145/98267.98294 */

_Float16 cr_logf16(_Float16 x){
	b32u32_u xf = {.f = x};
	if (xf.u == 0) return neginf.f;
	else if (xf.u >> 23 >= 0xff) {
		if (xf.u == 0x80000000) return neginf.f;
		else if (xf.u >> 31) return 0.0f / 0.0f;
		else return x + x;
	}
	static const float log2 = 0x1.62e430p-1;
	static const float tb[] = // tabulate values of log(1 + i2^-5) for i in [0, 31]
		{0x0p+0f, 0x7.e0a6cp-8f, 0xf.85186p-8f, 0x1.6f0d28p-4f,  
		 0x1.e27076p-4f, 0x2.52aa6p-4f, 0x2.bfe61p-4f, 0x3.2a4b54p-4f,  
		 0x3.91fef8p-4f, 0x3.f7230cp-4f, 0x4.59d728p-4f, 0x4.ba38ccp-4f,  
		 0x5.1862fp-4f, 0x5.746f7p-4f, 0x5.ce76p-4f, 0x6.268cep-4f,  
		 0x6.7cc8f8p-4f, 0x6.d13dcp-4f, 0x7.23fdfp-4f, 0x7.751a8p-4f,  
		 0x7.c4a3d8p-4f, 0x8.12a95p-4f, 0x8.5f397p-4f, 0x8.aa61fp-4f,  
		 0x8.f42fbp-4f, 0x9.3caf1p-4f, 0x9.83ebap-4f, 0x9.c9f07p-4f,  
		 0xa.0ec7fp-4f, 0xa.527c3p-4f, 0xa.95169p-4f, 0xa.d6a02p-4f};
	static const float tl[] = // tabulate values of 1 / (1 + i2^-5) for i in [0, 31]
		{0x1p-23f, 0xf.83e1p-27f, 0xf.0f0f1p-27f, 0xe.a0ea1p-27f,  
		 0xe.38e39p-27f, 0xd.d67c9p-27f, 0xd.79436p-27f, 0xd.20d21p-27f,  
		 0xc.ccccdp-27f, 0xc.7ce0cp-27f, 0xc.30c31p-27f, 0xb.e82fap-27f,  
		 0xb.a2e8cp-27f, 0xb.60b61p-27f, 0xb.21643p-27f, 0xa.e4c41p-27f,  
		 0xa.aaaabp-27f, 0xa.72f05p-27f, 0xa.3d70ap-27f, 0xa.0a0a1p-27f,  
		 0x9.d89d9p-27f, 0x9.a90e8p-27f, 0x9.7b426p-27f, 0x9.4f209p-27f,  
		 0x9.24925p-27f, 0x8.fb824p-27f, 0x8.d3dcbp-27f, 0x8.ad8f3p-27f,  
		 0x8.88889p-27f, 0x8.64b8ap-27f, 0x8.42108p-27f, 0x8.208168p-27f};
	int expo = (xf.u >> 23) - 127; // used float instead of flaot16 to avoid working with subnormal numbers
	int i = (xf.u & 0x007c0000) >> 18;
	xf.f = (xf.u & 0x0003ffff) * tl[i];
	// We have, x = 2^expo * (1 + i2^-5 + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2^-5) + log(1 + xf.f / (1 + i2^-5))
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.555554p-2f, xf.f, -0.5f), xf.f, 1.0f);
	return __builtin_fmaf(log2, (float) expo, tb[i] + xf.f);
}

// dummy function since GNU libc does not provide it
_Float16 logf16 (_Float16 x) {
  return (_Float16) logf ((float) x);
}
