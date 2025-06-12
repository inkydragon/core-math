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

<<<<<<< HEAD

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;
static const b16u16_u neginf = {.u = 0xfc00};

/* This code is largely inspired by TANG, P. T. P. Table-driven
implementation of the logarithm function in IEEE floating-point
arithmetic. ACM Trans. Math. Softw. 16, 4 (1990), 378–400.
https://dl.acm.org/doi/10.1145/98267.98294 */

_Float16 cr_logf16(_Float16 x){
	b16u16_u t = {.f = x};
	if (t.u == 0) return neginf.f;
	else if (t.u >> 10 >= 0x1f) {
		if (t.u == 0x8000) return neginf.f;
=======
typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_logf16(_Float16 x){
	b16u16_u t = {.f = x};
	if (t.u == 0) return -0xffcp15f - 0x1p5f;
	else if (t.u >> 10 >= 0x1f) {
		if (t.u == 0x8000) return -0xffcp+15f - 0x1p5f;
>>>>>>> 1fe06f90 (First implementation of logf16 with 4000 wrong cases)
		else if (t.u >> 15) return 0.0 / 0.0;
		else return x + x;
	}
	float log2 = 0x1.62e430p-1;
<<<<<<< HEAD
	static const float tb[] = // tabulate values of log(1 + i2^-5)
		{0x0p+0, 0x7.e0a6cp-8, 0xf.85186p-8, 0x1.6f0d28p-4,  
		 0x1.e27076p-4, 0x2.52aa6p-4, 0x2.bfe61p-4, 0x3.2a4b54p-4,  
		 0x3.91fef8p-4, 0x3.f7230cp-4, 0x4.59d728p-4, 0x4.ba38ccp-4,  
		 0x5.1862fp-4, 0x5.746f7p-4, 0x5.ce76p-4, 0x6.268cep-4,  
		 0x6.7cc8f8p-4, 0x6.d13dcp-4, 0x7.23fdfp-4, 0x7.751a8p-4,  
		 0x7.c4a3d8p-4, 0x8.12a95p-4, 0x8.5f397p-4, 0x8.aa61fp-4,  
		 0x8.f42fbp-4, 0x9.3caf1p-4, 0x9.83ebap-4, 0x9.c9f07p-4,  
		 0xa.0ec7fp-4, 0xa.527c3p-4, 0xa.95169p-4, 0xa.d6a02p-4};
	static const float tl[] = // tabulate values of 1 / (1 + i2^-5)
		{0x1p+0, 0xf.83e1p-4, 0xf.0f0f1p-4, 0xe.a0ea1p-4,  
		 0xe.38e39p-4, 0xd.d67c9p-4, 0xd.79436p-4, 0xd.20d21p-4,  
		 0xc.ccccdp-4, 0xc.7ce0cp-4, 0xc.30c31p-4, 0xb.e82fap-4,  
		 0xb.a2e8cp-4, 0xb.60b61p-4, 0xb.21643p-4, 0xa.e4c41p-4,  
		 0xa.aaaabp-4, 0xa.72f05p-4, 0xa.3d70ap-4, 0xa.0a0a1p-4,  
		 0x9.d89d9p-4, 0x9.a90e8p-4, 0x9.7b426p-4, 0x9.4f209p-4,  
		 0x9.24925p-4, 0x8.fb824p-4, 0x8.d3dcbp-4, 0x8.ad8f3p-4,  
		 0x8.88889p-4, 0x8.64b8ap-4, 0x8.42108p-4, 0x8.208168p-4};
	b32u32_u xf = {.f = x};
	int expo = (xf.u >> 23) - 127; // used float instead of flaot16 to avoid working with subnormal numbers
	int i = (xf.u & 0x007c0000) >> 18;
	xf.f = 0x1p-23 * (xf.u & 0x0003ffff) * tl[i];
	// We have, x = 2^expo * (1 + i2^-5 + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2^-5) + log(1 + xf.f / (1 + i2^-5))
=======
	static const float tb[] = 
		{0x0p+0, 0x3.f81518p-8, 0x7.e0a6cp-8, 0xb.ba2c8p-8,  
		 0xf.85186p-8, 0x1.341d7ap-4, 0x1.6f0d28p-4, 0x1.a926d4p-4,  
		 0x1.e27076p-4, 0x2.1aefdp-4, 0x2.52aa6p-4, 0x2.89a56cp-4,  
		 0x2.bfe61p-4, 0x2.f5712p-4, 0x3.2a4b54p-4, 0x3.5e7928p-4,  
		 0x3.91fef8p-4, 0x3.c4e0ecp-4, 0x3.f7230cp-4, 0x4.28c938p-4,  
		 0x4.59d728p-4, 0x4.8a508p-4, 0x4.ba38bp-4, 0x4.e99318p-4,  
		 0x5.1862fp-4, 0x5.46ab6p-4, 0x5.746f7p-4, 0x5.a1b208p-4,  
		 0x5.ce76p-4, 0x5.fabe1p-4, 0x6.268cep-4, 0x6.51e508p-4,  
		 0x6.7cc8f8p-4, 0x6.a73b28p-4, 0x6.d13dep-4, 0x6.fad368p-4,  
		 0x7.23fdfp-4, 0x7.4cbfap-4, 0x7.751a8p-4, 0x7.9d1098p-4,  
		 0x7.c4a3d8p-4, 0x7.ebd62p-4, 0x8.12a95p-4, 0x8.391f3p-4,  
		 0x8.5f397p-4, 0x8.84f9dp-4, 0x8.aa61fp-4, 0x8.cf736p-4,  
		 0x8.f42fbp-4, 0x9.18987p-4, 0x9.3caf1p-4, 0x9.6074fp-4,  
		 0x9.83ebap-4, 0x9.a7145p-4, 0x9.c9f07p-4, 0x9.ec813p-4,  
		 0xa.0ec7fp-4, 0xa.30c5ep-4, 0xa.527c3p-4, 0xa.73ec1p-4,  
		 0xa.95169p-4, 0xa.b5fcfp-4, 0xa.d6a02p-4, 0xa.f7015p-4};
	b32u32_u xf = {.f = x};
	int expo = (xf.u >> 23) - 127; // used float instead of flaot16 to avoid working with subnormalized
	int i = (t.u & 0x03f0) >> 4;
	xf.f = 0x1p-10 * (t.u & 0x000f);
	// We have, x = 2^expo * (1 + i2⁻⁶ + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2⁻⁶) + log(1 + xf.f / (1 + i2⁻⁶))
	xf.f /= __builtin_fmaf(0x1p-6, (float) i, (float) 1);
>>>>>>> 1fe06f90 (First implementation of logf16 with 4000 wrong cases)
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.555554p-2, xf.f, -0.5f), xf.f, 1.0f);
	return __builtin_fmaf(log2, (float) expo, tb[i] + xf.f);
}

// dummy function since GNU libc does not provide it
_Float16 logf16 (_Float16 x) {
  return (_Float16) logf ((float) x);
}
