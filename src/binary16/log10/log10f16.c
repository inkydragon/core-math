/* Correctly-rounded radix-10 logarithm function for binary16 value.

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

_Float16 cr_log10f16(_Float16 x){
	b16u16_u t = {.f = x};
	if (t.u == 0) return neginf.f;
	else if (t.u >> 10 >= 0x1f) {
		if (t.u == 0x8000) return neginf.f;
		else if (t.u >> 15) return 0.0f / 0.0f;
		else return x + x;
	}
	static const float log10_2 = 0x1.344136p-2;
	static const float tb[] = // tabulate values of log10(1 + i2^-5) for i in [0, 31]
		{0x0p+0f, 0x3.6bd21p-8f, 0x6.bd7e48p-8f, 0x9.f688dp-8f,  
		 0xd.1854fp-8f, 0x1.02428cp-4f, 0x1.31b306p-4f, 0x1.5fe804p-4f,  
		 0x1.8cf178p-4f, 0x1.b8de4ep-4f, 0x1.e3bc1ap-4f, 0x2.0d97acp-4f,  
		 0x2.367ce4p-4f, 0x2.5e76d4p-4f, 0x2.858fdp-4f, 0x2.abd18cp-4f,  
		 0x2.d1451p-4f, 0x2.f5f2e8p-4f, 0x3.19e2f2p-4f, 0x3.3d1cf8p-4f,  
		 0x3.5fa7c8p-4f, 0x3.818a28p-4f, 0x3.a2ca6p-4f, 0x3.c36e68p-4f,  
		 0x3.e37bdcp-4f, 0x4.02f818p-4f, 0x4.21e82p-4f, 0x4.4050c8p-4f,  
		 0x4.5e3698p-4f, 0x4.7b9dep-4f, 0x4.988acp-4f, 0x4.b5013p-4f};
	static const float tl[] = // tabulate values of 1 / (1 + i2^-5) for i in [0, 31]
		{0x1p-23f, 0xf.83e1p-27f, 0xf.0f0f1p-27f, 0xe.a0ea1p-27f,  
		 0xe.38e39p-27f, 0xd.d67c9p-27f, 0xd.79436p-27f, 0xd.20cdep-27f,  
		 0xc.ccccdp-27f, 0xc.7ce0cp-27f, 0xc.30c31p-27f, 0xb.e82fap-27f,  
		 0xb.a2cbc8p-27f, 0xb.60b61p-27f, 0xb.21643p-27f, 0xa.e4c41p-27f,  
		 0xa.aaaabp-27f, 0xa.72ed28p-27f, 0xa.3d70ap-27f, 0xa.0a0a1p-27f,  
		 0x9.d89d9p-27f, 0x9.a90e8p-27f, 0x9.7b426p-27f, 0x9.4f209p-27f,  
		 0x9.24925p-27f, 0x8.fb824p-27f, 0x8.d3dcbp-27f, 0x8.ad8f3p-27f,  
		 0x8.88938ap-27f, 0x8.64bfd8p-27f, 0x8.4218b8p-27f, 0x8.208136p-27f};

        // deal with some exceptional cases
        if (t.u == 0x38e5) return -0x1.b5p-3f + 0x1p-15f; // x=1253/2048
        if (t.u == 0x70e2) return 0x1p+2f; // x=10000
        if (t.u == 0x57e1) return 0x1.0ccp+1f + 0x1p-11f; // x=2017/16
        if (t.u == 0x63d0) return 0x1.8p+1f; // x=1000
	b32u32_u xf = {.f = x};
	int expo = (xf.u >> 23) - 127; // used float instead of float16 to avoid working with subnormal numbers
	uint32_t i = (xf.u & 0x007c0000) >> 18;
	xf.f = (xf.u & 0x0003ffff) * tl[i];
	// We have, x = 2^expo * (1 + i2^-5 + xf.f)
	// Thus, log10(x) = expo log10(2) + log10(1 + i2^-5) + log10(1 + xf.f / (1 + i2^-5))
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.287a76p-3f, xf.f, -0x1.bcb7b2p-3f), xf.f, 0x1.bcb7b2p-2f);
	return __builtin_fmaf(log10_2, (float) expo, xf.f + tb[i]); 
}

// dummy function since GNU libc does not provide it
_Float16 log10f16 (_Float16 x) {
  return (_Float16) log10f ((float) x);
}
