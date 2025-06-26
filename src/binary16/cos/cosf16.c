/* Correctly-rounded cosine for binary16 value.

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
typedef union {double f; uint64_t u;} b64u64_u;

_Float16 cr_cosf16(_Float16 x){
	b16u16_u t = {.f = x};
	if ((t.u & 0x7c00) == 0x7c00) return 0.0f / 0.0f;
	static const double sixteen_over_pi = 0x1.45f306dc9c883p+2;
	static const double minus_pi_over_sixteen = -0x1.921fb54442d18p-3;
	static const double tb_cos[] = 
		{0x1p+0, 0xf.b14be7fbae58p-4, 0xe.c835e79946a3p-4, 0xd.4db3148750d18p-4,  
		 0xb.504f333f9de68p-4, 0x8.e39d9cd734648p-4, 0x6.1f78a9abaa59p-4, 0x3.1f17078d34c1ap-4,  
		 0x4.69898cc51701cp-56, -0x3.1f17078d34c1p-4, -0x6.1f78a9abaa588p-4, -0x8.e39d9cd73463p-4,  
		 -0xb.504f333f9de6p-4, -0xd.4db3148750d2p-4, -0xe.c835e79946a3p-4, -0xf.b14be7fbae58p-4,  
		 -0x1p+0, -0xf.b14be7fbae58p-4, -0xe.c835e79946a38p-4, -0xd.4db3148750d28p-4,  
		 -0xb.504f333f9de7p-4, -0x8.e39d9cd73464p-4, -0x6.1f78a9abaa5b4p-4, -0x3.1f17078d34c32p-4,  
		 -0xd.3c9ca64f4505p-56, 0x3.1f17078d34c18p-4, 0x6.1f78a9abaa59cp-4, 0x8.e39d9cd734628p-4,  
		 0xb.504f333f9de58p-4, 0xd.4db3148750d18p-4, 0xe.c835e79946a2p-4, 0xf.b14be7fbae578p-4};
	static const double tb_sin[] = 
		{0x0p+0, 0x3.1f17078d34c14p-4, 0x6.1f78a9abaa58cp-4, 0x8.e39d9cd73464p-4,  
		 0xb.504f333f9de6p-4, 0xd.4db3148750d18p-4, 0xe.c835e79946a3p-4, 0xf.b14be7fbae58p-4,  
		 0x1p+0, 0xf.b14be7fbae58p-4, 0xe.c835e79946a3p-4, 0xd.4db3148750d28p-4,  
		 0xb.504f333f9de68p-4, 0x8.e39d9cd73464p-4, 0x6.1f78a9abaa594p-4, 0x3.1f17078d34c2ep-4,  
		 0x8.d313198a2e038p-56, -0x3.1f17078d34c1cp-4, -0x6.1f78a9abaa584p-4, -0x8.e39d9cd73463p-4,  
		 -0xb.504f333f9de6p-4, -0xd.4db3148750d18p-4, -0xe.c835e79946a2p-4, -0xf.b14be7fbae578p-4,  
		 -0x1p+0, -0xf.b14be7fbae58p-4, -0xe.c835e79946a28p-4, -0xd.4db3148750d28p-4,  
		 -0xb.504f333f9de7p-4, -0x8.e39d9cd73464p-4, -0x6.1f78a9abaa5b8p-4, -0x3.1f17078d34c36p-4};
	double xd = x;
	double j = __builtin_roundeven(sixteen_over_pi * xd);
	int i = (uint64_t) j & 0x1f;
	double xp = __builtin_fma(minus_pi_over_sixteen, j, xd);
	double xp2 = xp*xp;
	return tb_cos[i] * (1.0-0.5*xp2+0x1.555556p-5*xp2*xp2) - tb_sin[i] * (xp - 0x1.555556p-3*xp*xp2);
}

// dummy function since GNU libc does not provide it
_Float16 cosf16 (_Float16 x) {
  return (_Float16) sinf ((float) x);
}
