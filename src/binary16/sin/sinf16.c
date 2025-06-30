/* Correctly-rounded sine for binary16 value.

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

_Float16 cr_sinf16(_Float16 x){
	b16u16_u t = {.f = x};
	if ((t.u & 0x7c00) == 0x7c00) return 0.0f / 0.0f;
	if (!(t.u & 0x7fff)) return x;
	static const double thirtytwo_over_pi = 0x1.45f306dc9c883p+3;
	static const double minus_pi_over_thirtytwo = -0x1.921fb54442d18p-4;
	static const double tb_cos[] = // tabulate value of cos(i*pi/32) for i in [0, 63]
		{0x1p+0, 0xf.ec46d1e89293p-4, 0xf.b14be7fbae58p-4, 0xf.4fa0ab6316edp-4,  
		 0xe.c835e79946a3p-4, 0xe.1c5978c05ed88p-4, 0xd.4db3148750d18p-4, 0xc.5e40358a8ba08p-4,  
		 0xb.504f333f9de68p-4, 0xa.267992848eebp-4, 0x8.e39d9cd734648p-4, 0x7.8ad74e01bd8f8p-4,  
		 0x6.1f78a9abaa59p-4, 0x4.a5018bb567c14p-4, 0x3.1f17078d34c1ap-4, 0x1.917a6bc29b438p-4,  
		 0x4.69898cc51701cp-56, -0x1.917a6bc29b42fp-4, -0x3.1f17078d34c1p-4, -0x4.a5018bb567c08p-4,  
		 -0x6.1f78a9abaa588p-4, -0x7.8ad74e01bd8fp-4, -0x8.e39d9cd73463p-4, -0xa.267992848eea8p-4,  
		 -0xb.504f333f9de6p-4, -0xc.5e40358a8ba08p-4, -0xd.4db3148750d2p-4, -0xe.1c5978c05ed8p-4,  
		 -0xe.c835e79946a3p-4, -0xf.4fa0ab6316edp-4, -0xf.b14be7fbae58p-4, -0xf.ec46d1e892928p-4,  
		 -0x1p+0, -0xf.ec46d1e89293p-4, -0xf.b14be7fbae58p-4, -0xf.4fa0ab6316ed8p-4,  
		 -0xe.c835e79946a38p-4, -0xe.1c5978c05ed88p-4, -0xd.4db3148750d28p-4, -0xc.5e40358a8ba1p-4,  
		 -0xb.504f333f9de7p-4, -0xa.267992848eedp-4, -0x8.e39d9cd73464p-4, -0x7.8ad74e01bd8fcp-4,  
		 -0x6.1f78a9abaa5b4p-4, -0x4.a5018bb567c1cp-4, -0x3.1f17078d34c32p-4, -0x1.917a6bc29b421p-4,  
		 -0xd.3c9ca64f4505p-56, 0x1.917a6bc29b407p-4, 0x3.1f17078d34c18p-4, 0x4.a5018bb567cp-4,  
		 0x6.1f78a9abaa59cp-4, 0x7.8ad74e01bd8e8p-4, 0x8.e39d9cd734628p-4, 0xa.267992848eeb8p-4,  
		 0xb.504f333f9de58p-4, 0xc.5e40358a8b9fp-4, 0xd.4db3148750d18p-4, 0xe.1c5978c05ed78p-4,  
		 0xe.c835e79946a2p-4, 0xf.4fa0ab6316edp-4, 0xf.b14be7fbae578p-4, 0xf.ec46d1e89293p-4};
	static const double tb_sin[] = // tabulate value of sin(i*pi/32) for i in [0, 63]
		{0x0p+0, 0x1.917a6bc29b42cp-4, 0x3.1f17078d34c14p-4, 0x4.a5018bb567c14p-4,  
		 0x6.1f78a9abaa58cp-4, 0x7.8ad74e01bd8ecp-4, 0x8.e39d9cd73464p-4, 0xa.267992848eebp-4,  
		 0xb.504f333f9de6p-4, 0xc.5e40358a8ba08p-4, 0xd.4db3148750d18p-4, 0xe.1c5978c05ed8p-4,  
		 0xe.c835e79946a3p-4, 0xf.4fa0ab6316ed8p-4, 0xf.b14be7fbae58p-4, 0xf.ec46d1e892928p-4,  
		 0x1p+0, 0xf.ec46d1e89293p-4, 0xf.b14be7fbae58p-4, 0xf.4fa0ab6316ed8p-4,  
		 0xe.c835e79946a3p-4, 0xe.1c5978c05ed88p-4, 0xd.4db3148750d28p-4, 0xc.5e40358a8ba1p-4,  
		 0xb.504f333f9de68p-4, 0xa.267992848eebp-4, 0x8.e39d9cd73464p-4, 0x7.8ad74e01bd8fcp-4,  
		 0x6.1f78a9abaa594p-4, 0x4.a5018bb567c18p-4, 0x3.1f17078d34c2ep-4, 0x1.917a6bc29b43cp-4,  
		 0x8.d313198a2e038p-56, -0x1.917a6bc29b42bp-4, -0x3.1f17078d34c1cp-4, -0x4.a5018bb567c04p-4,  
		 -0x6.1f78a9abaa584p-4, -0x7.8ad74e01bd8ecp-4, -0x8.e39d9cd73463p-4, -0xa.267992848eeap-4,  
		 -0xb.504f333f9de6p-4, -0xc.5e40358a8b9fp-4, -0xd.4db3148750d18p-4, -0xe.1c5978c05ed8p-4,  
		 -0xe.c835e79946a2p-4, -0xf.4fa0ab6316edp-4, -0xf.b14be7fbae578p-4, -0xf.ec46d1e89293p-4,  
		 -0x1p+0, -0xf.ec46d1e89293p-4, -0xf.b14be7fbae58p-4, -0xf.4fa0ab6316ed8p-4,  
		 -0xe.c835e79946a28p-4, -0xe.1c5978c05ed88p-4, -0xd.4db3148750d28p-4, -0xc.5e40358a8bap-4,  
		 -0xb.504f333f9de7p-4, -0xa.267992848eedp-4, -0x8.e39d9cd73464p-4, -0x7.8ad74e01bd9p-4,  
		 -0x6.1f78a9abaa5b8p-4, -0x4.a5018bb567c2p-4, -0x3.1f17078d34c36p-4, -0x1.917a6bc29b425p-4};
	double xd = x;
	double j = __builtin_roundeven(thirtytwo_over_pi * xd);
	int i = (uint64_t) j & 0x3f;
	double xp = __builtin_fma(minus_pi_over_thirtytwo, j, xd);
	// xd = j*pi/32 + xp = 2kpi + i*pi/32 + xp with k an integer
	// so sin(xd) = sin(i*pi/32 + xp) = sin(i*pi/32)cos(xp) + cos(i*pi/32)sin(xp)
	double xp2 = xp*xp;
	return tb_sin[i] * (1.0 + xp2 * (-0.5 + xp2 * 0x1.5555p-5)) + tb_cos[i] * (xp - 0x1.5555p-3 * xp * xp2);
}

// dummy function since GNU libc does not provide it
_Float16 sinf16 (_Float16 x) {
  return (_Float16) sinf ((float) x);
}
