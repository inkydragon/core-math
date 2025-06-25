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
	if ((t.u & 0x7fff) < 0x2800) return (float) x * (1.0f - 0x1.555556p-3f * x * x);
	static const double inv_2pi = 0x1.45f306dc9c883p-3;
	static const double minus_2pi = -0x1.921fb54442d18p+2;
	static const float tb_sin[] = // tabulate values of sin(i2^-4) for i in [-2^4pi, 2^4pi]
		{-0x4.3f5dcp-8f, -0x1.43a038p-4f, -0x2.42070cp-4f, -0x3.3e2c0cp-4f,  
		 -0x4.371328p-4f, -0x5.2bc388p-4f, -0x6.1b489p-4f, -0x7.04b2dp-4f,  
		 -0x7.e718f8p-4f, -0x8.c198bp-4f, -0x9.93578p-4f, -0xa.5b83dp-4f,  
		 -0xb.19558p-4f, -0xb.cc0ebp-4f, -0xc.72fcep-4f, -0xd.0d792p-4f,  
		 -0xd.9ae8fp-4f, -0xe.1abfp-4f, -0xe.8c7b7p-4f, -0xe.efacbp-4f,  
		 -0xf.43ef7p-4f, -0xf.88ef9p-4f, -0xf.be681p-4f, -0xf.e4237p-4f,  
		 -0xf.f9fcp-4f, -0xf.ffdbfp-4f, -0xf.f5bd5p-4f, -0xf.dbaa5p-4f,  
		 -0xf.b1bcfp-4f, -0xf.781f2p-4f, -0xf.2f0a8p-4f, -0xe.d6c81p-4f,  
		 -0xe.6fb02p-4f, -0xd.fa29bp-4f, -0xd.76aa4p-4f, -0xc.e5b55p-4f,  
		 -0xc.47dbbp-4f, -0xb.9dbb4p-4f, -0xa.e7fe1p-4f, -0xa.2759cp-4f,  
		 -0x9.5c8efp-4f, -0x8.88686p-4f, -0x7.abba2p-4f, -0x6.c760cp-4f,  
		 -0x5.dc4098p-4f, -0x4.eb44a8p-4f, -0x3.f55ddcp-4f, -0x2.fb8204p-4f,  
		 -0x1.feaaeep-4f, -0xf.fd557p-8f, 0x0p+0f, 0xf.fd557p-8f,  
		 0x1.feaaeep-4f, 0x2.fb8204p-4f, 0x3.f55ddcp-4f, 0x4.eb44a8p-4f,  
		 0x5.dc4098p-4f, 0x6.c760cp-4f, 0x7.abba2p-4f, 0x8.88686p-4f,  
		 0x9.5c8efp-4f, 0xa.2759cp-4f, 0xa.e7fe1p-4f, 0xb.9dbb4p-4f,  
		 0xc.47dbbp-4f, 0xc.e5b55p-4f, 0xd.76aa4p-4f, 0xd.fa29bp-4f,  
		 0xe.6fb02p-4f, 0xe.d6c81p-4f, 0xf.2f0a8p-4f, 0xf.781f2p-4f,  
		 0xf.b1bcfp-4f, 0xf.dbaa5p-4f, 0xf.f5bd5p-4f, 0xf.ffdbfp-4f,  
		 0xf.f9fcp-4f, 0xf.e4237p-4f, 0xf.be681p-4f, 0xf.88ef9p-4f,  
		 0xf.43ef7p-4f, 0xe.efacbp-4f, 0xe.8c7b7p-4f, 0xe.1abfp-4f,  
		 0xd.9ae8fp-4f, 0xd.0d792p-4f, 0xc.72fcep-4f, 0xb.cc0ebp-4f,  
		 0xb.19558p-4f, 0xa.5b83dp-4f, 0x9.93578p-4f, 0x8.c198bp-4f,  
		 0x7.e718f8p-4f, 0x7.04b2dp-4f, 0x6.1b489p-4f, 0x5.2bc388p-4f,  
		 0x4.371328p-4f, 0x3.3e2c0cp-4f, 0x2.42070cp-4f, 0x1.43a038p-4f,  
		 0x4.3f5dcp-8f};
	static const float tb_cos[] = // tabulate values of cos(i2^-4) for i in [-2^4pi, 2^4pi]
		{-0xf.ff6fbp-4f, -0xf.f331fp-4f, -0xf.d7026p-4f, -0xf.aafd2p-4f,  
		 -0xf.6f4e3p-4f, -0xf.24313p-4f, -0xe.c9f15p-4f, -0xe.60e8bp-4f,  
		 -0xd.e9806p-4f, -0xd.642fbp-4f, -0xc.d17bfp-4f, -0xc.31f7dp-4f,  
		 -0xb.8642bp-4f, -0xa.cf085p-4f, -0xa.0cffdp-4f, -0x9.40eb1p-4f,  
		 -0x8.6b964p-4f, -0x7.8dd688p-4f, -0x6.a88998p-4f, -0x5.bc94a8p-4f,  
		 -0x4.cae3a8p-4f, -0x3.d46824p-4f, -0x2.da1888p-4f, -0x1.dcef14p-4f,  
		 -0xd.de8d8p-8f, 0x2.1fb3acp-8f, 0x1.21bd54p-4f, 0x2.205dccp-4f,  
		 0x3.1cde1p-4f, 0x4.1641b8p-4f, 0x5.0b8f78p-4f, 0x5.fbd21p-4f,  
		 0x6.e61958p-4f, 0x7.c97b18p-4f, 0x8.a514p-4f, 0x9.7808ap-4f,  
		 0xa.4185fp-4f, 0xb.00c29p-4f, 0xb.b4ff6p-4f, 0xc.5d883p-4f,  
		 0xc.f9b47p-4f, 0xd.88e82p-4f, 0xe.0a94p-4f, 0xe.7e368p-4f,  
		 0xe.e35bfp-4f, 0xf.399f5p-4f, 0xf.80aa5p-4f, 0xf.b835fp-4f,  
		 0xf.e00abp-4f, 0xf.f800bp-4f, 0x1p+0f, 0xf.f800bp-4f,  
		 0xf.e00abp-4f, 0xf.b835fp-4f, 0xf.80aa5p-4f, 0xf.399f5p-4f,  
		 0xe.e35bfp-4f, 0xe.7e368p-4f, 0xe.0a94p-4f, 0xd.88e82p-4f,  
		 0xc.f9b47p-4f, 0xc.5d883p-4f, 0xb.b4ff6p-4f, 0xb.00c29p-4f,  
		 0xa.4185fp-4f, 0x9.7808ap-4f, 0x8.a514p-4f, 0x7.c97b18p-4f,  
		 0x6.e61958p-4f, 0x5.fbd21p-4f, 0x5.0b8f78p-4f, 0x4.1641b8p-4f,  
		 0x3.1cde1p-4f, 0x2.205dccp-4f, 0x1.21bd54p-4f, 0x2.1fb3acp-8f,  
		 -0xd.de8d8p-8f, -0x1.dcef14p-4f, -0x2.da1888p-4f, -0x3.d46824p-4f,  
		 -0x4.cae3a8p-4f, -0x5.bc94a8p-4f, -0x6.a88998p-4f, -0x7.8dd688p-4f,  
		 -0x8.6b964p-4f, -0x9.40eb1p-4f, -0xa.0cffdp-4f, -0xa.cf085p-4f,  
		 -0xb.8642bp-4f, -0xc.31f7dp-4f, -0xc.d17bfp-4f, -0xd.642fbp-4f,  
		 -0xd.e9806p-4f, -0xe.60e8bp-4f, -0xe.c9f15p-4f, -0xf.24313p-4f,  
		 -0xf.6f4e3p-4f, -0xf.aafd2p-4f, -0xf.d7026p-4f, -0xf.f331fp-4f,  
		 -0xf.ff6fbp-4f};
	double xd = x;
	double k = __builtin_roundeven(xd * inv_2pi);
	float xp = __builtin_fma(minus_2pi, k, xd);
	// x = 2kpi + xp
	float i = __builtin_roundevenf(0x1p4f * xp);
	float xpp = __builtin_fmaf(i, -0x1p-4f, xp);
	float xpp2 = xpp*xpp;
	// x = 2kpi + i2^-4 + xpp
	// sin(x) = sin(i2^-4 + xpp) = sin(i2^-4)cos(xpp) + sin(xpp)cos(i2^-4)
	return __builtin_fmaf(tb_sin[(int)i+50], __builtin_fmaf(-0.5f, xpp2, 1.0f), __builtin_fmaf(-0x1.555554p-3f, xpp*xpp2, xpp) * tb_cos[(int)i+50]);
}

// dummy function since GNU libc does not provide it
_Float16 sinf16 (_Float16 x) {
  return (_Float16) sinf ((float) x);
}
