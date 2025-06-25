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

// convert x to double, assuming x is not NaN, Inf or 0
static double to_double (_Float16 x) {
  b16u16_u v = {.f = x};
  b64u64_u w;
  uint64_t u = v.u, au = u & 0x7fff;
  if (au < 0x400) { // subnormal case
    uint64_t nbZero = __builtin_clz(au); // number of leading 0 on 32 bits
    w.u = (au << (21 + nbZero)) + (0x3f00000000000000 - ((nbZero - 21) << 52));
  } else
    // 1 is encoded by 0x3c00 in _Float16, and by 0x3ff0 0000 0000 0000 in double 
    w.u = (au << 42) + 0x3f00000000000000;
  w.u |= (u >> 15) << 63; // sign
  return w.f;
}

_Float16 cr_sinf16(_Float16 x){
	b16u16_u t = {.f = x};
	if ((t.u & 0x7c00) == 0x7c00) return 0.0f / 0.0f;
	if (__builtin_expect(!(t.u & 0x7fff), 0)) return x;
	static const double inv_2pi = 0x1.45f306dc9c883p-3;
	static const double minus_2pi = -0x1.921fb54442d18p+2;
	static const double tb_sin[] = // tabulate values of sin(i2^-4) for i in [-2^4pi, 2^4pi]
		{-0x4.3f5dc280f96a8p-8, -0x1.43a0378fadb65p-4, -0x2.42070db6daab6p-4, -0x3.3e2c0cf9f6e7p-4,  
		 -0x4.3713255c21484p-4, -0x5.2bc38482335b4p-4, -0x6.1b488e7060fecp-4, -0x7.04b2d211d289p-4,  
		 -0x7.e718f895c631cp-4, -0x8.c198aeb2ef95p-4, -0x9.935786e7e5588p-4, -0xa.5b83d3d66f668p-4,  
		 -0xb.195579df6c9d8p-4, -0xb.cc0eb72dc5888p-4, -0xc.72fce16871d48p-4, -0xd.0d79184dee53p-4,  
		 -0xd.9ae8ec8696ebp-4, -0xe.1abefa072012p-4, -0xe.8c7b7568da23p-4, -0xe.efacabaa7222p-4,  
		 -0xf.43ef73d979c98p-4, -0xf.88ef923316d68p-4, -0xf.be680c58c123p-4, -0xf.e4236e44e1d78p-4,  
		 -0xf.f9fbffba64d38p-4, -0xf.ffdbe9f9d1298p-4, -0xf.f5bd4d9636c58p-4, -0xf.dbaa48541e16p-4,  
		 -0xf.b1bceb0c9b4e8p-4, -0xf.781f1f9ea151p-4, -0xf.2f0a7f08a0c6p-4, -0xe.d6c817d456438p-4,  
		 -0xe.6fb0250e56e7p-4, -0xd.fa29b6125dc28p-4, -0xd.76aa47848677p-4, -0xc.e5b54dde73258p-4,  
		 -0xc.47dbb205c6d1p-4, -0xb.9dbb406f52bcp-4, -0xa.e7fe0b5fc7868p-4, -0xa.2759c0e79c358p-4,  
		 -0x9.5c8ef544210fp-4, -0x8.8868625b4e1d8p-4, -0x7.abba1d12c17cp-4, -0x6.c760c14c8585cp-4,  
		 -0x5.dc40955d9085p-4, -0x4.eb44a5da74f6p-4, -0x3.f55dda9e62aeep-4, -0x2.fb8205f75e56ap-4,  
		 -0x1.feaaeee86ee36p-4, -0xf.fd557776a76d8p-8, 0x0p0, 0xf.fd557776a76d8p-8,  
		 0x1.feaaeee86ee36p-4, 0x2.fb8205f75e56ap-4, 0x3.f55dda9e62aeep-4, 0x4.eb44a5da74f6p-4,  
		 0x5.dc40955d9085p-4, 0x6.c760c14c8585cp-4, 0x7.abba1d12c17cp-4, 0x8.8868625b4e1d8p-4,  
		 0x9.5c8ef544210fp-4, 0xa.2759c0e79c358p-4, 0xa.e7fe0b5fc7868p-4, 0xb.9dbb406f52bcp-4,  
		 0xc.47dbb205c6d1p-4, 0xc.e5b54dde73258p-4, 0xd.76aa47848677p-4, 0xd.fa29b6125dc28p-4,  
		 0xe.6fb0250e56e7p-4, 0xe.d6c817d456438p-4, 0xf.2f0a7f08a0c6p-4, 0xf.781f1f9ea151p-4,  
		 0xf.b1bceb0c9b4e8p-4, 0xf.dbaa48541e16p-4, 0xf.f5bd4d9636c58p-4, 0xf.ffdbe9f9d1298p-4,  
		 0xf.f9fbffba64d38p-4, 0xf.e4236e44e1d78p-4, 0xf.be680c58c123p-4, 0xf.88ef923316d68p-4,  
		 0xf.43ef73d979c98p-4, 0xe.efacabaa7222p-4, 0xe.8c7b7568da23p-4, 0xe.1abefa072012p-4,  
		 0xd.9ae8ec8696ebp-4, 0xd.0d79184dee53p-4, 0xc.72fce16871d48p-4, 0xb.cc0eb72dc5888p-4,  
		 0xb.195579df6c9d8p-4, 0xa.5b83d3d66f668p-4, 0x9.935786e7e5588p-4, 0x8.c198aeb2ef95p-4,  
		 0x7.e718f895c631cp-4, 0x7.04b2d211d289p-4, 0x6.1b488e7060fecp-4, 0x5.2bc38482335b4p-4,  
		 0x4.3713255c21484p-4, 0x3.3e2c0cf9f6e7p-4, 0x2.42070db6daab6p-4, 0x1.43a0378fadb65p-4,  
		 0x4.3f5dc280f96a8p-8};
	static const double tb_cos[] = // tabulate values of cos(i2^-4) for i in [-2^4pi, 2^4pi]
		{-0xf.ff6fa88a0b1ap-4, -0xf.f331f2c377a9p-4, -0xf.d7025f42f2e9p-4, -0xf.aafd1b42c52p-4,
		 -0xf.6f4e285bf2c78p-4, -0xf.243130884e3a8p-4, -0xe.c9f14a7d768ap-4, -0xe.60e8ae9c638ep-4,
		 -0xd.e9805cc089618p-4, -0xd.642fb348bc9fp-4, -0xc.d17bf7c2c5be8p-4, -0xc.31f7d1b0ee128p-4,
		 -0xb.8642b7eeb5b38p-4, -0xa.cf08514741768p-4, -0xa.0cffc8dcdd368p-4, -0x9.40eb170d1c9cp-4,
		 -0x8.6b963f88a709p-4, -0x7.8dd6856086ae4p-4, -0x6.a88995d4dc814p-4, -0x5.bc94aaba1896cp-4,
		 -0x4.cae3a5523f384p-4, -0x3.d468227f4e51ep-4, -0x2.da18893a7d314p-4, -0x1.dcef1441cb33cp-4,
		 -0xd.de8d7f21b4f78p-8, 0x2.1fb3ab8127beap-8, 0x1.21bd54fc5f9a7p-4, 0x2.205dca0ffebep-4,
		 0x3.1cde0eb530facp-4, 0x4.1641b7b14df68p-4, 0x5.0b8f7622f6548p-4, 0x5.fbd210bc2f334p-4,
		 0x6.e61958e741658p-4, 0x7.c97b1ae14bf08p-4, 0x8.a51407da8346p-4, 0x9.780899321072p-4,
		 0xa.4185ebea67598p-4, 0xb.00c2937ab1ef8p-4, 0xb.b4ff632a908f8p-4, 0xc.5d882d2ee48p-4,
		 0xc.f9b476c897c28p-4, 0xd.88e820b15263p-4, 0xe.0a94032dbea8p-4, 0xe.7e367d2956cf8p-4,
		 0xe.e35bf5ccac89p-4, 0xf.399f500c9eap-4, 0xf.80aa4fbef7508p-4, 0xf.b835efcf670ep-4,
		 0xf.e00aa93eade98p-4, 0xf.f800aaa4fa698p-4, 0x1p+0, 0xf.f800aaa4fa698p-4,
		 0xf.e00aa93eade98p-4, 0xf.b835efcf670ep-4, 0xf.80aa4fbef7508p-4, 0xf.399f500c9eap-4,
		 0xe.e35bf5ccac89p-4, 0xe.7e367d2956cf8p-4, 0xe.0a94032dbea8p-4, 0xd.88e820b15263p-4,
		 0xc.f9b476c897c28p-4, 0xc.5d882d2ee48p-4, 0xb.b4ff632a908f8p-4, 0xb.00c2937ab1ef8p-4,
		 0xa.4185ebea67598p-4, 0x9.780899321072p-4, 0x8.a51407da8346p-4, 0x7.c97b1ae14bf08p-4,
		 0x6.e61958e741658p-4, 0x5.fbd210bc2f334p-4, 0x5.0b8f7622f6548p-4, 0x4.1641b7b14df68p-4,
		 0x3.1cde0eb530facp-4, 0x2.205dca0ffebep-4, 0x1.21bd54fc5f9a7p-4, 0x2.1fb3ab8127beap-8,
		 -0xd.de8d7f21b4f78p-8, -0x1.dcef1441cb33cp-4, -0x2.da18893a7d314p-4, -0x3.d468227f4e51ep-4,
		 -0x4.cae3a5523f384p-4, -0x5.bc94aaba1896cp-4, -0x6.a88995d4dc814p-4, -0x7.8dd6856086ae4p-4,
		 -0x8.6b963f88a709p-4, -0x9.40eb170d1c9cp-4, -0xa.0cffc8dcdd368p-4, -0xa.cf08514741768p-4,
		 -0xb.8642b7eeb5b38p-4, -0xc.31f7d1b0ee128p-4, -0xc.d17bf7c2c5be8p-4, -0xd.642fb348bc9fp-4,
		 -0xd.e9805cc089618p-4, -0xe.60e8ae9c638ep-4, -0xe.c9f14a7d768ap-4, -0xf.243130884e3a8p-4,
		 -0xf.6f4e285bf2c78p-4, -0xf.aafd1b42c52p-4, -0xf.d7025f42f2e9p-4, -0xf.f331f2c377a9p-4,
		 -0xf.ff6fa88a0b1ap-4};
	double xd = to_double (x);
	double k = __builtin_roundeven(xd * inv_2pi);
	double xp = __builtin_fma(minus_2pi, k, xd);
	// x = 2kpi + xp
	double i = __builtin_roundeven(0x1p4 * xp);
	double xpp = __builtin_fma(i, -0x1p-4, xp);
	double xpp2 = xpp*xpp;
	// x = 2kpi + i2^-4 + xpp
	// sin(x) = sin(i2^-4 + xpp) = sin(i2^-4)cos(xpp) + sin(xpp)cos(i2^-4)
	return __builtin_fma(tb_sin[(int)i+50], __builtin_fma(-0.5, xpp2, 1.0), __builtin_fma(-0x1.5555p-3, xpp*xpp2, xpp) * tb_cos[(int)i+50]);
}

// dummy function since GNU libc does not provide it
_Float16 sinf16 (_Float16 x) {
  return (_Float16) sinf ((float) x);
}
