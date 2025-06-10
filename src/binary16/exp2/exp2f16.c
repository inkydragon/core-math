/* Correctly-rounded natural exponential function for binary32 value.

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

_Float16 cr_exp2f16(_Float16 x){
	uint16_t x0 = 0xce40; // binary representation of x0 in order to compare uint16_t rather than flt16
 	_Float16 x1 = 0x1.ffcp3; // largest _Float16 such that exp(x1) <= MAX_FLOAT16 < exp(x1+)
	static const float tb[] = // tabulate value of 2^(i/2^6) for i in [-2^5, 2^5], size(tb) = 2^6 + 1
		{0xb.504f3p-4, 0xb.6fd92p-4, 0xb.8fbafp-4, 0xb.aff5bp-4,
		 0xb.d08a4p-4, 0xb.f179ap-4, 0xc.12c4dp-4, 0xc.346cdp-4,
		 0xc.5672ap-4, 0xc.78d75p-4, 0xc.9b9bep-4, 0xc.bec15p-4,
		 0xc.e248cp-4, 0xd.06334p-4, 0xd.2a81ep-4, 0xd.4f35bp-4,
		 0xd.744fdp-4, 0xd.99d16p-4, 0xd.bfbb8p-4, 0xd.e60f5p-4,
		 0xe.0ccdfp-4, 0xe.33f89p-4, 0xe.5b907p-4, 0xe.8396ap-4,
		 0xe.ac0c7p-4, 0xe.d4f3p-4, 0xe.fe4bap-4, 0xf.281778p-4,
		 0xf.5257dp-4, 0xf.7d0dfp-4, 0xf.a83b3p-4, 0xf.d3e0cp-4,
		 0x1p+0, 0x1.02c9a4p+0, 0x1.059b0ep+0, 0x1.087452p+0,
		 0x1.0b5586p+0, 0x1.0e3ec4p+0, 0x1.11301ep+0, 0x1.1429aap+0,
		 0x1.172b84p+0, 0x1.1a35bep+0, 0x1.1d4874p+0, 0x1.2063bap+0,
		 0x1.2387a6p+0, 0x1.26b456p+0, 0x1.29e9ep+0, 0x1.2d285ap+0,
		 0x1.306fep+0, 0x1.33c08cp+0, 0x1.371a74p+0, 0x1.3a7db4p+0,
		 0x1.3dea64p+0, 0x1.4160a2p+0, 0x1.44e086p+0, 0x1.486a2cp+0,
		 0x1.4bfdaep+0, 0x1.4f9b28p+0, 0x1.5342b6p+0, 0x1.56f474p+0,
		 0x1.5ab07ep+0, 0x1.5e76f2p+0, 0x1.6247ecp+0, 0x1.662388p+0,
		 0x1.6a09e6p+0};
	b16u16_u v = {.f = x};
	if (v.u > x0) return (_Float16) 0x1p-25f;
	else if (x > x1) return (_Float16) 0x1.ffcp15f + 0x1p4f; 
	else if (v.u == 0x11c5) return (_Float16) 0x1.004p+0 - 0x1p-12; // only one wrong case
	else {		
		float k = __builtin_roundevenf((float) x); 
		float xp = -k + x;
		int i = 0x40 * xp;
		float xpp = __builtin_fmaf((float) i, -0x1p-6, xp); // x = k + i/2^6 + xpp
																		   									// So, 2^x = 2^k * 2^(i/2^6) * 2^xpp

		// result
    xpp = 1.0 + xpp * (0x1.62e43p-1 + xpp * (0x1.ebfbep-3 + xpp * 0x1.c64d84p-5));
		b32u32_u w = {.f = xpp * tb[i + 32]};
    w.u += (int) k << 23;
  	return w.f;
	}
}

// dummy function since GNU libc does not provide it
_Float16 exp2f16 (_Float16 x) {
  return (_Float16) exp2f ((float) x);
}
