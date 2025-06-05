/* Correctly-rounded natural exponential function for binary32 value.

Copyright (c) 2025 Maxence Ponsardin.

This file is part of the CORE-MATH project
(https://core-math.gitlabpages.inria.fr/).

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software isref=0x1.02p-15
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
#include <math.h>

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;

_Float16 cr_expf16(_Float16 x){
	// _Float16 x0 = -0x1.154p4; // smallest _Float16 such that exp(x0-) < MIN_FLOAT16 / 2 <= exp(x0)
	uint16_t x0 = 0xcc55; // binaryrepresentation of x0 in order to compare uint16_t rather than flt16
 	_Float16 x1 = 0x1.62bp3; // largest _Float16 such that exp(x1) <= MAX_FLOAT16 < exp(x1+)
	static const float tb[] = // tabulate value of exp(i/2^6) for i in [-2^6, 2^6], size(tb) = 2^7 + 1
		{0x5.e2d59p-4, 0x5.fa9038p-4, 0x6.12aa9p-4, 0x6.2b261p-4, 0x6.44044p-4, 0x6.5d46b8p-4, 0x6.76efp-4, 0x6.90feb8p-4,
		 0x6.ab778p-4, 0x6.c65bp-4, 0x6.e1aae8p-4, 0x6.fd68fp-4, 0x7.1996c8p-4, 0x7.363638p-4, 0x7.53491p-4, 0x7.70d12p-4,
		 0x7.8ed038p-4, 0x7.ad484p-4, 0x7.cc3b2p-4, 0x7.ebaacp-4, 0x8.0b992p-4, 0x8.2c083p-4, 0x8.4cfa1p-4, 0x8.6e70cp-4,
		 0x8.906e5p-4, 0x8.b2f4ep-4, 0x8.d606bp-4, 0x8.f9a5dp-4, 0x9.1dd4ap-4, 0x9.42954p-4, 0x9.67eap-4, 0x9.8dd55p-4,
		 0x9.b4598p-4, 0x9.db78fp-4, 0xa.03361p-4, 0xa.2b937p-4, 0xa.54939p-4, 0xa.7e38fp-4, 0xa.a8863p-4, 0xa.d37ep-4,
		 0xa.ff231p-4, 0xb.2b78p-4, 0xb.587fcp-4, 0xb.863cfp-4, 0xb.b4b29p-4, 0xb.e3e38p-4, 0xc.13d2bp-4, 0xc.44832p-4,
		 0xc.75f7dp-4, 0xc.a833ep-4, 0xc.db3a8p-4, 0xd.0f0edp-4, 0xd.43b41p-4, 0xd.792d8p-4, 0xd.af7e9p-4, 0xd.e6aaap-4,
		 0xe.1eb51p-4, 0xe.57a17p-4, 0xe.91735p-4, 0xe.cc2e4p-4, 0xf.07d6p-4, 0xf.446e3p-4, 0xf.81fabp-4, 0xf.c07f5p-4,
		 0x1p+0, 0x1.04080ap+0, 0x1.082056p+0, 0x1.0c4924p+0, 0x1.1082b6p+0, 0x1.14cd5p+0, 0x1.192938p+0, 0x1.1d96bp+0,
		 0x1.221604p+0, 0x1.26a77ap+0, 0x1.2b4b58p+0, 0x1.3001ecp+0, 0x1.34cb82p+0, 0x1.39a862p+0, 0x1.3e98dep+0, 0x1.439d44p+0,
		 0x1.48b5e4p+0, 0x1.4de30ep+0, 0x1.532518p+0, 0x1.587c54p+0, 0x1.5de918p+0, 0x1.636bbap+0, 0x1.690492p+0, 0x1.6eb3fcp+0,
		 0x1.747a52p+0, 0x1.7a57eep+0, 0x1.804d3p+0, 0x1.865a78p+0, 0x1.8c8024p+0, 0x1.92be9ap+0,0x1.99163ap+0, 0x1.9f876ep+0,
		 0x1.a61298p+0, 0x1.acb826p+0, 0x1.b3787ep+0, 0x1.ba540ep+0, 0x1.c14b44p+0, 0x1.c85e8ep+0, 0x1.cf8e5ep+0, 0x1.d6db26p+0,
		 0x1.de455ep+0, 0x1.e5cd7ap+0, 0x1.ed73f2p+0, 0x1.f53942p+0, 0x1.fd1de6p+0, 0x2.05225cp+0, 0x2.0d4724p+0, 0x2.158ccp+0,
		 0x2.1df3b8p+0, 0x2.267c8cp+0, 0x2.2f27c8p+0, 0x2.37f5f8p+0, 0x2.40e7a8p+0, 0x2.49fd64p+0, 0x2.5337c4p+0, 0x2.5c9754p+0,
		 0x2.661cbp+0, 0x2.6fc87p+0, 0x2.799b28p+0, 0x2.83957cp+0, 0x2.8db808p+0, 0x2.980374p+0, 0x2.a2785cp+0, 0x2.ad177p+0, 0x2.b7e15p+0};

	b16u16_u v = {.f = x};
	if (v.u > x0) return (_Float16) 0x1p-25f;
	else if (x > x1) return (_Float16) 0x1.ffcp15f + 0x1p4f; 
	else {
		float minus_log2 = -0x1.62e430p-1;
		float inv_log2 = 0x1.715476p0f;
		float k = __builtin_roundevenf(inv_log2 * x); 
		float xp = __builtin_fmaf(k, minus_log2, x); // xp is a float such that |xp| is minimal and x = klog(2) + xp
		int i = 0x40 * xp;
		float xpp = xp - (float) i / 0x40; // x = klog(2) + i/2^6 + xpp
																		   // So, exp(x) = 2^k * exp(i/2^6) * exp(xpp)

		// result
		xpp = 1.0 + xpp * (1 + xpp * (0.5 + xpp * (0x1.555644p-3 + xpp * 0x1.555bep-5)));
		if (k >= 0) return (_Float16) (xpp * tb[i + (1<<6)] * (1 << (int)k));
		else return (_Float16) (xpp * tb[i + (1<<6)] / (1 << (int) -k));
	}
}



// dummy function since GNU libc does not provide it
_Float16 expf16 (_Float16 x) {
  return (_Float16) expf ((float) x);
}
