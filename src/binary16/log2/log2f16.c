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
#include <math.h> // only used during performance tests

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

<<<<<<< HEAD
<<<<<<< HEAD

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;
<<<<<<< HEAD
static const b16u16_u neginf = {.u = 0xfc00};

/* This code is largely inspired by TANG, P. T. P. Table-driven
implementation of the logarithm function in IEEE floating-point
arithmetic. ACM Trans. Math. Softw. 16, 4 (1990), 378–400.
https://dl.acm.org/doi/10.1145/98267.98294 */

_Float16 cr_log2f16(_Float16 x){
	b16u16_u t = {.f = x};
	if (t.u == 0) return neginf.f;
	else if (t.u >> 10 >= 0x1f) {
		if (t.u == 0x8000) return neginf.f;
=======
=======

>>>>>>> 9f181131 (Fixed almost 4k wrong cases, it was due to subnormalized : I was looking for fraction right after the exponent but in subnormalized number, the fraction iskinda shifted.\nStill 7 wrong cases)
typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {float f; uint32_t u;} b32u32_u;
=======
b16u16_u neginf = {.u = 0xfc00};
>>>>>>> 9a3bb25c (Added logf16 succesfully with 0 wrong case)

_Float16 cr_logf16(_Float16 x){
	b16u16_u t = {.f = x};
	if (t.u == 0) return neginf.f;
	else if (t.u >> 10 >= 0x1f) {
<<<<<<< HEAD
		if (t.u == 0x8000) return -0xffcp+15f - 0x1p5f;
>>>>>>> 1fe06f90 (First implementation of logf16 with 4000 wrong cases)
=======
		if (t.u == 0x8000) return neginf.f;
>>>>>>> 9a3bb25c (Added logf16 succesfully with 0 wrong case)
		else if (t.u >> 15) return 0.0 / 0.0;
		else return x + x;
	}
<<<<<<< HEAD:src/binary16/log/logf16.c
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
	float log2 = 0x1.62e430p-1;
>>>>>>> 9f181131 (Fixed almost 4k wrong cases, it was due to subnormalized : I was looking for fraction right after the exponent but in subnormalized number, the fraction iskinda shifted.\nStill 7 wrong cases)
=======
>>>>>>> f2e5401c (Added log2f16):src/binary16/log2/log2f16.c
	static const float tb[] = 
		{0x0p+0, 0xb.5d69cp-8, 0x1.663f7p-4, 0x2.118b1p-4,  
		 0x2.b80348p-4, 0x3.59ebc4p-4, 0x3.f782d8p-4, 0x4.9101e8p-4,  
		 0x5.269e1p-4, 0x5.b8887p-4, 0x6.46eeap-4, 0x6.d1fbp-4,  
		 0x7.59d4f8p-4, 0x7.dea158p-4, 0x8.60828p-4, 0x8.df989p-4,  
		 0x9.5c01ap-4, 0x9.d5dap-4, 0xa.4d3c2p-4, 0xa.c2411p-4,  
		 0xb.35004p-4, 0xb.a58ffp-4, 0xc.1404fp-4, 0xc.80731p-4,  
		 0xc.eaedp-4, 0xd.53848p-4, 0xd.ba4a4p-4, 0xe.1f4e5p-4,  
		 0xe.829fbp-4, 0xe.e44cdp-4, 0xf.44636p-4, 0xf.a2f04p-4};
	static const float tl[] =
		{0x1p+0, 0xf.83e1p-4, 0xf.0f0f1p-4, 0xe.a0ea1p-4,  
		 0xe.38e39p-4, 0xd.d67c9p-4, 0xd.79436p-4, 0xd.20d21p-4,  
		 0xc.ccccdp-4, 0xc.7ce0cp-4, 0xc.30c31p-4, 0xb.e82fap-4,  
		 0xb.a2e8cp-4, 0xb.60b61p-4, 0xb.21643p-4, 0xa.e4c41p-4,  
		 0xa.aaaabp-4, 0xa.72f05p-4, 0xa.3d70ap-4, 0xa.0a0a1p-4,  
		 0x9.d89d9p-4, 0x9.a90e8p-4, 0x9.7b426p-4, 0x9.4f209p-4,  
		 0x9.24925p-4, 0x8.fb824p-4, 0x8.d3dcbp-4, 0x8.ad8f3p-4,  
		 0x8.88889p-4, 0x8.64ba18p-4, 0x8.42108p-4, 0x8.20821p-4};
	b32u32_u xf = {.f = x};
	int expo = (xf.u >> 23) - 127; // used float instead of flaot16 to avoid working with subnormalized
<<<<<<< HEAD:src/binary16/log/logf16.c
<<<<<<< HEAD
=======
	double log2 = 0x1.62e42fefa39efp-1;
	static const double tb[] = 
		{0x0p+0, 0x3.f815161f807c8p-8, 0x7.e0a6c39e0ccp-8, 0xb.ba2c7b196e7ep-8,  
     0xf.85186008b153p-8, 0x1.341d7961bd1d1p-4, 0x1.6f0d28ae56b4cp-4, 0x1.a926d3a4ad563p-4,  
     0x1.e27076e2af2e6p-4, 0x2.1aefcf9a11cb2p-4, 0x2.52aa5f03fea46p-4, 0x2.89a56d996fa3cp-4,  
     0x2.bfe60e14f27a8p-4, 0x2.f57120421b212p-4, 0x3.2a4b539e8ad68p-4, 0x3.5e7929d017fe6p-4,  
     0x3.91fef8f353444p-4, 0x3.c4e0edc55e5ccp-4, 0x3.f7230dabc7c56p-4, 0x4.28c9389ce438cp-4,  
     0x4.59d72aeae9838p-4, 0x4.8a507ef3de598p-4, 0x4.ba38aeb8474c4p-4, 0x4.e993155a517a8p-4,  
     0x5.1862f08717b08p-4, 0x5.46ab61cb7e0b4p-4, 0x5.746f6fd602728p-4, 0x5.a1b207a6c52bcp-4,  
     0x5.ce75fdaef401cp-4, 0x5.fabe0ee0abf0cp-4, 0x6.268ce1b05096cp-4, 0x6.51e5070845becp-4,  
     0x6.7cc8fb2fe613p-4, 0x6.a73b26a682128p-4, 0x6.d13ddef323d8cp-4, 0x6.fad36769c6dfp-4,  
     0x7.23fdf1e6a6888p-4, 0x7.4cbf9f803af54p-4, 0x7.751a813071284p-4, 0x7.9d109875a1e2p-4,  
     0x7.c4a3d7ebc1bb4p-4, 0x7.ebd623de3cc7cp-4, 0x8.12a952d2e87f8p-4, 0x8.391f2e0e6fap-4,  
     0x8.5f39721295418p-4, 0x8.84f9cf16a64b8p-4, 0x8.aa61e97a6af5p-4, 0x8.cf735a33e4b78p-4,  
     0x8.f42faf382068p-4, 0x9.18986bdf5fa18p-4, 0x9.3caf0944d88d8p-4, 0x9.6074f6a24746p-4,  
     0x9.83eb99a7885fp-4, 0x9.a7144ece70e98p-4, 0x9.c9f069ab150dp-4, 0x9.ec813538ab7d8p-4,  
     0xa.0ec7f4233957p-4, 0xa.30c5e10e2f61p-4, 0xa.527c2ed81f5d8p-4, 0xa.73ec08dbadd88p-4,  
     0xa.9516932de2d58p-4, 0xa.b5fcead9f9cc8p-4, 0xa.d6a0261acf968p-4, 0xa.f70154920b3a8p-4};
	b64u64_u xf = {.f = x};
	int expo = (xf.u >> 52) - 1023; // used double instead of flaot16 to avoid working with subnormalized
>>>>>>> c7042ccb (Tried to use double instead of float to have a better precision but it didn't work, still ~4000 wrong cases)
	int i = (t.u & 0x03f0) >> 4;
	xf.f = 0x1p-10 * (t.u & 0x000f);
	// We have, x = 2^expo * (1 + i2⁻⁶ + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2⁻⁶) + log(1 + xf.f / (1 + i2⁻⁶))
<<<<<<< HEAD
	xf.f /= __builtin_fmaf(0x1p-6, (float) i, (float) 1);
>>>>>>> 1fe06f90 (First implementation of logf16 with 4000 wrong cases)
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.555554p-2, xf.f, -0.5f), xf.f, 1.0f);
	return __builtin_fmaf(log2, (float) expo, tb[i] + xf.f);
=======
	xf.f /= __builtin_fma(0x1p-6, (double) i, 1.0);
	xf.f = __builtin_fma(__builtin_fma(__builtin_fma(__builtin_fma(-0x1.000ccd7266ca8p-2, xf.f, 0x1.556222b13ce3fp-2), xf.f, -0x1.ffffffd10d907p-2), xf.f, 0x1.fffffff332571p-1), xf.f, -0x1p-58);
	return __builtin_fma(log2, (double) expo, tb[i] + xf.f);
>>>>>>> c7042ccb (Tried to use double instead of float to have a better precision but it didn't work, still ~4000 wrong cases)
=======
	int i = (xf.u & 0x007e0000) >> 17;
	xf.f = 0x1p-23 * (xf.u & 0x0001ffff);
	// We have, x = 2^expo * (1 + i2⁻⁶ + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2⁻⁶) + log(1 + xf.f / (1 + i2⁻⁶))
	xf.f /= __builtin_fmaf(0x1p-6, (float) i, (float) 1);
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.555554p-2, xf.f, -0.5f), xf.f, 1.0f);
	return __builtin_fmaf(log2, (float) expo, tb[i] + xf.f);
>>>>>>> 9f181131 (Fixed almost 4k wrong cases, it was due to subnormalized : I was looking for fraction right after the exponent but in subnormalized number, the fraction iskinda shifted.\nStill 7 wrong cases)
=======
	int i = (xf.u & 0x007c0000) >> 18;
	xf.f = 0x1p-23 * (xf.u & 0x0003ffff) * tl[i];
	// We have, x = 2^expo * (1 + i2⁻⁵ + xf.f)
	// Thus, log(x) = expo log(2) + log(1 + i2⁻⁵) + log(1 + xf.f / (1 + i2⁻⁵))
	xf.f *= __builtin_fmaf(__builtin_fmaf(0x1.555554p-2, xf.f, -0.5f), xf.f, 1.0f) * 0x1.715476p+0;
	return (float) expo + tb[i] + xf.f;
>>>>>>> f2e5401c (Added log2f16):src/binary16/log2/log2f16.c
}

// dummy function since GNU libc does not provide it
_Float16 log2f16 (_Float16 x) {
  return (_Float16) log2f ((float) x);
}
