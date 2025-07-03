/* Correctly-rounded power function for binary16 value.

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
#include <errno.h>
#include <math.h> // only used during performance tests

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {_Float16 f; uint16_t u;} b16u16_u;
typedef union {double f; uint64_t u;} b64u64_u;
static const b16u16_u poszero = {.u = 0x0000};
static const b16u16_u negzero = {.u = 0x8000};
static const b16u16_u one = {.u = 0x3c00};
static const b16u16_u neginf = {.u = 0xfc00};
static const b16u16_u posinf = {.u = 0x7c00};

static inline int isint(b16u16_u v) {
	if (v.f == 0.0f) return 1;
	return (v.u & 0x7fff) >> 10 >= 25 - __builtin_ctz(0x400 | v.u);
}

static inline int isodd(b16u16_u v) {
	if (v.f == 0.0f) return 0;
	return (v.u & 0x7fff) >> 10 == 25 - __builtin_ctz(0x400 | v.u);
}

static inline double log_in_pow(double x) {
	b64u64_u xd = {.f = x};
	static const double log2 = 0x1.62e42fefa39efp-1;
	static const double tb[] = // tabulate value of log(1 + i2^-5) for i in [0, 31]
		{0x0p+0, 0x7.e0a6c39e0ccp-8, 0xf.85186008b153p-8, 0x1.6f0d28ae56b4cp-4,
		 0x1.e27076e2af2e6p-4, 0x2.52aa5f03fea46p-4, 0x2.bfe60e14f27a8p-4, 0x3.2a4b539e8ad68p-4,
		 0x3.91fef8f353444p-4, 0x3.f7230dabc7c56p-4, 0x4.59d72aeae9838p-4, 0x4.ba38aeb8474c4p-4,
		 0x5.1862f08717b08p-4, 0x5.746f6fd602728p-4, 0x5.ce75fdaef401cp-4, 0x6.268ce1b05096cp-4,
		 0x6.7cc8fb2fe613p-4, 0x6.d13ddef323d8cp-4, 0x7.23fdf1e6a6888p-4, 0x7.751a813071284p-4,
		 0x7.c4a3d7ebc1bb4p-4, 0x8.12a952d2e87f8p-4, 0x8.5f39721295418p-4, 0x8.aa61e97a6af5p-4,
		 0x8.f42faf382068p-4, 0x9.3caf0944d88d8p-4, 0x9.83eb99a7885fp-4, 0x9.c9f069ab150dp-4,
		 0xa.0ec7f4233957p-4, 0xa.527c2ed81f5d8p-4, 0xa.9516932de2d58p-4, 0xa.d6a0261acf968p-4};
	static const double tl[] = // tabulate value of 1 / (1 + i2^-5) for i in [0, 31]
		{0x1p-52, 0xf.83e0f83e0f84p-56, 0xf.0f0f0f0f0f0fp-56, 0xe.a0ea0ea0ea0e8p-56,
		 0xe.38e38e38e38ep-56, 0xd.d67c8a60dd68p-56, 0xd.79435e50d794p-56, 0xd.20d20d20d20dp-56,
		 0xc.cccccccccccdp-56, 0xc.7ce0c7ce0c7dp-56, 0xc.30c30c30c30cp-56, 0xb.e82fa0be82fap-56,
		 0xb.a2e8ba2e8ba3p-56, 0xb.60b60b60b60b8p-56, 0xb.21642c8590b2p-56, 0xa.e4c415c9882b8p-56,
		 0xa.aaaaaaaaaaaa8p-56, 0xa.72f05397829c8p-56, 0xa.3d70a3d70a3d8p-56, 0xa.0a0a0a0a0a0ap-56,
		 0x9.d89d89d89d8ap-56, 0x9.a90e7d95bc608p-56, 0x9.7b425ed097b4p-56, 0x9.4f2094f2094fp-56,
		 0x9.249249249249p-56, 0x8.fb823ee08fb8p-56, 0x8.d3dcb08d3dcbp-56, 0x8.ad8f2fba93868p-56,
		 0x8.8888888888888p-56, 0x8.64b8a7de6d1d8p-56, 0x8.421084210842p-56, 0x8.208208208208p-56};
	int expo = (xd.u >> 52) - 1023;
	int i = (xd.u & (0x1full << 47)) >> 47;
	xd.f = (xd.u * 0x7fffffffffff) * tl[i];
	xd.f *= __builtin_fma(__builtin_fma(0x1.5555555555555p-2, xd.f, -0.5), xd.f, 1.0);
	return __builtin_fma(log2, (double) expo, tb[i] + xd.f);
}

static inline double exp_in_pow(double x) {
	b64u64_u xd = {.f = x};
	static const b64u64_u x0 = {.f = -0x1.154p+4};
	static const b64u64_u x1 = {.f = 0x1.62cp+3};
	static const double tb[] = // tabulate value of exp(log(2)*i/32) for i in [0, 31]
		{0x1p+0, 0x1.059b0d3158574p+0, 0x1.0b5586cf9890fp+0, 0x1.11301d0125b51p+0,
		 0x1.172b83c7d517bp+0, 0x1.1d4873168b9aap+0, 0x1.2387a6e756238p+0, 0x1.29e9df51fdee1p+0,
		 0x1.306fe0a31b715p+0, 0x1.371a7373aa9cbp+0, 0x1.3dea64c123422p+0, 0x1.44e086061892dp+0,
		 0x1.4bfdad5362a27p+0, 0x1.5342b569d4f82p+0, 0x1.5ab07dd485429p+0, 0x1.6247eb03a5585p+0,
		 0x1.6a09e667f3bccp+0, 0x1.71f75e8ec5f74p+0, 0x1.7a11473eb0187p+0, 0x1.82589994cce13p+0,
		 0x1.8ace5422aa0dbp+0, 0x1.93737b0cdc5e5p+0, 0x1.9c49182a3f09p+0, 0x1.a5503b23e255dp+0,
		 0x1.ae89f995ad3adp+0, 0x1.b7f76f2fb5e47p+0, 0x1.c199bdd85529cp+0, 0x1.cb720dcef9069p+0,
		 0x1.d5818dcfba487p+0, 0x1.dfc97337b9b5fp+0, 0x1.ea4afa2a490d9p+0, 0x1.f50765b6e4541p+0};
#ifdef CORE_MATH_SUPPORT_ERRNO
	if (xd.f > x1.f || xd.f < -0x1.368p+3)
		errno = ERANGE;
#endif
	if (xd.u & (0x7ffull << 52)) { // if x is NaN or x is inf
		if (xd.u == 0xffull << 25) return 0.0; // x is -Inf
		else return xd.f + xd.f;
	}
	else if (xd.u > x0.u) return 0x1p-25;
	else if (xd.f > x1.f) return 0x1.ffcp+15 + 0x1p+5;
	static const double thirtytwo_over_log2 = 0x1.71547652b82fep+5;
	static const double minus_log2_over_thirtytwo = -0x1.62e42fefa39efp-6;
	double j = __builtin_roundeven(thirtytwo_over_log2 * xd.f);
	int64_t jint = j;
	int i = jint & 0x1f;
	double xp = __builtin_fma(minus_log2_over_thirtytwo, j, xd.f);
	xp = __builtin_fma(__builtin_fma(0.5, xp, 1.0), xp, 1.0);
	b64u64_u ret = {.f = xp * tb[i]};
	ret.u += (jint >> 5) * (1ull << 52);
	return ret.f;
}

_Float16 cr_powf16(_Float16 x, _Float16 y){
	b16u16_u vx = {.f = x}, vy = {.f = y};
	uint64_t sign = 0;
	if ((vx.u & 0x7fff) == 0x3c00) { // |x| = 1
		if (vx.u >> 15) { // x = -1
			if ((vy.u & 0x7fff) > 0x7c00) return y + y; // y = NaN
			if (isint(vy)) return (isodd(vy)) ? vx.f : -vx.f;
#ifdef CORE_MATH_SUPPORT_ERRNO
			errno = EDOM;
#endif
			return 0.0f / 0.0f;
		}
		return ((vy.u & 0x7fff) > 0x7c000 && !(vy.u & 0x0200)) ? y + y : x;
		// 1^y = y except if y = sNaN
	}
	if (!(vy.u & 0x7fff)) // y = 0
		return ((vy.u & 0x7fff) > 0x7c000 && !(vy.u & 0x0200)) ? x + x : one.f;
		// x^0 = 1 except if x = sNaN
	if ((vy.u & 0x7fff) >= 0x7c00) { // y = Inf/NaN
		// the case |x| = 1 was checked above
		if ((vx.u & 0x7fff) > 0x7c00) return x + x; // x = NaN
		if ((vy.u & 0x7fff) == 0x7c00) { // y = +/-Inf
			if (((vx.u & 0x7fff) < 0x3c00) ^ (vy.u >> 15)) {
				return posinf.f; // |x| < 1 && y = -Inf or |x| > 1 && y = +Inf
			} else {
				return poszero.f; // |x| < 1 && y = +Inf or |x| > 1 && y = -Inf
			}
		}
		return y + y; // y = NaN
	}
	if (!(vx.u & 0x7fff)) { // if x = 0
		if (vy.u >> 15) { // y < 0
#ifdef CORE_MATH_SUPPORT_ERRNO
			errno = ERANGE;
#endif
			if (isodd(vy) && vx.u >> 1) return neginf.f;
			else return posinf.f;
		} else { // y > 0
			if (isodd(vy) && vx.u >> 1) return negzero.f;
			else return poszero.f;
		}
	}
	if (vx.u >= 0x7c00) { // x = Inf or x = NaN or x <= 0
		if ((vx.u & 0x7fff) == 0x7c00) { // x = +/-Inf
			if (!isodd(vy)) vx.u &= 0x7fff; // y even -> ret will be positive
			if (vy.u >> 15) vx.u &= 0x8000; // y < 0 -> ret will be +/-0
			return vx.f;
		}
		if ((vx.u & 0x7fff) > 0x7c00) return x + x; // x is NaN
		// x < 0
		if (!isint(vy)) {
#ifdef CORE_MATH_SUPPORT_ERRNO
			errno = EDOM;
#endif
			return 0.0f / 0.0f;
		} 
		else if (isodd(vy)) sign = 1ull << 63;
		vx.u &= 0x7fff;
	}
	b64u64_u ret = {.f = exp(log(x) * y)};
	ret.u += sign;
	return ret.f;
}

// dummy function since GNU libc does not provide it
_Float16 powf16 (_Float16 x, _Float16 y) {
	return (_Float16) powf ((float) x, (float) y);
}
