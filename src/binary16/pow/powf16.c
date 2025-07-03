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

int isint(b16u16_u v) {
	if (v.f == 0.0f) return 1;
	return (v.u & 0x7fff) >> 10 >= 25 - __builtin_ctz(0x400 | v.u);
}

int isodd(b16u16_u v) {
	if (v.f == 0.0f) return 0;
	return (v.u & 0x7fff) >> 10 == 25 - __builtin_ctz(0x400 | v.u);
}

_Float16 cr_powf16(_Float16 x, _Float16 y){
	b16u16_u vx = {.f = x}, vy = {.f = y};
	int sign = 0;
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
		else if (isodd(vy)) sign = 1;
		vx.u &= 0x7fff;
	}
	return x + y;
}

// dummy function since GNU libc does not provide it
_Float16 powf16 (_Float16 x, _Float16 y) {
	return (_Float16) powf ((float) x, (float) y);
}
