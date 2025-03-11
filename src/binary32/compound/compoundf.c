/* Correctly-rounded compound function for binary32 values.

Copyright (c) 2025 Alexei Sibidanov and Paul Zimmermann

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
#include <fenv.h>
#ifdef __x86_64__
#include <x86intrin.h>
#define FLAG_T uint32_t
#else
#define FLAG_T fexcept_t
#endif

// Warning: clang also defines __GNUC__
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#pragma STDC FENV_ACCESS ON

typedef union {float f; uint32_t u;} b32u32_u;
typedef union {double f; uint64_t u;} b64u64_u;

static inline int issignalingf(float x) {
  b32u32_u u = {.f = x};
  /* To keep the following comparison simple, toggle the quiet/signaling bit,
   so that it is set for sNaNs.  This is inverse to IEEE 754-2008 (as well as
   common practice for IEEE 754-1985).  */
  u.u ^= 0x00400000;
  /* We have to compare for greater (instead of greater or equal), because x's
     significand being all-zero designates infinity not NaN.  */
  return (u.u & 0x7fffffff) > 0x7fc00000;
}

static inline int isodd(float y){
  b32u32_u wy = {.f = y};
  int ey = ((wy.u>>23) & 0xff) - 127, s = ey + 9, odd = 0;
  if(ey>=0){
    if(s<32 && !(wy.u<<s)) odd = (wy.u>>(32-s))&1;
    if(s==32) odd = wy.u&1;
  }
  return odd;
}

static inline int isint(float y){
  b32u32_u wy = {.f = y};
  int ey = ((wy.u>>23) & 0xff) - 127, s = ey + 9;
  if(ey>=0){
    if(s>=32) return 1;
    return !(wy.u<<s);
  }
  if(!(wy.u<<1)) return 1;
  return 0;
}

/* __builtin_roundeven was introduced in gcc 10:
   https://gcc.gnu.org/gcc-10/changes.html,
   and in clang 17 */
#if ((defined(__GNUC__) && __GNUC__ >= 10) || (defined(__clang__) && __clang_major__ >= 17)) && (defined(__aarch64__) || defined(__x86_64__) || defined(__i386__))
# define roundeven_finite(x) __builtin_roundeven (x)
#else
/* round x to nearest integer, breaking ties to even */
static double
roundeven_finite (double x)
{
  double ix;
# if (defined(__GNUC__) || defined(__clang__)) && (defined(__AVX__) || defined(__SSE4_1__) || (__ARM_ARCH >= 8))
#  if defined __AVX__
   __asm__("vroundsd $0x8,%1,%1,%0":"=x"(ix):"x"(x));
#  elif __ARM_ARCH >= 8
   __asm__ ("frintn %d0, %d1":"=w"(ix):"w"(x));
#  else /* __SSE4_1__ */
   __asm__("roundsd $0x8,%1,%0":"=x"(ix):"x"(x));
#  endif
# else
  ix = __builtin_round (x); /* nearest, away from 0 */
  if (__builtin_fabs (ix - x) == 0.5)
  {
    /* if ix is odd, we should return ix-1 if x>0, and ix+1 if x<0 */
    union { double f; uint64_t n; } u, v;
    u.f = ix;
    v.f = ix - __builtin_copysign (1.0, x);
    if (__builtin_ctz (v.n) > __builtin_ctz (u.n))
      ix = v.f;
  }
# endif
  return ix;
}
#endif

// INVLOG2 = 1/log(2) * (1 + eps1) with |eps1| < 2^-55.976
#define INVLOG2 0x1.71547652b82fep+0 
// LOG2 = log(2) * (1 + eps) with |eps| < 2^-54.730
#define LOG2 0x1.62e42fefa39efp-1

/* for x > -1, approximates log2(1+x)
   Note: assume x is representable in binary32
   For |x| < 2^-29, the relative error is bounded by 2^-50.950.
*/
static double _log2p1 (double x)
{
  /* for x > 0, 1+x is exact when 2^-29 <=  x < 2^53
     for x < 0, 1+x is exact when -1 < x <= 2^-30 */

  if (__builtin_expect (__builtin_fabs (x) < 0x1p-29, 0)) { // |x| < 2^-29
    // |log2(1+x) - 1/log(2) * (x - x^2/2)| < 2^-59.584 * |log2(1+x)|
    double t = x - (x * x) * 0.5;
    /* since x is epresentable in binary32, x*x is exact, and so is (x * x) * 0.5.
       Thus the only error in the computation of t is the final rounding, which
       is bounded by ulp(t): t = (x - x^2/2) * (1 + eps2) with |eps2| < 2^-52
    */
    return INVLOG2 * t;
    /* since INVLOG2 = 1/log(2) * (1 + eps1) and
       and   t = (x - x^2/2) * (1 + eps2)
       let u = o(INVLOG2 * t) then u = INVLOG2 * t * (1 + eps3) with |eps3|<2^-53
       thus u = 1/log(2) * (x - x^2/2) * (1 + eps1)*(1 + eps2)*(1 + eps3)
              = 1/log(2) * (x - x^2/2) * (1 + eps4) with |eps4| < 2^-50.954
       Now Sollya says the relative error by approximating log2(1+x) by
       1/log(2) * (x - x^2/2) for |x| < 2^-29 is bounded by 2^-59.584
       (file compoundf.sollya), thus:
       u = log2(1+x) * (1+eps4)*(1+eps5) with |eps5| < 2^-59.584
         = log2(1+x) * (1+eps6) with |eps6| < 2^-50.950 */
  }

  double u = (x >= 0x1p53) ? x : 1.0 + x;
  /* For x < 0x1p53, x + 1 is exact thus u = x+1.
     For x >= 2^53, we estimate log2(x) instead of log2(1+x),
     since log2(1+x) = log2(x) + log2(1+1/x),
     log2(x) >= 53 and |log2(1+1/x)| < 2^-52.471, the additional relative
     error is bounded by 2^-52.471/53 < 2^-58.198 */

  b64u64_u v = {.f = u};
  uint64_t m = v.u & 0xfffffffffffffull;
  int e = (v.u >> 52) - 0x3ff + (m >= 0x6a09e667f3bcdull);
  // 2^e/sqrt(2) < u < 2^e*sqrt(2)
  v.u = ((0x3ffull + e) << 52) | ml;
  double t = v.f;
  // u = 2^e*t with 1/sqrt(2) < t < sqrt(2)
  // thus log2(u) = e + log2(t)
  int i = (v.u >> 48) & 0x1f; // 6 <= i <= 22
  /* for 6 <= i < 16, inv[i-6] approximates 1/t for 1/2+i/32 <= t < 1/2+(i+1)/32,
     for 16 <= i <= 22, inv[i-6] approximates 1/t for i/16 <= t < (i+1)/16. */
  static const double inv[] = {
    0x1.7p+0, 0x1.6p+0, 0x1.5p+0, 0x1.4p+0, 0x1.38p+0, 0x1.28p+0, 0x1.2p+0,
    0x1.18p+0, 0x1.1p+0, 0x1.08p+0, 0x1.fp-1, 0x1.dp-1, 0x1.cp-1, 0x1.ap-1,
    0x1.9p-1, 0x1.8p-1, 0x1.7p-1,
  };
  double r = inv[i-6];
  double z = __builtin_fma (r, t, -1.0); // exact, -19/512 <= z <= 20/512
}

float cr_compoundf (float x, float y)
{
  /* Rules from IEEE 754-2019 for compound (x, n) with n integer:
     compound (x, 0) is 1 for x >= −1 or quiet NaN
     compound (−1, n) is +Inf and signals the divideByZero exception for n < 0
     compound (−1, n) is +0 for n > 0
     compound (+/-0, n) is 1
     compound (+Inf, n) is +Inf for n > 0
     compound (+Inf, n) is +0 for n < 0
     compound (x, n) is qNaN and signals the invalid exception for x < −1
     compound (qNaN, n) is qNaN for n <> 0.
  */
  double xd = x, yd = y;
  b64u64_u tx = {.f = xd}, ty = {.f = yd};

  if(__builtin_expect (x == 0, 0)) // compound(0,y) = 1 except for y = sNaN
    return issignalingf (y) ? x + y : 1.0f;

  if(__builtin_expect (ty.u<<1 == 0, 0)) // compound(x,0) = 1 except for x = sNaN
    return issignalingf (x) ? x + y : 1.0f;

  if(__builtin_expect ((ty.u<<1) >= (uint64_t)0x7ff<<53, 0)){ // y=Inf/NaN
    // the case x=0 was already checked above
    if((tx.u<<1) > (uint64_t)0x7ff<<53) return x + y; // x=NaN
    if((ty.u<<1) == (uint64_t)0x7ff<<53) { // y=+/-Inf
      int sy = ty.u>>63; // sign bit of y
      if (-2.0f < x && x < 0 && sy == 0)
	return 0;
      if ((x < -2.0f || 0 < x) && sy != 0)
        return 0;
      if (x == -2.0f)
        return 1.0f;
      return __builtin_inf();
    }
    return x + y; // case y=NaN
  }

  uint64_t minus_two = 0xc000000000000000ull; // encoding of -2
  if (__builtin_expect (tx.u >= minus_two || tx.u >> 52 == 0x7ff, 0)){
    // x is Inf, NaN or less than -1
    if ((tx.u<<1) == (uint64_t)0x7ff<<53){ // x is +Inf or -Inf
      if (!isodd(y)) x = __builtin_fabsf (x);
      return (ty.u>>63) ? 1.0f/x : x;
    }
    if ((tx.u<<1) > (uint64_t)0x7ff<<53) return x + y; // x is NaN
    if (x < -1.0) // x <= 0
      if (!isint(y) && x != 0) {
#ifdef CORE_MATH_SUPPORT_ERRNO
        errno = EDOM;
#endif
	return (x - x) / (x - x);  // NaN, should raise 'Invalid operation' exception.
      }
  }

  if (__builtin_expect (x == -1.0f, 0)) { // x=-1
    if (ty.u>>63){ // y < 0
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE;
#endif
      return (isodd(y)) ? 1.0f/__builtin_copysignf(0.0f,x) : 1.0f/0.0f;
    } else // y > 0
      return (isodd(y)) ? __builtin_copysignf(1.0f,x)*0.0f : 0.0f;
  }

  // now x > -1

  double l = _log2p1 (tx.f); // approximates log2(1+x)

  return x;
}
