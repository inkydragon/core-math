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

#include <stdio.h>
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

/* for -1/32 <= z <= 17/512, returns an approximation of log2(1+z)
   with relative error < 2^-43.69 (see p1.sollya) */
static double p1 (double z)
{
  static const double P[] = {0, 0x1.71547652b835bp0, -0x1.71547652e79f8p-1,
    0x1.ec709dbbc7569p-2, -0x1.71546b1ea8be3p-2, 0x1.2777226bea632p-2,
    -0x1.ed1a8e430fc2bp-3, 0x1.a43ba4572b135p-3};
    
  double z2 = z * z;
  double c6 = __builtin_fma (P[7], z, P[6]);
  double c4 = __builtin_fma (P[5], z, P[4]);
  double c2 = __builtin_fma (P[3], z, P[2]);
  double c0 = P[1] * z; // P[0] = 0
  double z4 = z2 * z2;
  c4 = __builtin_fma (c6, z2, c4);
  c0 = __builtin_fma (c2, z2, c0);
  return __builtin_fma (c4, z4, c0);
}

/* for |z| <= 2^-6, returns an approximation of 2^z
   with absolute error < 2^-33.225 (see q1.sollya) */
static double q1 (double z)
{
  //  int bug = z == 0x1.715474e163b99p-21;
  int bug = 0;
  static const double Q[] = {1.0, 0x1.62e42fef46c6bp-1,
                             0x1.ebfce69bff861p-3, 0x1.c6b19adeb965dp-5};
  double z2 = z * z;
  double c2 = __builtin_fma (Q[3], z, Q[2]);
  double c0 = __builtin_fma (Q[1], z, Q[0]);
  return __builtin_fma (c2, z2, c0);
}

/* for x > -1, approximates log2(1+x)
   Note: assume x is representable in binary32
   For |x| < 2^-29, the relative error is bounded by 2^-50.950.
*/
static double _log2p1 (double x)
{
  // int bug = x == 1.0;
  int bug = 0;
  /* for x > 0, 1+x is exact when 2^-29 <=  x < 2^53
     for x < 0, 1+x is exact when -1 < x <= 2^-30 */

  if (__builtin_expect (__builtin_fabs (x) < 0x1p-29, 0)) { // |x| < 2^-29
    /* |log2(1+x) - 1/log(2) * (x - x^2/2)| < 2^-59.584 * |log2(1+x)|
       (cf compoundf.sollya) */
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
  if (bug) printf ("u=%la\n", u);
  /* For x < 0x1p53, x + 1 is exact thus u = x+1.
     For x >= 2^53, we estimate log2(x) instead of log2(1+x),
     since log2(1+x) = log2(x) + log2(1+1/x),
     log2(x) >= 53 and |log2(1+1/x)| < 2^-52.471, the additional relative
     error is bounded by 2^-52.471/53 < 2^-58.198 */

  b64u64_u v = {.f = u};
  uint64_t m = v.u & 0xfffffffffffffull;
  int64_t e = (v.u >> 52) - 0x3ff + (m >= 0x6a09e667f3bcdull);
  if (bug) printf ("e=%ld\n", e);
  // 2^e/sqrt(2) < u < 2^e*sqrt(2)
  // v.u = ((0x3ffull + e) << 52) | m;
  v.u -= e << 52;
  double t = v.f;
  if (bug) printf ("t=%la\n", t);
  // u = 2^e*t with 1/sqrt(2) < t < sqrt(2)
  // thus log2(u) = e + log2(t)
  v.f += 2.0; // add 2 so that v.f is always in the binade [2, 4)
  int i = (v.u >> 46) - 0x10016; // 0 <= i <= 23
  if (bug) printf ("i=%d\n", i);
  // for 0<=i<24, inv[i] approximates 1/t for 1/2+(i+6)/32 <= t < 1/2+(i+7)/32
  static const double inv[] = {
    0x1.7p+0, 0x1.6p+0, 0x1.5p+0, 0x1.4p+0, 0x1.38p+0, 0x1.28p+0, 0x1.2p+0,
    0x1.18p+0, 0x1.1p+0, 0x1p+0, 0x1p+0, 0x1.e8p-1, 0x1.d8p-1, 0x1.dp-1,
    0x1.cp-1, 0x1.b4p-1, 0x1.a8p-1, 0x1.ap-1, 0x1.94p-1, 0x1.8cp-1, 0x1.8p-1,
    0x1.78p-1, 0x1.7p-1, 0x1.68p-1,
  };
  // log2inv[i] approximates -log2(inv[i])
  static const double log2inv[] = {
    -0x1.0c10500d63aa6p-1, -0x1.d6753e032ea0fp-2, -0x1.91bba891f1709p-2,
    -0x1.49a784bcd1b8bp-2, -0x1.24407ab0e073ap-2, -0x1.acf5e2db4ec94p-3,
    -0x1.5c01a39fbd688p-3, -0x1.08c588cda79e4p-3, -0x1.663f6fac91316p-4,
    0x0p+0, 0x0p+0, 0x1.1bb32a600549dp-4, 0x1.e0b1ae8f2fd56p-4,
    0x1.22dadc2ab3497p-3, 0x1.8a8980abfbd32p-3, 0x1.dac22d3e441d3p-3,
    0x1.169c05363f158p-2, 0x1.32bfee370ee68p-2, 0x1.5dfdcf1eeae0ep-2,
    0x1.7b89f02cf2aadp-2, 0x1.a8ff971810a5ep-2, 0x1.c819dc2d45fe4p-2,
    0x1.e7df5fe538ab3p-2, 0x1.042bd4b9a7c99p-1,
  };
  double r = inv[i];
  double z = __builtin_fma (r, t, -1.0); // exact, -1/32 <= z <= 17/512
  if (bug) printf ("r=%la z=%la\n", r, z);
  // we approximates log2(t) by -log2(r) + log2(r*t)
  double p = p1 (z); // approximates log2(r*t)
  return (double) e + (log2inv[i] + p);
}

// return the correct rounding of (1+x)^y for 0x1.7154759a0df53p-24 <= |t| <= 150
// otherwise -1.0
static double _exp2_1 (double t)
{
  // int bug = t == -0x1.fffffffa3aae2p+6;
  int bug = 0;
  double k = roundeven_finite (t); // 0 <= |k| <= 150
  if (bug) printf ("k=%la\n", k);
  double r = t - k; // |r| <= 1/2
  if (bug) printf ("r=%la\n", r);
  b64u64_u v = {.f = 3.015625 + r}; // 2.5 <= v <= 3.5015625
  // we add 2^-6 so that i is rounded to nearest
  int i = (v.u >> 46) - 0x10010; // 0 <= i <= 32
  if (bug) printf ("i=%d\n", i);
  // r is near (i-16)/2^5
  static const double exp2_T[] = { // exp2_T[i] = (i-16)/2^5
    -0x1p-1, -0x1.ep-2, -0x1.cp-2, -0x1.ap-2, -0x1.8p-2, -0x1.6p-2, -0x1.4p-2,
    -0x1.2p-2, -0x1p-2, -0x1.cp-3, -0x1.8p-3, -0x1.4p-3, -0x1p-3, -0x1.8p-4,
    -0x1p-4, -0x1p-5, 0x0p+0, 0x1p-5, 0x1p-4, 0x1.8p-4, 0x1p-3, 0x1.4p-3,
    0x1.8p-3, 0x1.cp-3, 0x1p-2, 0x1.2p-2, 0x1.4p-2, 0x1.6p-2, 0x1.8p-2,
    0x1.ap-2, 0x1.cp-2, 0x1.ep-2, 0x1p-1,
  };
  r -= exp2_T[i];
  if (bug) printf ("exp2_T[i]=%la r=%la\n", exp2_T[i], r);
  // now |r| <= 2^-6
  static const double exp2_U[] = {
    0x1.6a09e667f3bcdp-1, 0x1.71f75e8ec5f74p-1, 0x1.7a11473eb0187p-1,
    0x1.82589994cce13p-1, 0x1.8ace5422aa0dbp-1, 0x1.93737b0cdc5e5p-1,
    0x1.9c49182a3f09p-1, 0x1.a5503b23e255dp-1, 0x1.ae89f995ad3adp-1,
    0x1.b7f76f2fb5e47p-1, 0x1.c199bdd85529cp-1, 0x1.cb720dcef9069p-1,
    0x1.d5818dcfba487p-1, 0x1.dfc97337b9b5fp-1, 0x1.ea4afa2a490dap-1,
    0x1.f50765b6e454p-1, 0x1p+0, 0x1.059b0d3158574p+0, 0x1.0b5586cf9890fp+0,
    0x1.11301d0125b51p+0, 0x1.172b83c7d517bp+0, 0x1.1d4873168b9aap+0,
    0x1.2387a6e756238p+0, 0x1.29e9df51fdee1p+0, 0x1.306fe0a31b715p+0,
    0x1.371a7373aa9cbp+0, 0x1.3dea64c123422p+0, 0x1.44e086061892dp+0,
    0x1.4bfdad5362a27p+0, 0x1.5342b569d4f82p+0, 0x1.5ab07dd485429p+0,
    0x1.6247eb03a5585p+0, 0x1.6a09e667f3bcdp+0,
  };
  // 2^t = 2^k * exp2_U[i] * 2^r
  v.f = exp2_U[i] * q1 (r);
  double err = 0x1.36p-33; // 2^-33.225 * exp2_U[32] < err
  float lb = v.f - err, rb = v.f + err;
  if (lb != rb) return -1.0f;
  if (bug) printf ("q1(r)=%la\n", q1 (r));
  if (bug) printf ("v=%la k=%la\n", v.f, k);
  v.u += (int64_t) k << 52; // scale v by 2^k
  if (__builtin_expect (v.f < 0x1p-126f, 0)) { // underflow
    feraiseexcept (FE_UNDERFLOW);
  }
  if (bug) printf ("v=%la\n", v.f);
  return v.f;
}

float cr_compoundf (float x, float y)
{
  // int bug = x == 0x1.fffffep+127f && y == -0x1p+0f;
  int bug = 0;
  /* Rules from IEEE 754-2019 for compound (x, n) with n integer:
     (a) compound (x, 0) is 1 for x >= −1 or quiet NaN
     (b) compound (−1, n) is +Inf and signals the divideByZero exception for n < 0
     (c) compound (−1, n) is +0 for n > 0
     (d) compound (+/-0, n) is 1
     (e) compound (+Inf, n) is +Inf for n > 0
     (f) compound (+Inf, n) is +0 for n < 0
     (g) compound (x, n) is qNaN and signals the invalid exception for x < −1
     (h) compound (qNaN, n) is qNaN for n <> 0.
  */
  double xd = x, yd = y;
  b64u64_u tx = {.f = xd}, ty = {.f = yd};

  if(__builtin_expect (x == 0, 0)) // compound(0,y) = 1 except for y = sNaN
    return issignalingf (y) ? x + y : 1.0f;

  if(__builtin_expect (ty.u<<1 == 0, 0)) { // compound (x, 0)
    if (issignalingf (x)) return x + y; // x = sNaN
    return (x < -1.0f) ? 0.0f / 0.0f : 1.0f; // rules (g) and (a)
  }

  if(__builtin_expect ((ty.u<<1) >= (uint64_t)0x7ff<<53, 0)){ // y=Inf/NaN
    // the case x=0 was already checked above
    if((tx.u<<1) > (uint64_t)0x7ff<<53) return x + y; // x=NaN
    if((ty.u<<1) == (uint64_t)0x7ff<<53) { // y=+/-Inf
      int sy = ty.u>>63; // sign bit of y
      if (x < -1.0f) return 0.0f / 0.0f; // rule (g)
      if (x < 0 && sy == 0)
	return 0;
      if (0 < x && sy != 0)
        return 0;
      return __builtin_inf();
    }
    return x + y; // case y=NaN
  }

  uint64_t minus_one = 0xbff0000000000000ull; // encoding of -1
  if (__builtin_expect (tx.u >= minus_one || tx.u >> 52 == 0x7ff, 0)){
    // x is Inf, NaN or <= -1
    if ((tx.u<<1) == (uint64_t)0x7ff<<53){ // x is +Inf or -Inf
      if (x < -1.0f) return 0.0f / 0.0f; // x = -Inf, rule (g)
      // (1 + Inf)^y = +Inf for y > 0, +0 for y < 0
      return (ty.u>>63) ? 1.0f/x : x;
    }
    if ((tx.u<<1) > (uint64_t)0x7ff<<53) return x + y; // x is NaN
    if (x < -1.0f) return 0.0f / 0.0f; // x < -1.0: rule (g)
    // now x = -1
    if (ty.u>>63) { // y < 0
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE;
#endif
      return 1.0f/0.0f;
    } else // y > 0
      return 0.0f;
  }

  // now x > -1

  double l = _log2p1 (tx.f); // approximates log2(1+x)
  if (bug) printf ("log2(1+x)=%la\n", l);
  /* 2^-149 <= |l| < 128 */
  
  b64u64_u t = {.f = l * ty.f};
  if (bug) printf ("y*log2(1+x)=%la\n", t.f);
  /* since 2^-149 <= |l| < 128 and 2^-149 <= |y| < 2^128, we have
     2^-298 <= |t| < 2^135, thus no underflow/overflow in double is possible */

  // detect overflow/underflow
  if (__builtin_expect ((t.u << 1) >= (0x1018bull<<47), 0)) { // |t| >= 150
    if (t.u >> 63) { // t <= -256
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE; // underflow
#endif
      return 0x1p-126f * 0x1p-126f;
    } else { // t >= 256: overflow
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE; // overflow
#endif
      return 0x1p126f * 0x1p126f;
    }
  }

  // 2^t rounds to 1 to nearest when |t| <= 0x1.7154759a0df53p-24
  if (__builtin_expect ((t.u << 1) <= 0x7cee2a8eb341bea6ull, 0))
    return (t.u >> 63) ? 1.0f - 0x1p-25f : 1.0f + 0x1p-25f;

  float res = _exp2_1 (t.f);
  if (res != -1.0f)
    return res;

  if (bug) printf ("fast path failed\n");
  // fast path failed
  return res;
}
