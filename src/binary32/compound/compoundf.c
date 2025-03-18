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

/* References:
   [1] FastTwoSum revisited, Claude-Pierre Jeannerod and Paul Zimmermann,
       Proceedings of Arith'2025, https://inria.hal.science/hal-04875749
*/

#include <stdio.h>
#include <stdint.h>
#include <errno.h>
#include <math.h> // needed to define compoundf since it is not in glibc
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

/* s + t <- a + b, assuming |a| >= |b|
   By Theorem 6 from [1], the error is bounded by 2u^2*ufp(s)
   where u = 2^-53, thus by 2^-105*ufp(s) <= 2^-105 |s|
*/
static void
fast_two_sum (double *s, double *t, double a, double b)
{
  *s = a + b;
  double e = *s - a;
  *t = b - e;
}

// Multiply exactly a and b, such that *hi + *lo = a * b.
static inline void a_mul(double *hi, double *lo, double a, double b) {
  *hi = a * b;
  *lo = __builtin_fma (a, b, -*hi);
}

// Multiply a double with a double double : a * (bh + bl)
static inline void s_mul (double *hi, double *lo, double a, double bh,
                          double bl) {
  a_mul (hi, lo, a, bh); /* exact */
  *lo = __builtin_fma (a, bl, *lo);
}

// Returns (ah + al) * (bh + bl) - (al * bl)
// We can ignore al * bl when assuming al <= ulp(ah) and bl <= ulp(bh)
static inline void d_mul(double *hi, double *lo, double ah, double al,
                         double bh, double bl) {
  double s, t;

  a_mul(hi, &s, ah, bh);
  t = __builtin_fma(al, bh, s);
  *lo = __builtin_fma(ah, bl, t);
}

// INVLOG2 = 1/log(2) * (1 + eps1) with |eps1| < 2^-55.976
#define INVLOG2 0x1.71547652b82fep+0 
// LOG2 = log(2) * (1 + eps) with |eps| < 2^-54.730
#define LOG2 0x1.62e42fefa39efp-1

/* for |z| <= 1/64, returns an approximation of log2(1+z)
   with relative error < 2^-49.642 (see analyze_p1() in compoundf.sage)
   and |p1(z)| < 2^-5.459 */
static double p1 (double z)
{
  // we include P[0] = 0 so that P[i] corresponds to degree i
  // this degree-8 polynomial generated by Sollya (cf p1.sollya)
  // has relative error < 2^-50.98
  static const double P[] = {0, 0x1.71547652b82fep0, -0x1.71547652b8d11p-1,
    0x1.ec709dc3a5014p-2, -0x1.715475b144983p-2, 0x1.2776c3fda300ep-2,
    -0x1.ec990162358cep-3, 0x1.a645337c29e27p-3};
    
  double z2 = z * z;
  double c5 = __builtin_fma (P[6], z, P[5]);
  double c3 = __builtin_fma (P[4], z, P[3]);
  double c1 = __builtin_fma (P[2], z, P[1]);
  double z4 = z2 * z2;
  c5 = __builtin_fma (P[7], z2, c5);
  c1 = __builtin_fma (c3, z2, c1);
  c1 = __builtin_fma (c5, z4, c1);
  return z * c1;
}

// put in h+l an approximation of log2(1+zh+zl)
static void p2 (double *h, double *l, double zh, double zl)
{
  /* degree-12 polynomial generated by Sollya which approximates
     log2(1+z) for |z| <= 1/64 with relative error < 2^-86.447
     (cf file p2.sollya)
  */
  static const double P2[] = {
    0x1.71547652b82fep+0, 0x1.777d0ffb8e4f8p-56,   // degree 1 (0, 1)
    -0x1.71547652b82fep-1, -0x1.777d0fd5e20dcp-57, // degree 2 (2, 3)
    0x1.ec709dc3a03fdp-2, 0x1.d28b18700358cp-56,   // degree 3 (4, 5)
    -0x1.71547652b82fep-2, -0x1.7802bd1628a8p-58,  // degree 4 (6, 7)
    0x1.2776c50ef9bfcp-2,                          // degree 5 (8)
    -0x1.ec709dc3a03f5p-3,                         // degree 6 (9)
    0x1.a61762a7c47f2p-3,                          // degree 7 (10)
    -0x1.71547652d5ef5p-3,                         // degree 8 (11)
    0x1.484b11004fbfap-3,                          // degree 9 (12)
    -0x1.2776c1df4182bp-3,                         // degree 10 (13)
    0x1.0cc4300023c2fp-3,                          // degree 11 (14)
    -0x1.ecc45ce2a71e9p-4,                         // degree 12 (15)
  };
  /* since we can't expect a relative accuracy better than 2^-86.447,
     the lower part of the double-double approximation only needs to
     have about 86-53 = 33 accurate bits. Since |p5*z^5/p1| < 2^-32,
     we evaluate terms of degree 5 or more in double precision only. */
  double t;
  *h = P2[12];
  for (int i = 11; i >= 5; i--)
    *h = __builtin_fma (*h, zh, P2[3+i]);
  // now evaluate in double-double
  s_mul (h, l, *h, zh, zl);
  fast_two_sum (h, &t, P2[6], *h);
  *l += t + P2[7];                 // h+l = h_old*z + P[6] + P[7]
  d_mul (h, l, *h, *l, zh, zl);
  fast_two_sum (h, &t, P2[4], *h);
  *l += t + P2[5];                 // h+l = (h_old+l_old)*z + P[4] + P[5]
  d_mul (h, l, *h, *l, zh, zl);
  fast_two_sum (h, &t, P2[2], *h);
  *l += t + P2[3];                 // h+l = (h_old+l_old)*z + P[2] + P[3]
  d_mul (h, l, *h, *l, zh, zl);
  fast_two_sum (h, &t, P2[0], *h);
  *l += t + P2[1];                 // h+l = (h_old+l_old)*z + P[0] + P[1]
  // final multiplication by z
  d_mul (h, l, *h, *l, zh, zl);
}

/* for |z| <= 2^-6, returns an approximation of 2^z
   with absolute error < 2^-43.540 (see analyze_q1() in compoundf.sage) */
static double q1 (double z)
{
  /* Q is a degree-4 polynomial generated by Sollya (cf q1.sollya)
     with absolute error < 2^-43.549 */
  static const double Q[] = {1.0, 0x1.62e42fef6d01ap-1,
             0x1.ebfbdff7feebap-3, 0x1.c6b167e579beep-5, 0x1.3b2b3428d06dep-7};
  double z2 = z * z;
  double c3 = __builtin_fma (Q[4], z, Q[3]);
  double c0 = __builtin_fma (Q[1], z, Q[0]);
  double c2 = __builtin_fma (c3  , z, Q[2]);
  return __builtin_fma (c2, z2, c0);
}

// for 0<=i<46, inv[i] approximates 1/t for 1/2+(i+13)/64 <= t < 1/2+(i+14)/64
static const double inv[] = {
  0x1.68p+0, 0x1.6p+0, 0x1.58p+0, 0x1.5p+0, 0x1.4cp+0, 0x1.44p+0, 0x1.4p+0,
  0x1.38p+0, 0x1.34p+0, 0x1.2cp+0, 0x1.28p+0, 0x1.2p+0, 0x1.1cp+0, 0x1.18p+0,
  0x1.14p+0, 0x1.1p+0, 0x1.0cp+0, 0x1.08p+0, 0x1p+0, 0x1p+0, 0x1.f4p-1,
  0x1.ecp-1, 0x1.e4p-1, 0x1.ep-1, 0x1.d8p-1, 0x1.dp-1, 0x1.cap-1, 0x1.c4p-1,
  0x1.bep-1, 0x1.b8p-1, 0x1.b2p-1, 0x1.acp-1, 0x1.a8p-1, 0x1.a2p-1,
  0x1.9cp-1, 0x1.98p-1, 0x1.92p-1, 0x1.8cp-1, 0x1.88p-1, 0x1.84p-1,
  0x1.8p-1, 0x1.7cp-1, 0x1.76p-1, 0x1.72p-1, 0x1.6ep-1, 0x1.6ap-1,
};

// log2inv[i] is a double-double approximation of -log2(inv[i])
static const double log2inv[][2] = {
  {-0x1.f7a8568cb06cfp-2,0x1.8f3673ffdd785p-57},
  {-0x1.d6753e032ea0fp-2,0x1.c141e66faaaadp-62},
  {-0x1.b47ebf73882a1p-2,0x1.6fae441c09d76p-56},
  {-0x1.91bba891f1709p-2,0x1.2d352bea51e59p-56},
  {-0x1.800a563161c54p-2,-0x1.9575b04fa6fbdp-57},
  {-0x1.5c01a39fbd688p-2,0x1.817fd3b7d7e5dp-56},
  {-0x1.49a784bcd1b8bp-2,0x1.b6d40900b2502p-62},
  {-0x1.24407ab0e073ap-2,0x1.f6e91ad16ecffp-56},
  {-0x1.11307dad30b76p-2,0x1.a7b47d2c352d9p-57},
  {-0x1.d49ee4c32597p-3,0x1.b85a54d7ee2fdp-58},
  {-0x1.acf5e2db4ec94p-3,0x1.01ee1343fe7cap-59},
  {-0x1.5c01a39fbd688p-3,0x1.817fd3b7d7e5dp-57},
  {-0x1.32ae9e278ae1ap-3,-0x1.f51f2c075a74cp-59},
  {-0x1.08c588cda79e4p-3,0x1.a7610e40bd6abp-57},
  {-0x1.bc84240adabbap-4,-0x1.8ecb169b9465fp-58},
  {-0x1.663f6fac91316p-4,-0x1.f3314e0985116p-58},
  {-0x1.0eb389fa29f9bp-4,0x1.30c22d15199b8p-58},
  {-0x1.6bad3758efd87p-5,-0x1.89b03784b5be1p-60},
  {0x0p+0,0x0p+0},
  {0x0p+0,0x0p+0},
  {0x1.184b8e4c56af8p-5,0x1.491f06c085bc2p-60},
  {0x1.d6ebd1f1febfep-5,0x1.155660710eb2ap-63},
  {0x1.4c560fe68af88p-4,0x1.c141e66faaaadp-61},
  {0x1.7d60496cfbb4cp-4,0x1.9ced1447e30adp-58},
  {0x1.e0b1ae8f2fd56p-4,0x1.92ce9636c90ap-58},
  {0x1.22dadc2ab3497p-3,-0x1.696e2866c718ep-58},
  {0x1.494f863b8df35p-3,-0x1.1562d61af73f8p-57},
  {0x1.7046031c79f85p-3,-0x1.0798d1aa21694p-57},
  {0x1.97c1cb13c7ec1p-3,-0x1.e95734abd2fccp-57},
  {0x1.bfc67a7fff4ccp-3,0x1.bc0af7b82e7d7p-61},
  {0x1.e857d3d361368p-3,-0x1.086fce864a1f6p-57},
  {0x1.08bce0d95fa38p-2,-0x1.3d56efe4338fep-58},
  {0x1.169c05363f158p-2,0x1.c8d43e017579bp-56},
  {0x1.2baa0c34be1ecp-2,-0x1.0132ae5e417cdp-58},
  {0x1.4106017c3eca3p-2,-0x1.c658d602e66bp-56},
  {0x1.4f6fbb2cec598p-2,0x1.e393a16b94b52p-56},
  {0x1.6552b49986277p-2,0x1.ac9080333c605p-56},
  {0x1.7b89f02cf2aadp-2,0x1.8f89e2eb553b2p-57},
  {0x1.8a8980abfbd32p-2,0x1.99aa6df8b7d83p-56},
  {0x1.99b072a96c6b2p-2,0x1.bca36fd02defp-56},
  {0x1.a8ff971810a5ep-2,0x1.817fd3b7d7e5dp-58},
  {0x1.b877c57b1b07p-2,-0x1.01d98c3531027p-58},
  {0x1.cffae611ad12bp-2,0x1.8a38b4175d665p-56},
  {0x1.dfdd89d586e2bp-2,0x1.38c8946414c6ap-59},
  {0x1.efec61b011f85p-2,0x1.6d261f1753e0bp-56},
  {0x1.0014332be0033p-1,-0x1.7398fe685f171p-55},
};

/* for -1 < x < 2^128, approximates log2(1+x)
   Note: assume x is representable in binary32.
   For |x| < 2^-29, the relative error is bounded by 2^-50.950,
   and for |x| >= 2^-29, it is bounded by 2^-49.623 (see analyze_log2p1()
   in compoundf.sage).
   Thus for all x, the relative error is bounded by 2^-49.623.
*/
static double _log2p1 (double x)
{
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
  /* For x < 0x1p53, x + 1 is exact thus u = x+1.
     For x >= 2^53, we estimate log2(x) instead of log2(1+x),
     since log2(1+x) = log2(x) + log2(1+1/x),
     log2(x) >= 53 and |log2(1+1/x)| < 2^-52.471, the additional relative
     error is bounded by 2^-52.471/53 < 2^-58.198 */

  b64u64_u v = {.f = u};
  uint64_t m = v.u & 0xfffffffffffffull;
  int64_t e = (v.u >> 52) - 0x3ff + (m >= 0x6a09e667f3bcdull);
  // 2^e/sqrt(2) < u < 2^e*sqrt(2), with -29 <= e <= 128
  v.u -= e << 52;
  double t = v.f;
  // u = 2^e*t with 1/sqrt(2) < t < sqrt(2)
  // thus log2(u) = e + log2(t)
  v.f += 2.0; // add 2 so that v.f is always in the binade [2, 4)
  int i = (v.u >> 45) - 0x2002d; // 0 <= i <= 45
  double r = inv[i];
  double z = __builtin_fma (r, t, -1.0); // exact, -1/64 <= z <= 1/64
  // we approximates log2(t) by -log2(r) + log2(r*t)
  double p = p1 (z);
  // p approximates log2(r*t) with rel. error < 2^-49.642, and |p| < 2^-5.459
  return (double) e + (log2inv[i][0] + p);
}

/* return the correct rounding of (1+x)^y, otherwise -1.0
   where t is an approximation of y*log2(1+x) with absolute error < 2^-42.139,
   assuming 0x1.7154759a0df53p-24 <= |t| <= 150 */
static double _exp2_1 (double t)
{
  double k = roundeven_finite (t); // 0 <= |k| <= 150
  double r = t - k; // |r| <= 1/2, exact
  b64u64_u v = {.f = 3.015625 + r}; // 2.5 <= v <= 3.5015625
  // we add 2^-6 so that i is rounded to nearest
  int i = (v.u >> 46) - 0x10010; // 0 <= i <= 32
  // r is near (i-16)/2^5
  static const double exp2_T[] = { // exp2_T[i] = (i-16)/2^5
    -0x1p-1, -0x1.ep-2, -0x1.cp-2, -0x1.ap-2, -0x1.8p-2, -0x1.6p-2, -0x1.4p-2,
    -0x1.2p-2, -0x1p-2, -0x1.cp-3, -0x1.8p-3, -0x1.4p-3, -0x1p-3, -0x1.8p-4,
    -0x1p-4, -0x1p-5, 0x0p+0, 0x1p-5, 0x1p-4, 0x1.8p-4, 0x1p-3, 0x1.4p-3,
    0x1.8p-3, 0x1.cp-3, 0x1p-2, 0x1.2p-2, 0x1.4p-2, 0x1.6p-2, 0x1.8p-2,
    0x1.ap-2, 0x1.cp-2, 0x1.ep-2, 0x1p-1,
  };
  r -= exp2_T[i]; // exact
  // now |r| <= 2^-6
  // exp2_U[i] approximates 2^exp2_T[i] with absolute error < 2^-53.092
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
  /* the absolute error on exp2_U[i] is bounded by 2^-53.092, with
     exp2_U[i] < 2^0.5, and that on q1(r) is bounded by 2^-43.540,
     with |q1(r)| < 1.011, thus |v| < 1.43, and the absolute error on v is
     bounded by ulp(v) + 2^0.5 * 2^-43.540 + 2^-53.092 * 1.011 < 2^-43.035.
     Now t approximates u := y*log2(1+x) with |t-u| < 2^-42.139 thus
     exp(u) = exp(t) * (1 + eps) with eps < 2^(2^-42.139)-1 < 2^-42.667.
     The total absolute error is thus bounded by 2^-43.035 + 2^-42.667
     < 2^-41.839. */
  b64u64_u err = {.f = 0x1.1fp-42}; // 2^-41.839 < 0x1.1fp-42
  v.u += (int64_t) k << 52; // scale v by 2^k
  err.u += (int64_t) k << 52; // scale the error by 2^k too
  float lb = v.f - err.f, rb = v.f + err.f;
  if (lb != rb) return -1.0f; // rounded test failed
#ifdef CORE_MATH_SUPPORT_ERRNO
  if (__builtin_expect (v.f < 0x1p-126, 0)) // underflow
    errno = ERANGE;
#endif

  return v.f;
}

// assume -1 < x < 2^128, and x is representable in binary32.
static void
log2p1_accurate (double *h, double *l, double x)
{
  if (1.0 >= x) // then 1.0 >= |x| since x > -1
    fast_two_sum (h, l, 1.0, x);
  else
    fast_two_sum (h, l, x, 1.0);
  // now 1 + x ~ h + l with error <= 2^-105 |h|

  b64u64_u v = {.f = *h}, w = {.f = *l};
  uint64_t m = v.u & 0xfffffffffffffull;
  int64_t e = (v.u >> 52) - 0x3ff + (m >= 0x6a09e667f3bcdull);
  // 2^e/sqrt(2) < h < 2^e*sqrt(2), with -29 <= e <= 128
  // divide h, l by 2^e
  v.u -= e << 52;
  w.u -= e << 52;
  *h = v.f;
  *l = w.f;

  // now 1 + x ~ 2^e * (h + l) thus log2(1+x) ~ e + log2(h+l)

  v.f = 2.0 + *h; // add 2 so that v.f is always in the binade [2, 4)
  int i = (v.u >> 45) - 0x2002d;
  double r = inv[i];
  double zh = __builtin_fma (r, *h, -1.0); // exact, -1/64 <= zh <= 1/64
  double zl = r * *l;
  // zh + zl = r*(h+l)-1
  // log2(h+l) = -log2(r) + log2(r*(h+l)) = -log2(r) + log2(1+zh+zl)
  double ph, pl;
  p2 (&ph, &pl, zh, zl);
  /* since |log2inv[i][0]| < 1 and e is integer, the precondition of
     fast_two_sum is fulfilled: either |e| >= 1, or e=0 and fast_two_sum
     is exact */
  fast_two_sum (h, l, (double) e, log2inv[i][0]);
  *l += log2inv[i][1];
  // e + log2(h+l) ~ rh + rl + ph + pl
  double t;
  fast_two_sum (h, &t, *h, ph);
  *l += t + pl;
}

static double
accurate_path (float x, float y)
{
  double h, l;

  log2p1_accurate (&h, &l, x);
  // h + l is a double-double approximation of log(1+x)
}

float cr_compoundf (float x, float y)
{
  /* Rules from IEEE 754-2019 for compound (x, n) with n integer:
     (a) compound (x, 0) is 1 for x >= -1 or quiet NaN
     (b) compound (-1, n) is +Inf and signals the divideByZero exception for n < 0
     (c) compound (-1, n) is +0 for n > 0
     (d) compound (+/-0, n) is 1
     (e) compound (+Inf, n) is +Inf for n > 0
     (f) compound (+Inf, n) is +0 for n < 0
     (g) compound (x, n) is qNaN and signals the invalid exception for x < -1
     (h) compound (qNaN, n) is qNaN for n <> 0.
  */
  double xd = x, yd = y;
  b64u64_u tx = {.f = xd}, ty = {.f = yd};

  if (__builtin_expect((tx.u & ty.u)<<1 == 0, 0)) {
    if (tx.u<<1 == 0)   // compound(0,y) = 1 except for y = sNaN
      return issignalingf (y) ? x + y : 1.0f;

    if (ty.u<<1 == 0) { // compound (x, 0)
      if (issignalingf (x)) return x + y; // x = sNaN
      return (x < -1.0f) ? 0.0f / 0.0f : 1.0f; // rules (g) and (a)
    }
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

  double l = _log2p1 (tx.f);
  /* l approximates log2(1+x) with relative error < 2^-49.623,
     and 2^-149 <= |l| < 128 */
  
  b64u64_u t = {.f = l * ty.f};
  /* since 2^-149 <= |l| < 128 and 2^-149 <= |y| < 2^128, we have
     2^-298 <= |t| < 2^135, thus no underflow/overflow in double is possible.
     The relative error is bounded by (1+2^-49.623)*(1+2^-52)-1 < 2^-49.368 */

  // detect overflow/underflow
  if (__builtin_expect ((t.u << 1) >= (0x406ull<<53), 0)) { // |t| >= 128
    if (t.u >= 0x3018bull<<46) { // t <= -150
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE; // underflow
#endif
      return 0x1p-126f * 0x1p-126f;
    } else if (!(t.u >> 63)) { // t >= 128: overflow
#ifdef CORE_MATH_SUPPORT_ERRNO
      errno = ERANGE; // overflow
#endif
      return 0x1p126f * 0x1p126f;
    }
  }

  /* since |t| < 150, the absolute error on t is bounded by
     150*2^-49.368 < 2^-42.139 */

  // 2^t rounds to 1 to nearest when |t| <= 0x1.715476ba97f14p-25
  if (__builtin_expect ((t.u << 1) <= 0x3e6715476ba97f14ull, 0))
    return (t.u >> 63) ? 1.0f - 0x1p-25f : 1.0f + 0x1p-25f;

  float res = _exp2_1 (t.f);
  if (res != -1.0f)
    return res;

  // fast path failed
  return res;
}

#ifndef SKIP_C_FUNC_REDEF // icx provides this function
/* just to compile since glibc does not contain this function */
float compoundf(float x, float y){
  return powf(1.0f+x,y);
}
#endif
