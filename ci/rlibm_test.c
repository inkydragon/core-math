// to check RLIBM function log2:
// gcc -DSTR=log2 -O3 test.c -lmpfr -lm rlibm.a -fopenmp
// ./a.out

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fenv.h>
#include <math.h>
#include <mpfr.h>
#include <omp.h>

#define CAT1(X,Y) X ## Y
#define CAT2(X,Y) CAT1(X,Y)
#define FOO CAT2(STR,f)
#define RLIBM_FOO CAT2(rlibm_,FOO)
#define MPFR_FOO CAT2(mpfr_,STR)

double rlibm_log2f(float);

typedef union {float x; uint32_t u;} b32u32_t;

static float ref (float x, mpfr_rnd_t rnd) {
  mpfr_t y;
  mpfr_init2 (y, 24);
  mpfr_set_emin (-148);
  mpfr_set_flt (y, x, MPFR_RNDN);
  int inex = MPFR_FOO (y, y, rnd);
  mpfr_subnormalize (y, inex, rnd);
  float ret = mpfr_get_flt (y, rnd);
  mpfr_clear (y);
  return ret;
}

int rnd[] = {FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD};
int rnd2[] = {MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDZ};

int is_equal (float y, float z) {
  if (isnan (y)) return isnan (z);
  if (isnan (z)) return isnan (y);
  return y == z;
}

static void check (float x) {
  for (int r = 0; r < 4; r++) {
    // recommended way (cf README.md)
    fesetround (FE_TONEAREST);
    double temp  = RLIBM_FOO (x);
    fesetround (rnd[r]);
    float y = (float) temp;
    float z = ref (x, rnd2[r]);
    if (!is_equal (y, z)) {
      printf ("FAIL %s x=%a y=%a ref=%a\n",
              mpfr_print_rnd_mode (rnd2[r]), x, y, z);
      exit (1);
    }
  }
}

int main() {
#pragma omp parallel for schedule(dynamic,1024)
  for (uint64_t u = 0; u < (1ull<<32); u++) {
    b32u32_t v = {.u = u};
    check (v.x);
  }
}
