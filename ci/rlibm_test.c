// check RLIBM log2

#include <stdio.h>
#include <stdint.h>
#include <fenv.h>
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

static void check (float x) {
  for (int r = 0; r < 4; r++) {
    // recommended way (cf README.md)
    fesetround (FE_TONEAREST);
    double temp  = RLIBM_FOO (x);
    fesetround (rnd[r]);
    float y = (float) temp;
    float z = ref (x, rnd2[r]);
  }
}

int main() {
#pragma omp parallel for schedule(dynamic,1024)
  for (uint64_t u = 0; u < (1ull<<32); u++) {
    b32u32_t v = {.u = u};
    check (v.x);
  }
}
