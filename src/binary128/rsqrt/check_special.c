/* Generate special cases for rsqrt testing.

Copyright (c) 2025 Alexei Sibidanov <sibid@uvic.ca>

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
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <math.h>
#include <unistd.h>
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#include <omp.h>
#endif

typedef union {
  unsigned __int128 a;
  __float128 f;
} b128u128_u;

int ref_fesetround(int);
void ref_init(void);

__float128 cr_rsqrtq(__float128);
__float128 ref_rsqrtq(__float128);

int rnd1[] = { FE_TONEAREST, FE_TOWARDZERO, FE_UPWARD, FE_DOWNWARD };

int rnd = 0;
int verbose = 0;

#define MAX_THREADS 192

static unsigned int Seed[MAX_THREADS];

/* define our own is_nan function to avoid depending from math.h */
static inline int is_nan(__float128 x){
  b128u128_u u = {.f = x};
  return (u.a<<1)>((unsigned __int128)0x7fff<<113);
}

static inline int is_equal(__float128 x, __float128 y){
  if(is_nan(x)) return is_nan(y);
  if(is_nan(y)) return is_nan(x);
  b128u128_u ux = {.f = x}, uy = {.f = y};
  return ux.a == uy.a;
}

static void check(__float128 x){
  __float128 y1 = ref_rsqrtq(x), y2 = cr_rsqrtq(x);
  if(!is_equal(y1, y2)) {
    printf("FAIL x=%Qa ref=%Qa z=%Qa\n", x, y1, y2);
    fflush(stdout);
  }
}

static __float128 get_random(int tid){
  b128u128_u v = {.a = 0};
  v.a |= (__int128)rand_r(Seed + tid) << 31*0;
  v.a |= (__int128)rand_r(Seed + tid) << 31*1;
  v.a |= (__int128)rand_r(Seed + tid) << 31*2;
  v.a |= (__int128)rand_r(Seed + tid) << 31*3;
  v.a |= (__int128)rand_r(Seed + tid) << 31*4;
  v.a &= ~((__int128)1<<127);
  return v.f;
}

int main(int argc, char *argv[]){
  while(argc >= 2){
    if(strcmp(argv[1], "--rndn") == 0){
      rnd = 0;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--rndz") == 0){
      rnd = 1;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--rndu") == 0){
      rnd = 2;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--rndd") == 0){
      rnd = 3;
      argc --;
      argv ++;
    } else if(strcmp(argv[1], "--verbose") == 0){
      verbose = 1;
      argc --;
      argv ++;
    } else {
      fprintf(stderr, "Error, unknown option %s\n", argv[1]);
      exit(1);
    }
  }
  
  ref_init();
  ref_fesetround(rnd);
  fesetround(rnd1[rnd]);
  
  printf("Checking random values\n");
#ifndef CORE_MATH_TESTS
#define CORE_MATH_TESTS 1000UL*1000*1000 /* total number of tests */
#endif
  
  unsigned int seed = getpid();
  for(int i = 0; i < MAX_THREADS; i++)
    Seed[i] = seed + i;
  
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
#pragma omp parallel for
#endif
  for(uint64_t n = 0; n < CORE_MATH_TESTS; n++){
    ref_init();
    ref_fesetround(rnd);
    int tid;
#if (defined(_OPENMP) && !defined(CORE_MATH_NO_OPENMP))
    tid = omp_get_thread_num();
#else
    tid = 0;
#endif
    __float128 x = get_random(tid);
    check(x);
  }
  return 0;
}
