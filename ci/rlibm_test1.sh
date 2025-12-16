#!/bin/bash
# updated 16 Dec 2025 (revision be69161 of RLIBM)
for f in log10 log2 log exp10 exp2 exp cosh sinh cospi sinpi sin cos; do
   echo Testing $f
   gcc -DSTR=$f -O3 -march=native ci/rlibm_test.c -lmpfr -lm /tmp/The-RLIBM-Project/libm/rlibm.a -fopenmp
   ./a.out
done
