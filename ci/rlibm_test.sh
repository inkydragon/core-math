#!/bin/bash
for f in log2 do;
   gcc -DSTR=f -O3 -march=native rlibm_test.c -lmpfr -lm /tmp/The-RLIBM-Project/libm/rlibm.a -fopenmp
   ./a.out
done
