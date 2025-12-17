#!/bin/bash
# updated 16 Dec 2025 (revision be69161 of RLIBM)
# for f in log10 log2 log exp10 exp2 exp cosh sinh; do
L=/tmp/The-RLIBM-Project/libm/rlibm.a
for f in log2; do
   echo Testing $f
   CORE_MATH_CHECK_STD=true LIBM=$L EXTRA_CFLAGS="-DCORE_MATH_SUPPORT_ERRNO -DCORE_MATH_CHECK_INEXACT" ./check.sh ${f}f
done
