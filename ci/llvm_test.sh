#!/bin/bash
# https://libc.llvm.org/headers/math/index.html#higher-math-functions
# update 10 Dec 2025

if [ "$FORCE" != "" ]; then
   FLOAT_FUNS=(acos acosh asin asinh atan atan2 atanh cbrt cos cosh cospi erf exp exp10 exp10m1 exp2 exp2m1 expm1 hypot log log10 log1p log2 pow rsqrt sin sincos sinh sinpi tan tanh tanpi)
else
   # disable by default binary32 functions since exhaustive test is long
   FLOAT_FUNS=""
fi
for f in ${FLOAT_FUNS}; do
   echo Testing ${f}f
   CORE_MATH_CHECK_STD=true LIBM=/tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a ./check.sh ${f}f
done
# atan, atan2, pow: 1 ulp
DOUBLE_FUNS="acos asin atan atan2 cbrt cos exp exp10 exp2 expm1 hypot log log10 log1p log2 pow sin sincos tan"
for f in ${DOUBLE_FUNS}; do
   echo Testing $f
   CORE_MATH_CHECK_STD=true LIBM=/tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a ./check.sh $f
done
LDOUBLE_FUNS="sqrt"
for f in ${LDOUBLE_FUNS}; do
   echo Testing ${f}l
   CORE_MATH_CHECK_STD=true LIBM=/tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a ./check.sh ${f}l
done
if [ "$FORCE" != "" ]; then
   FLOAT16_FUNS="acos acosh acospi asin asinh asinpi atan atanh atanpi cos cosh cospi exp exp10 exp10m1 exp2 exp2m1 expm1 hypot log log10 log2 rsqrt sin sinh sinpi sqrt tan tanh tanpi"
else
   # disable by default hypotf16 test which is long
   FLOAT16_FUNS="acos acosh acospi asin asinh asinpi atan atanh atanpi cos cosh cospi exp exp10 exp10m1 exp2 exp2m1 expm1 log log10 log2 rsqrt sin sinh sinpi sqrt tan tanh tanpi"
fi
for f in ${FLOAT16_FUNS}; do
   echo Testing ${f}f16
   CORE_MATH_CHECK_STD=true LIBM=/tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a ./check.sh ${f}f16
done
# atan2: 1 ulp
FLOAT128_FUNS="atan2 sqrt"
for f in ${FLOAT128_FUNS}; do
   echo Testing ${f}f128
   CORE_MATH_CHECK_STD=true LIBM=/tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a ./check.sh ${f}f128
done
# log: ?
BFLOAT16_FUNS="log sqrt"
for f in ${BFLOAT16_FUNS}; do
   echo Testing ${f}bf16
   CORE_MATH_CHECK_STD=true LIBM=/tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a ./check.sh ${f}bf16
done
