#!/bin/bash
CORE_MATH_CHECK_STD=true LIBM=/tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a ./check.sh hypotf
