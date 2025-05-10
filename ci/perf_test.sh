#!/bin/sh
export CORE_MATH_PERF_MODE=rdtsc CFLAGS="-march=x86-64-v2 -I/usr/local/include" LDFLAGS="-L/usr/local/lib"
./perf.sh expf
./perf.sh exp
./perf.sh expl
./perf.sh rsqrtq
