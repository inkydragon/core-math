#!/bin/bash
# Examples of usage:
# DRY=--dry ./ci.sh to only try compilation (of last modified functions)
# FORCE=true DRY=--dry ./ci.sh to only try compilation (of all functions)
# FORCE_FUNCTIONS="xxx yyy" ./ci.sh to force checking xxx and yyy
# CC=clang CFLAGS=-Werror ./ci.sh
# SKIP16=1 ./ci.sh to avoid _Float16 tests
# SKIPB16=1 ./ci.sh to avoid __bf16 tests

set -e # We want the above command to potentially fail, only set -e now.

if [ -z "$LAST_COMMIT" ]; then
    LAST_COMMIT="HEAD~"
    if [ -n "$CI_COMMIT_BEFORE_SHA" ] && [ "$CI_COMMIT_BEFORE_SHA" != "0000000000000000000000000000000000000000" ]; then
        LAST_COMMIT="$CI_COMMIT_BEFORE_SHA"
    fi
fi

# use the same order as on https://core-math.gitlabpages.inria.fr/
FUNCTIONS_EXHAUSTIVE=(acosf acosf16 acoshf acoshf16 acospif acospif16 asinf asinf16 asinhf asinhf16 asinpif atanf atanf16 atan2f16 atan2pif16 atanhf atanhf16 atanpif atanpif16 cbrtf cbrtf16 hcbrt compoundf16 cosf cosf16 coshf coshf16 cospif cospif16 erff erff16 erfcf erfcf16 expf expf16 hexp exp10f hexp10 exp10f16 exp10m1f exp10m1f16 exp2f hexp2 exp2f16 exp2m1f exp2m1f16 expm1f expm1f16 hypotf16 lgammaf lgammaf16 logf logf16 hlog log10f log10f16 hlog10 log10p1f log10p1f16 log1pf log1pf16 log2f log2f16 hlog2 log2p1f log2p1f16 powf16 rsqrtf rsqrtf16 hrsqrt sincosf sincosf16 sinf sinf16 sinhf sinhf16 sinpif sinpif16 sqrtf16 hsqrt tanf tanf16 tanhf tanhf16 tanpif tanpif16 tgammaf tgammaf16)
FUNCTIONS_WORST=(acos acosh acospi asin asinh asinpi atan atan2 atan2f atan2pi atan2pif atanh atanpi cbrt cbrtl cbrtq compoundf cos cosh cospi erf erfc exp expl expq exp10 exp10q exp10m1 exp2 exp2l exp2q exp2m1 expm1 expm1q hypot hypotf hypotl hypotq lgamma log log10 log10p1 log1p log2 log2l log2p1 pow powf powl rsqrt rsqrtl rsqrtq sin sincos sinh sinpi sqrtq tan tanh tanpi tgamma)
FUNCTIONS_SPECIAL=(acos acosf acosh acospi acospif asin asinh asinpi asinpif atan atanf atan2 atan2f atan2pi atan2pif atanh atanpi atanpif cbrt cbrtl compoundf cos cosh cospi cospif erf erfc erfcf exp expf expl expq exp10 exp10q exp10m1 exp2 exp2l exp2q exp2m1 exp2m1f expm1 expm1q hypot hypotf hypotl hypotq lgamma lgammaf log log10 log10p1 log1p log2 log2l log2p1 pow powf powl rsqrt rsqrtl rsqrtq sin sincos sinh sinpi tan tanh tanpi tanpif tgamma)

echo "Reference commit is $LAST_COMMIT"

check () {
    KIND="$1"
    if [ "$FORCE" != "" ]; then
        doit=1
    elif ! { echo "$FORCE_FUNCTIONS" | tr ' ' '\n' | grep --quiet '^'"$FUNCTION"'$'; } && git diff --quiet "$LAST_COMMIT".. -- src/*/*/$FUNCTION.c; then
        doit=0
    else
        doit=1
    fi
    if [ "$doit" == "1" ] && [ "$SKIP128" == "1" ] && $CC -E src/*/*/$FUNCTION.c | grep -q  __int128; then
        echo "__int128 support is needed for" $FUNCTION "but is not available"
        doit=0
    fi
    if [ "$doit" == "1" ] && [ "$SKIP80" == "1" ] && echo src/*/*/$FUNCTION.c | grep -q binary80; then
        echo "binary80 support is needed for" $FUNCTION "but is not available"
        doit=0
    fi
    if [ "$doit" == "1" ] && [ "$SKIP16" == "1" ] && echo src/*/*/$FUNCTION.c | grep -q binary16; then
        echo "With SKIP16, skipping " $FUNCTION
        doit=0
    fi
    if [ "$doit" == "1" ] && [ "$SKIPB16" == "1" ] && echo src/*/*/$FUNCTION.c | grep -q binaryb16; then
        echo "With SKIPB16, skipping " $FUNCTION
        doit=0
    fi
    if [ "$doit" == "1" ] && [ "$SKIPQ" == "1" ] && [ "`basename $FUNCTION q`" != "$FUNCTION" ]; then
        echo "libquadmath is needed for" $FUNCTION "but is not available"
        doit=0
    fi

    if [ "$doit" == "0" ]; then
        echo "Skip $FUNCTION"
    else
        echo "Checking $FUNCTION..."				
	# we want to detect compiler warnings
	if [ "$EXTRA_CFLAGS" != "" ]; then
	    EXTRA_CFLAGS=-Werror
	fi
        EXTRA_CFLAGS=$EXTRA_CFLAGS ./check.sh $DRY "$KIND" "$FUNCTION"
    fi
}

if [ -z "$CC" ]; then
	CC="cc"
fi

if $CC -E $CFLAGS ci/int128test.c -o /dev/null &> /dev/null; then
    echo "Compiler supports __int128"
else
    echo "Compiler lacks __int128 support"
    SKIP128=1
fi

if $CC -E $CFLAGS ci/ldbl80test.c -o /dev/null &> /dev/null; then
    echo "long double is binary80"
else
    echo "long double is not binary80"
    SKIP80=1
fi

if $CC -c $CFLAGS ci/quadmath_test.c -o /dev/null &> /dev/null; then
   echo "Compiler supports libquadmath"
else
   echo "Compiler lacks libquadmath support"
   SKIPQ=1
fi

for FUNCTION in "${FUNCTIONS_EXHAUSTIVE[@]}"; do
    check --exhaustive
done

for FUNCTION in "${FUNCTIONS_WORST[@]}"; do
    check --worst
done

for FUNCTION in "${FUNCTIONS_SPECIAL[@]}"; do
    check --special
done
