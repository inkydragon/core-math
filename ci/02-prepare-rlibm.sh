#!/bin/bash
apt-get update -qq && apt-get install -qq cmake ninja-build clang
cp ci/rlibm_wrapper.c /tmp
cp ci/rlibm.patch /tmp
cd /tmp
git clone https://github.com/rutgers-apl/The-RLIBM-Project.git
cd The-RLIBM-Project/libm
cp /tmp/rlibm_wrapper.c wrapper.c
patch -i /tmp/rlibm.patch
make
# library located at /tmp/The-RLIBM-Project/libm/rlibm.a
