#!/bin/bash
apt-get update -qq && apt-get install -qq cmake ninja-build clang
cd /tmp
git clone https://github.com/rutgers-apl/The-RLIBM-Project.git
cd The-RLIBM-Project/libm
make
# library located at /tmp/The-RLIBM-Project/libm/rlibm.a
