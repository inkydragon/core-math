#!/bin/bash
apt-get update -qq && apt-get install -qq cmake ninja-build clang
cd /tmp
git clone https://github.com/llvm/llvm-project.git
cd llvm-project
mkdir build
cd build
cmake ../llvm -GNinja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DLLVM_ENABLE_PROJECTS="libc"
ninja libc
# library located at /tmp/llvm-project/build/projects/libc/lib/libllvmlibc.a
