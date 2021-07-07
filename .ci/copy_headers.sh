#!/usr/bin/env bash

set -ex

# We need to put in place all relevant headers before running clang-tidy.
mkdir ../build
cd ../build
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -Dall=On -Dtesting=On -Dx11=Off -Dalien=Off \
      -Dcuda=Off -Dtmva-gpu=Off -Dveccore=Off ../root
# We need to prebuild a minimal set of targets which are responsible for header copy
# or generation.
make -j4 move_headers intrinsics_gen clang-tablegen-targets ClangDriverOptions \
         googletest Dictgen BaseTROOT BUILTIN_cppzmq
ln -s $PWD/compile_commands.json $PWD/../root/

