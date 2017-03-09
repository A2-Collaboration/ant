#!/bin/bash
set -e
NCPU=2
mkdir build && cd build
cmake -DCTEST_PARALLEL_JOBS=$NCPU -DAnt_COVERAGE=On ..
make -j$NCPU
make build_and_test
