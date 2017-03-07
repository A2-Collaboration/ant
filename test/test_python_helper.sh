#!/bin/bash

TEST_DIR=$1
TOP_DIR=$2

pushd $PWD

CMD=python
if [ ! -z "$3" ]; then
    CMD='coverage run --parallel-mode'
fi

cd $TOP_DIR
$CMD -m unittest discover test/$TEST_DIR

popd

