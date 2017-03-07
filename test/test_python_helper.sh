#!/bin/bash

TEST_DIR=$1
TOP_DIR=$2
CWD=$PWD

pushd $PWD

CMD=python
if [ ! -z "$3" ]; then
    CMD='coverage run --parallel-mode'
fi

cd $TOP_DIR
$CMD -m unittest discover test/$TEST_DIR
#[ -f .coverage.* ] && mv .coverage.* $CWD

popd

#[ ! -z "$3" ] && coverage combine
