#!/bin/bash

TEST_DIR=$1
TOP_DIR=$2
CWD=$PWD

pushd $PWD

cd $TOP_DIR
coverage run -m unittest discover test/$TEST_DIR
[ -f .coverage ] || mv .coverage $CWD

popd
