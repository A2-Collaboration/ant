#!/bin/bash

set -e

CACHE=$HOME/cache

pushd $PWD

if [ -f $CACHE/.good ]; then
    cd $CACHE/root && source ./bin/thisroot.sh
else
    ## APLCON
    rm -rf $CACHE/APLCON
    git clone https://github.com/A2-Collaboration-dev/APLCON $CACHE/APLCON
    cd $CACHE/APLCON && mkdir build && cd build && cmake .. && make -j2

    ## CERN ROOT
    rm -rf $CACHE/root
    wget http://root.cern.ch/download/root_v5.34.32.source.tar.gz -O $CACHE/root.tar.gz
    tar xf $CACHE/root.tar.gz -C $CACHE
    cd $CACHE/root && ./configure --minimal && make -j2 && source ./bin/thisroot.sh

    ## GSI HADES PLUTO
    rm -rf $CACHE/pluto
    wget http://web-docs.gsi.de/~hadeshyp/pluto/v5.42/pluto_v5.42.tar.gz -O $CACHE/pluto.tar.gz
    # pluto installs in some version dependent directory, account for this in tar command
    mkdir $CACHE/pluto && tar xf $CACHE/pluto.tar.gz -C $CACHE/pluto --strip-components=1
    cd $CACHE/pluto && make -j2

    touch $CACHE/.good
fi

popd

export APLCONSYS=$CACHE/APLCON
export PLUTOSYS=$CACHE/pluto

mkdir build && cd build && cmake .. && make -j2 && make build_and_test
