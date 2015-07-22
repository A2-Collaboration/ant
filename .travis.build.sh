#!/bin/bash

set -e

CACHE=$HOME/cache

pushd $PWD

if [ -f $CACHE/.good ]; then
    export PATH=$CACHE/cmake/bin:$PATH
    export PATH=$CACHE/gcc/bin:$PATH
    export LIBRARY_PATH=$CACHE/gcc/lib64:$LIBRARY_PATH
    export LD_LIBRARY_PATH=$CACHE/gcc/lib64:$LD_LIBRARY_PATH
    export CPLUS_INCLUDE_PATH=$CACHE/gcc/include/c++/4.8.2:$CPLUS_INCLUDE_PATH    
    cd $CACHE/root && source ./bin/thisroot.sh
else
    # Get CMake 3.1
    wget https://github.com/Viq111/travis-container-packets/releases/download/cmake-3.1.2/cmake.tar.bz2 -O $CACHE/cmake.tar.bz2
    tar -xjf $CACHE/cmake.tar.bz2 -C $CACHE
    rm $CACHE/cmake.tar.bz2
    export PATH=$CACHE/cmake/bin:$PATH
    
    # Get GCC 4.8
    wget https://github.com/Viq111/travis-container-packets/releases/download/gcc-4.8.2/gcc.tar.bz2 -O $CACHE/gcc.tar.bz2
    tar -xjf $CACHE/gcc.tar.bz2 -C $CACHE
    rm $CACHE/gcc.tar.bz2
    export PATH=$CACHE/gcc/bin:$PATH
    export LIBRARY_PATH=$CACHE/gcc/lib64:$LIBRARY_PATH
    export LD_LIBRARY_PATH=$CACHE/gcc/lib64:$LD_LIBRARY_PATH
    export CPLUS_INCLUDE_PATH=$CACHE/gcc/include/c++/4.8.2:$CPLUS_INCLUDE_PATH

    ## APLCON
    rm -rf $CACHE/APLCON
    git clone https://github.com/A2-Collaboration-dev/APLCON $CACHE/APLCON
    cd $CACHE/APLCON && mkdir build && cd build && cmake .. && make -j2

    ## CERN ROOT
    rm -rf $CACHE/root
    wget http://root.cern.ch/download/root_v5.34.32.source.tar.gz -O $CACHE/root.tar.gz
    tar -xf $CACHE/root.tar.gz -C $CACHE
    cd $CACHE/root && ./configure --minimal && make -j2
    cd $CACHE/root && source ./bin/thisroot.sh
    rm $CACHE/root.tar.gz

    ## GSI HADES PLUTO
    rm -rf $CACHE/pluto
    wget http://web-docs.gsi.de/~hadeshyp/pluto/v5.42/pluto_v5.42.tar.gz -O $CACHE/pluto.tar.gz
    # pluto installs in some version dependent directory,
    # so account for this in tar command
    mkdir $CACHE/pluto && tar -xf $CACHE/pluto.tar.gz -C $CACHE/pluto --strip-components=1
    cd $CACHE/pluto && make -j2
    rm $CACHE/pluto.tar.gz

    touch $CACHE/.good
fi

popd

export APLCONSYS=$CACHE/APLCON
export PLUTOSYS=$CACHE/pluto

mkdir build && cd build
cmake -DCOVERALLS=On -DCMAKE_BUILD_TYPE=Debug ..
make -j2
make build_and_test
