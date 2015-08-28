#!/bin/bash

set -e

CACHE=$HOME/cache
CACHE_REV=$CACHE/.rev2
NCPU=2

pushd $PWD

if [ -f $CACHE_REV ]; then
    export PATH=$CACHE/cmake/bin:$PATH
    export PATH=$CACHE/gcc/bin:$PATH
    export PATH=$CACHE/doxygen/bin:$PATH
    export LIBRARY_PATH=$CACHE/gcc/lib64:$LIBRARY_PATH
    export LD_LIBRARY_PATH=$CACHE/gcc/lib64:$LD_LIBRARY_PATH
    export CPLUS_INCLUDE_PATH=$CACHE/gcc/include/c++/4.8.2:$CPLUS_INCLUDE_PATH    
    cd $CACHE/root && source ./bin/thisroot.sh
else
    # Clean the cache
    echo "Cleaning the cache..."
    rm -rf $CACHE
    mkdir -p $CACHE
    
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

    # Get Doxygen
    wget http://ftp.stack.nl/pub/users/dimitri/doxygen-1.8.10.linux.bin.tar.gz -O $CACHE/doxygen.tar.gz
    mkdir $CACHE/doxygen && tar -xzf $CACHE/doxygen.tar.gz -C $CACHE/doxygen --strip-components=1
    rm $CACHE/doxygen.tar.gz
    export PATH=$CACHE/doxygen/bin:$PATH    
    
    ## APLCON
    rm -rf $CACHE/APLCON
    git clone https://github.com/A2-Collaboration-dev/APLCON $CACHE/APLCON
    cd $CACHE/APLCON && mkdir build && cd build && cmake .. && make -j$NCPU

    ## CERN ROOT
    rm -rf $CACHE/root
    wget http://root.cern.ch/download/root_v5.34.32.source.tar.gz -O $CACHE/root.tar.gz
    tar -xf $CACHE/root.tar.gz -C $CACHE
    cd $CACHE/root && ./configure --minimal && make -j$NCPU
    cd $CACHE/root && source ./bin/thisroot.sh
    rm $CACHE/root.tar.gz

    ## GSI HADES PLUTO
    rm -rf $CACHE/pluto
    wget http://web-docs.gsi.de/~hadeshyp/pluto/v5.42/pluto_v5.42.tar.gz -O $CACHE/pluto.tar.gz
    # pluto installs in some version dependent directory,
    # so account for this in tar command
    mkdir $CACHE/pluto && tar -xf $CACHE/pluto.tar.gz -C $CACHE/pluto --strip-components=1
    cd $CACHE/pluto && make -j$NCPU
    rm $CACHE/pluto.tar.gz

    touch $CACHE_REV
fi

popd

export APLCONSYS=$CACHE/APLCON
export PLUTOSYS=$CACHE/pluto

mkdir build && cd build
cmake -DCTEST_PARALLEL_JOBS=$NCPU -DCOVERAGE=On ..
make -j$NCPU
make build_and_test
