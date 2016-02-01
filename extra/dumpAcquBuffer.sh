#!/bin/bash -e

BLEN=$((0x8000 * 10))

FILE="${1}"
shift

function dumpBuffer() {
    dd if="${FILE}" bs=$BLEN count=1 skip=$1
}

dumpBuffer 0

for buffer in $*; do    
    dumpBuffer $buffer
done

