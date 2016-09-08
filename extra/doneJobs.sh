#!/bin/bash -e

function die() {
    echo $* >&2
    exit 1
}

base=$1

test -d $base/log || die "log directory not found in ${base}"

find $base/log -name "*.log" | while read f; do
    fname=$(basename $f | sed s/log$/root/)
    rootname=$base/root/$fname

    if [ -f "${rootname}" ]; then
        echo $rootname
    else
        echo "No root file found for log file ${f}" >&2
    fi
done
