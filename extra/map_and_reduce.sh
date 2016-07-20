#!/bin/bash -e

if [ $# -ne 3 ]; then
    echo "usage: ${0} <Jobs> <command> <outputfile>"
    echo "Pipe in file names"
    exit 1
fi

JOBS=$1
COMMAND=$2
OUT=$3
PREFIX="__MAPREDUCE"
JOBFILES=$(eval echo ${PREFIX}{00..$(($JOBS -1))})

echo "Splitting files into ${JOBS} groups"
split -d -a 2 --number=r/${JOBS} - ${PREFIX}

echo "Creating chains..."
for i in $JOBFILES;
do
    Ant-chain -o chain_${i}.root $(cat ${i} | xargs)
done

echo "Running ${COMMAND} on all chains..."
echo $JOBFILES | tr ' ' '\n' | parallel --results mapreduce_logs $COMMAND -i chain_{}.root -o out_{}.root

echo "Merging results..."
echo Ant-hadd $OUT $(for i in $JOBFILES; do echo out_$i.root; done | xargs)
