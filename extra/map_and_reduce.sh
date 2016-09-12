#!/bin/bash -e

if [ $# -ne 2 ]; then
    echo "usage: ${0} <command> <outputfile>"
    echo "Pipe in file names"
    exit 1
fi

function tryget() { ${1} || echo ${2}; }

JOBS=$(tryget nproc 1)
COMMAND=$1
OUT=$2
PREFIX="__MAPREDUCE"
JOBFILES=$(eval echo ${PREFIX}{00..$(($JOBS -1))})

PARALLEL=$(tryget "wich parallel" "")

echo "Splitting files into ${JOBS} groups"
split -d -a 2 --number=r/${JOBS} - ${PREFIX}

echo "Creating chains..."
for i in $JOBFILES;
do
    Ant-chain --ignoretreeevents --macrooverwrite -o chain_${i}.root $(cat ${i} | xargs)
done

echo "Running ${COMMAND} on all chains..."

if [ $PARALLEL ]; then
    echo $JOBFILES | tr ' ' '\n' | parallel --results mapreduce_logs $COMMAND -i chain_{}.root -o out_{}.root
else
    echo "GNU Parallel not found. Falling back to just backgrounding the jobs."
    for f in $JOBFILES; do
       $COMMAND -i chain_${f}.root -o out_${f}.root &
    done
    wait
fi

echo "Merging results..."
Ant-hadd -f $OUT $(for i in $JOBFILES; do echo out_$i.root; done | xargs)

rm $JOBFILES
rm $(for i in $JOBFILES; do echo out_$i.root; done | xargs)
rm $(for i in $JOBFILES; do echo chain_$i.root; done | xargs)

