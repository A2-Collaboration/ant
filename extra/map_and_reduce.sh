#!/bin/bash -e

if [ $# -ne 2 ]; then
    echo "usage: ${0} <command> <outputfile>"
    echo "Pipe in file names"
    exit 1
fi

function tryget() { ${1} 2>/dev/null || echo ${2}; }

JOBS=$(tryget nproc 1)
COMMAND=$1
OUT=$2
PREFIX="__MAPREDUCE"
PARALLEL=$(tryget "which parallel" "")


WORK_DIR=`mktemp -d -p .`
echo $WORK_DIR

STDIN=$(cat)

NFILES=$(echo "$STDIN" | wc -l)
NFILES_PER_JOB=$(($NFILES/$JOBS+1))

echo "Splitting $NFILES files into ${JOBS} groups, $NFILES_PER_JOB each"

cd $WORK_DIR
echo "$STDIN" | split -d -a 2 -l $NFILES_PER_JOB - ${PREFIX}
JOBFILES=$(ls -1 $PREFIX*)
cd ..

echo "Creating chains..."
for i in $JOBFILES;
do
    Ant-chain --ignoretreeevents -o $WORK_DIR/chain_${i}.root $(cat $WORK_DIR/${i} | xargs)
done

echo "Running ${COMMAND} on all chains..."

if [ $PARALLEL ]; then
    echo $JOBFILES | tr ' ' '\n' | parallel --results mapreduce_logs $COMMAND -i $WORK_DIR/chain_{}.root -o $WORK_DIR/out_{}.root
else
    echo "GNU Parallel not found. Falling back to just backgrounding the jobs."
    for f in $JOBFILES; do
       $COMMAND -i $WORK_DIR/chain_${f}.root -o $WORK_DIR/out_${f}.root &
    done
    wait
fi

echo "Merging results..."
Ant-hadd $OUT $(for i in $JOBFILES; do echo $WORK_DIR/out_$i.root; done | xargs)

echo "Cleaning up..."
rm -rfv $WORK_DIR
