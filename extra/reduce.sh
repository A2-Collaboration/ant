#~/bin/bash -e
#
# Usage example:
#  ls *.root | reduce.sh Ant-hadd out_%x.root %f
# always hadds 10 files and calls them out_0.root out_1.root ...
#

N=10
FILES=""
X=0

CMD=$*
echo $CMD

while read f; do
    FILES+=" ${f}"
    i=$((i+1))

    if [ $i -eq $N ]; then
        echo $FILES
        eval $(echo $CMD | sed s/%x/${X}/ | sed s/%f/"\${FILES}"/)
        FILES=""
        i=0
        X=$((X+1))
        
    fi
done
eval $(echo $CMD | sed s/%x/${X}/ | sed s/%f/"\${FILES}"/)
