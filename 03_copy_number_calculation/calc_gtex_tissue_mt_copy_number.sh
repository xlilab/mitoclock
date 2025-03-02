#!/usr/bin/env bash

set -e -u

OUT=mt_copy_number.txt

calc_copy_number(){
    local FILE=${1}
    local mt=`awk '{sum+=$3} END {print sum/NR}' ${FILE}_mt.depth`
    local chr=`awk '{sum+=$1} END {print sum/NR}' ${FILE}_chr.depth`
    if [ $(echo "$mt != 0" |bc) -eq 1 ] && [ $(echo "$chr != 0" |bc) -eq 1 ]
    then
        echo -ne $FILE "\t" 
        echo "scale=4; $mt*2/$chr" |bc 
    fi
}

#for FILE in `ls -l *mt.depth |grep "May  3" |awk '{print $9}' |sed 's/_mt.depth//g'`
export -f calc_copy_number
#parallel --jobs 10 --xapply --group calc_copy_number ::: $(for i in `ls *mt.depth |awk -F '_' '{print $1}' |uniq `; do echo $i ; done) ::: $OUT
for i in `ls *mt.depth |awk -F '_' '{print $1}' |uniq `; do echo $i ; done |parallel --jobs 16 --group calc_copy_number >$OUT