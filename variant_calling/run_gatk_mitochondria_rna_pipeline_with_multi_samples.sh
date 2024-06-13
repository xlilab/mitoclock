#!/usr/bin/env bash
#Description: 该脚本需放置于gatk-workflows目录下运行
#Author: Wang Zhenguo
#Created Time: 2021/08/31 11:16

set -e -u -o pipefail
ulimit -u 65535

echo "Usage: sh $0 [ bam_sample_name.list tissue ]"

INPUT=$1
TISSUE=$2
if [ ! -d json ]; then mkdir json; fi

mt_pipeline(){
    local NAME=$1
    local TISSUE=$2
    sed "s#\(/.\{2,\}/\).\{2,\}/GTEX-.\{2,\}\.Aligned#\1$TISSUE/$NAME\.Aligned#g" mitochondria_m2_RNA_wdl/MitochondriaPipeline.inputs.json > json/$NAME.json
    java -XX:ActiveProcessorCount=1 -XX:ConcGCThreads=1 -jar cromwell-52.jar run mitochondria_m2_RNA_wdl/MitochondriaPipeline.wdl --inputs json/$NAME.json
}

export -f mt_pipeline
parallel --jobs 3 --xapply mt_pipeline ::: $(cat $INPUT) ::: $TISSUE
