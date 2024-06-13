#!/usr/bin/env bash
#Description: 该脚本需放置于gatk-workflows目录下运行
#Author: Wang Zhenguo
#Created Time: 2020/04/15 20:24

set -e -u -o pipefail
ulimit -u 65535

echo "Usage: sh $0"

DIR=`pwd`/cromwell-executions/MitochondriaPipeline
if [ ! -d output_vcf ]; then mkdir output_vcf; fi
if [ ! -d output_coverage ]; then mkdir output_coverage; fi
OUT_vcf=`pwd`/output_vcf
OUT_coverage=`pwd`/output_coverage

rename_and_extract(){
    local NAME=undefined
    local FILE=$1
    local DIR=$2
    local OUT_vcf=$3
    local OUT_coverage=$4

    cd $DIR
    if [ -d "$FILE/call-SplitMultiAllelicSitesConsensus/" ]
    then
        cd $DIR/$FILE/call-SplitMultiAllelicSitesConsensus/execution
        local NAME=`ls *split.vcf |sed 's/.Aligned.sortedByCoord.out.patched.md.*vcf//'`
        cd $DIR
        rename $FILE $NAME $FILE
        ln $NAME/call-SplitMultiAllelicSitesOriginal/execution/*vcf $OUT_vcf/$NAME.original.vcf
        ln $NAME/call-SplitMultiAllelicSitesConsensus/execution/*vcf $OUT_vcf/$NAME.consensus.vcf
        ln $NAME/call-CoverageAtEveryBase/execution/per_base_coverage.tsv $OUT_coverage/$NAME.coverage.tsv
        echo "extract and rename $NAME"
    else
        rm -r $FILE
        echo "delete $FILE"
    fi
}

export -f rename_and_extract
cd $DIR
parallel rename_and_extract ::: $(for FILE in `ls |egrep "\w{8}-\w{4}-\w{4}-\w{4}-\w{12}"`; do echo $FILE; done) ::: $DIR ::: $OUT_vcf ::: $OUT_coverage

cd ../../
echo -e \
    "SampleID\tContamination\tSampleHomoplasmies\tSampleHeteroplasmies\tSampleMeanCoverage\tHgMajor\tHgQualityMajor\tHgMinor\tHgQualityMinor\tHomoplasmiesMajor\tHomoplasmiesMinor\tHeteroplasmiesMajor\tHeteroplasmiesMinor\tMeanHetLevelMajor\tMeanHetLevelMinor\tHG_Distance\tClusters" > haplocheck.txt
find $DIR -name output-data |xargs cat |sort |sed 's/\(GTEX-[0-9A-Z]*\)-[-0-9A-Z]*/\1/' >> haplocheck.txt

echo "done!"
