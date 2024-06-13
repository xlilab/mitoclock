tempfifo=$$.fifo
trap "exec 1000>&-;exec 1000<&-;exit 0" 2
mkfifo $tempfifo
exec 1000<>$tempfifo
rm -rf $tempfifo

for ((i=1; i<=25; i++))
do
    echo >&1000
done

#for i in `ls ~/data/mitochondrial_genome/gatk_mitochondria_pipeline/output |grep vcf$ |sed 's/.vcf//g' |grep GTEX |grep -vf ~/data/mitochondrial_genome/gatk_mitochondria_pipeline/temp`
for i in `ls mt_DNA_bam/whole_blood |awk -F '_' '{print $1}' |uniq`
do
    read -u1000
    {
        #if [[ ! `ls seq_depth/ |grep $i` ]]
        #then
        samtools depth -d 0 mt_DNA_bam/whole_blood/${i}_chrM.cram >seq_depth/${i}_mt.depth
        id=`echo $i |sed 's/-[0-9]\{4\}-SM-[0-9A-Z]\{5\}//g'`
        bcftools query -s $id -f '[%DP]\n' GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz |grep -v "\." >seq_depth/${i}_chr.depth
        #fi
        
        echo >&1000
    }&
done
