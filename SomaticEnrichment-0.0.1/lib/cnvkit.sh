#!/usr/bin/bash
set -euo pipefail


seqId=$1
sampleId=$2
panel=IlluminaTruSightCancer
ROI=./TEST/IlluminaTruSightCancer_ROI_b37.bed
FASTA=/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta

cnvkit=/share/apps/anaconda2/bin/cnvkit.py
bams=(`ls ./TEST/$seqId/IlluminaTruSightCancer/*M*/*.bam`)
samples=(`for i in ./TEST/$seqId/IlluminaTruSightCancer/*M*/; do basename $i;done`)

#$cnvkit access \
#    $FASTA \
#    -x ./resources/wgEncodeDacMapabilityConsensusExcludable.bed \
#    -o ./resources/access-excludes.hg19.bed

#$cnvkit autobin $bams -t $ROI -g ./resources/access-excludes.hg19.bed --annotate ./resources/refFlat.txt 

# iterate over bam files and calculate coverage
#for i in ./TEST/$seqId/IlluminaTruSightCancer/*M*/
#    do

#    sample=`basename $i`
#    echo $sample

#    $cnvkit coverage $i/*$sample.bam *.target.bed -o "$sample".targetcoverage.cnn
#    $cnvkit coverage $i/*$sample.bam *.antitarget.bed -o "$sample".antitargetcoverage.cnn

#    done        


for i in ${samples[@]}
    do

    test_sample=$i
    normal_samples=( ${samples[@]/$i} )
    
    tc="${normal_samples[@]/%/.targetcoverage.cnn}"
    atc="${normal_samples[@]/%/.antitargetcoverage.cnn}"

    echo "generating references"
    echo "TEST: "$test_sample
    echo "NORMALS: "${normal_samples[@]}
    

    $cnvkit reference $tc $atc --fasta $FASTA  -o "$test_sample"_reference.cnn
    echo "fixing ratios"
    $cnvkit fix "$test_sample".targetcoverage.cnn "$test_sample".antitargetcoverage.cnn "$test_sample"_reference.cnn -o "$test_sample".cnr
    echo "seqgmentation"
    $cnvkit segment "$test_sample".cnr -m cbs -o "$test_sample".cns
    echo "CNV calling"
    $cnvkit call "$test_sample".cns -o "$test_sample".call.cns

    $cnvkit metrics "$test_sample".targetcoverage.cnn "$test_sample".antitargetcoverage.cnn "$test_sample".cnr -s "$test_sample".cns > "$test_sample".metrics
    $cnvkit scatter "$test_sample".cnr -s "$test_sample".cns -o "$test_sample"-scatter.pdf
    $cnvkit breaks "$test_sample".cnr "$test_sample".cns > "$test_sample".breaks
    $cnvkit genemetrics "$test_sample".cnr -s "$test_sample".cns > "$test_sample".genemetrics
    $cnvkit sex "$test_sample".cnn "$test_sample".cnr "$test_sample".cns > "$test_sample".sex
    $cnvkit segmetrics "$test_sample".cnr -s "$test_sample".cns --ci --pi > "$test_sample".segmetrics


    done
