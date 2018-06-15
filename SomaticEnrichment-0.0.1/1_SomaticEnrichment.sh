#!/bin/bash
set -euo pipefail

# Description: Somatic Enrichment Pipeline
# Author:      Christopher Medway, All Wales Medical Genetics Service
# Mode:        BY_SAMPLE
# Use:         XXXX

version="0.0.1"

# load sample variables
. ./*.variables

# load pipeline variables
. /home/cm/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/"$panel"/"$panel".variables

# define fastq variables
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq)
    do
    
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    # cutadapt
#    ./lib/cutadapt.sh $fastqPair

    # fastqc
#    ./lib/fastqc.sh $seqId $sampleId $laneId

    # fastq to ubam
#    ./lib/fastq_to_ubam.sh $seqId $sampleId $laneId $worklistId $panel $expectedInsertSize

    # bwa
#    ./lib/bwa.sh $seqId $sampleId $laneId
    
    done

# merge & mark duplicate reads
#./lib/mark_duplicates.sh $seqId $sampleId 

# basequality recalibration
#if [ "$includeBQSR = true" ] ; then
#    echo "performing base quality recalibration"
#    ./lib/bqsr.sh $seqId $sampleId $pipelineName $version $panel
#else
#    echo "skipping base quality recalibration"
#    cp "$seqId"_"$sampleId"_rmdup.bam "$seqId"_"$sampleId".bam
#    cp "$seqId"_"$sampleId"_rmdup.bai "$seqId"_"$sampleId".bai
#fi


# variant calling
#./lib/mutect2.sh $seqId $sampleId $pipelineName $version $panel

# variant filter
#./lib/variant_filter.sh $seqId $sampleId

# annotation
./lib/annotation.sh $seqId $sampleId
