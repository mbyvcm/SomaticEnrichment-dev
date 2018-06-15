#!/bin/bash

echo "running cutadapt.sh"

fastqPair=$1

# load sample variables
. ./*.variables

# load pipeline variables
. /home/cm/"$pipelineName"/"$pipelineName"-"$pipelineVersion"/"$panel"/"$panel".variables

laneId=$(echo "$fastqPair" | cut -d_ -f3)
read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

/share/apps/anaconda2/bin/cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -m 50 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    --quiet \
    "$read1Fastq" \
    "$read2Fastq"
