#!/bin/bash
set -euo pipefail

seqId=$1
sampleId=$2
pipeline=$3
version=$4
panel=$5

gatk=/share/apps/GATK-distros/GATK_4.0.4.0/gatk

$gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g" \
    Mutect2 \
    --reference /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    --input "$seqId"_"$sampleId".bam \
    --tumor $sampleId \
    --min-base-quality-score 10 \
    --genotype-germline-sites true \
    --genotyping-mode DISCOVERY \
    --intervals /home/cm/"$pipeline"/"$pipeline"-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    --interval-padding 100 \
    --max-population-af 0.5 \
    --output-mode EMIT_ALL_SITES \
    --germline-resource /data/db/human/gnomad/gnomad.exomes.r2.0.1.sites.vcf.gz \
    --af-of-alleles-not-in-resource 0.01 \
    --output "$seqId"_"$sampleId".vcf.gz \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --bamout "$seqId"_"$sampleId"_mutect_realigned.bam \
    --verbosity ERROR \
    --QUIET true
