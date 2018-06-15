#!/bin/bash
set -euo pipefail

echo "performing base quality recalibraton"

seqId=$1
sampleId=$2
pipeline=$3
version=$4
panel=$5

gatk=/share/apps/GATK-distros/GATK_4.0.4.0/gatk

$gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g" \
    BaseRecalibrator \
    --reference /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    --known-sites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
    --known-sites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
    --known-sites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
    --input "$seqId"_"$sampleId"_rmdup.bam \
    --intervals /home/cm/"$pipeline"/"$pipeline"-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    --interval-padding 100 \
    --output "$seqId"_"$sampleId"_recal_data.table \
    --verbosity ERROR \
    --QUIET true

$gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g" \
    ApplyBQSR \
    --reference /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    --bqsr-recal-file "$seqId"_"$sampleId"_recal_data.table \
    --input "$seqId"_"$sampleId"_rmdup.bam \
    --output "$seqId"_"$sampleId".bam \
    --QUIET true \
    --verbosity ERROR

$gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g" \
    BaseRecalibrator \
    --reference /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    --known-sites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
    --known-sites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
    --known-sites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
    --input "$seqId"_"$sampleId".bam \
    --intervals /home/cm/"$pipeline"/"$pipeline"-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    --interval-padding 100 \
    --output "$seqId"_"$sampleId"_post_recal_data.table \
    --verbosity ERROR \
    --QUIET true
