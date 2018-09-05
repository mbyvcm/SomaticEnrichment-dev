#!/bin/bash
set -euo pipefail

seqId=$1
sampleId=$2
panel=$3
minimumCoverage=$4
vendorCaptureBed=$5
padding=$6
minBQS=$7
minMQS=$8

gatk3=/share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar


# Generate per-base coverage: variant detection sensitivity
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar $gatk3 \
    -T DepthOfCoverage \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -I "$seqId"_"$sampleId".bam \
    -L $vendorCaptureBed \
    -o "$seqId"_"$sampleId"_DepthOfCoverage \
    --countType COUNT_FRAGMENTS \
    --minMappingQuality $minMQS \
    --minBaseQuality $minBQS \
    -ct $minimumCoverage \
    --omitLocusTable \
    -rf MappingQualityUnavailable \
    -dt NONE

# zip per-base coverage file
#awk -F'[\t|:]' '{if(NR>1) print $1"\t"$2"\t"$3}' "$seqId"_"$sampleId"_DepthOfCoverage | \
#     /share/apps/htslib-distros/htslib-1.4.1/bgzip > "$seqId"_"$sampleId"_DepthOfCoverage.gz

# tabix index
#/share/apps/htslib-distros/htslib-1.4.1/tabix -b2 -e2 -s1 "$seqId"_"$sampleId"_DepthOfCoverage.gz

#Make PASS BED
#/share/apps/htslib-distros/htslib-1.4.1/tabix -R "$panel"_ClinicalCoverageTargetsHotspots.bed \
#"$seqId"_"$sampleId"_DepthOfCoverage.gz | \
#awk -v minimumCoverage="$minimumCoverage" '$3 >= minimumCoverage { print $1"\t"$2-1"\t"$2 }' | \
#sort -k1,1V -k2,2n -k3,3n | \
#/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge > "$seqId"_"$sampleId"_PASS.bed

#Calculate overlap between PASS BED and ClinicalCoverageTargets
#/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools coverage \
#-a "$panel"_ClinicalCoverageTargetsHotspots.bed \
#-b "$seqId"_"$sampleId"_PASS.bed | \
#tee "$seqId"_"$sampleId"_ClinicalCoverageTargetMetrics.txt | \
#awk '{pass[$4]+=$6; len[$4]+=$7} END { for(i in pass) printf "%s\t %.2f%\n", i, (pass[i]/len[i]) * 100 }' | \
#sort -k1,1 > "$seqId"_"$sampleId"_ClinicalCoverageGeneCoverage.txt

#Make GAP BED
#/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools subtract \
#-a "$panel"_ClinicalCoverageTargetsHotspots.bed \
#-b "$seqId"_"$sampleId"_PASS.bed | \
#sort -k1,1V -k2,2n -k3,3n \
#> "$seqId"_"$sampleId"_Gaps.bed
