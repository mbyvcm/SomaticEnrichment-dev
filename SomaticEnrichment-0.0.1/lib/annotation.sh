#!/bin/bash
set -euo pipefail

seqId=$1
sampleId=$2

perl /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl \
    --input "$seqId"_"$sampleId"_filteredStr.vcf.gz \
    --format vcf \
    --output "$seqId"_"$sampleId"_filteredStr_annotated.vcf.gz \
    --vcf \
    --no_progress \
    --cache \
    --cache_version 86 \
    --force_overwrite \
    --no_stats \
    --offline \
    --dir /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --fasta /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --species homo_sapiens \
    --refseq
