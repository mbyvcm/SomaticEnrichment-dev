#!/bin.bash
set -euo pipefail

seqId="160629_M02641_0115_000000000-ARCF9"
panel="IlluminaTruSightCancer"
sampleId="14M08325"
ROI="./TEST/IlluminaTruSightCancer_ROI_b37.bed"

bash GATK4_somaticCNV.sh $seqId $sampleId $panel $ROI
