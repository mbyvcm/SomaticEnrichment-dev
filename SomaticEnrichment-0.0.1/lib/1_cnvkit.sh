#!/bin/bash
#PBS -l walltime=20:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//"`)
cd $PBS_O_WORKDIR

cnvkit=$1
seqId=$2
panel=$3
sample=$4

# cd /data/results/$seqId/$panel/$sample

# mkdir -p CNVKit

$cnvkit coverage /data/results/$seqId/$panel/$sample/"$seqId"_"$sample".bam /data/results/$seqId/$panel/*.target.bed -o "$sample".targetcoverage.cnn
$cnvkit coverage /data/results/$seqId/$panel/$sample/"$seqId"_"$sample".bam /data/results/$seqId/$panel/*.antitarget.bed -o "$sample".antitargetcoverage.cnn

if [ -e /data/results/$seqId/$panel/"$sample".antitargetcoverage.cnn ]
then
    echo $sample >> /data/results/$seqId/$panel/samplesCNVKit.txt
fi
