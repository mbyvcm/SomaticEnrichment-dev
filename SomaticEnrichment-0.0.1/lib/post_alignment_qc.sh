#!/bin/bash
set -euo pipefail

seqId=$1
sampleId=$2
panel=$3
minimumCoverage=$4
ROI=$5

gatk=/share/apps/GATK-distros/GATK_4.0.4.0/gatk


#Convert BED to interval_list for later
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g \
    -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar BedToIntervalList \
    I=$ROI \
    O="$panel"_ROI.interval_list \
    SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict

#Alignment metrics: library sequence similarity
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g \
    -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar CollectAlignmentSummaryMetrics \
    R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    I="$seqId"_"$sampleId".bam \
    O="$seqId"_"$sampleId"_AlignmentSummaryMetrics.txt \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/state/partition1/tmpdir

#Calculate insert size: fragmentation performance
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g \
    -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar CollectInsertSizeMetrics \
    I="$seqId"_"$sampleId".bam \
    O="$seqId"_"$sampleId"_InsertMetrics.txt \
    H="$seqId"_"$sampleId"_InsertMetrics.pdf \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/state/partition1/tmpdir

#HsMetrics: capture & pooling performance
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g \
     -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar CollectHsMetrics \
     I="$seqId"_"$sampleId".bam \
     O="$seqId"_"$sampleId"_HsMetrics.txt \
     R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
     BAIT_INTERVALS="$panel"_ROI.interval_list \
     TARGET_INTERVALS="$panel"_ROI.interval_list \
     MAX_RECORDS_IN_RAM=2000000 \
     TMP_DIR=/state/partition1/tmpdir

#Generate per-base coverage: variant detection sensitivity
$gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g" \
    DepthOfCoverage \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -o "$seqId"_"$sampleId"_DepthOfCoverage \
    -I "$seqId"_"$sampleId".bam \
    -L $ROI \
    --countType COUNT_FRAGMENTS \
    --minMappingQuality 20 \
    --minBaseQuality 10 \
    -ct "$minimumCoverage" \
    --omitLocusTable \
    -rf MappingQualityUnavailable \
    -dt NONE

# zip per-base coverage file
awk -F'[\t|:]' '{if(NR>1) print $1"\t"$2"\t"$3}' "$seqId"_"$sampleId"_DepthOfCoverage | \
     /share/apps/htslib-distros/htslib-1.4.1/bgzip > "$seqId"_"$sampleId"_DepthOfCoverage.gz

# tabix index    
/share/apps/htslib-distros/htslib-1.4.1/tabix -b2 -e2 -s1 "$seqId"_"$sampleId"_DepthOfCoverage.gz

#add hotspot regions
cat /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/GermlineEnrichment_Hotspots.bed "$panel"_ClinicalCoverageTargets.bed | \
sort -k1,1V -k2,2n -k3,3n > "$panel"_ClinicalCoverageTargetsHotspots.bed

#Make PASS BED
/share/apps/htslib-distros/htslib-1.4.1/tabix -R "$panel"_ClinicalCoverageTargetsHotspots.bed \
"$seqId"_"$sampleId"_DepthOfCoverage.gz | \
awk -v minimumCoverage="$minimumCoverage" '$3 >= minimumCoverage { print $1"\t"$2-1"\t"$2 }' | \
sort -k1,1V -k2,2n -k3,3n | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge > "$seqId"_"$sampleId"_PASS.bed

#Calculate overlap between PASS BED and ClinicalCoverageTargets
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools coverage \
-a "$panel"_ClinicalCoverageTargetsHotspots.bed \
-b "$seqId"_"$sampleId"_PASS.bed | \
tee "$seqId"_"$sampleId"_ClinicalCoverageTargetMetrics.txt | \
awk '{pass[$4]+=$6; len[$4]+=$7} END { for(i in pass) printf "%s\t %.2f%\n", i, (pass[i]/len[i]) * 100 }' | \
sort -k1,1 > "$seqId"_"$sampleId"_ClinicalCoverageGeneCoverage.txt

#Make GAP BED
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools subtract \
-a "$panel"_ClinicalCoverageTargetsHotspots.bed \
-b "$seqId"_"$sampleId"_PASS.bed | \
sort -k1,1V -k2,2n -k3,3n \
> "$seqId"_"$sampleId"_Gaps.bed


#Gather QC metrics
meanInsertSize=$(head -n8 "$seqId"_"$sampleId"_InsertMetrics.txt | tail -n1 | cut -s -f5) #mean insert size
sdInsertSize=$(head -n8 "$seqId"_"$sampleId"_InsertMetrics.txt | tail -n1 | cut -s -f6) #insert size standard deviation
duplicationRate=$(head -n8 "$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt | tail -n1 | cut -s -f9) #The percentage of mapped sequence that is marked as duplicate.
totalReads=$(head -n8 "$seqId"_"$sampleId"_HsMetrics.txt | tail -n1 | cut -s -f6) #The total number of reads in the SAM or BAM file examine.
pctSelectedBases=$(head -n8 "$seqId"_"$sampleId"_HsMetrics.txt | tail -n1 | cut -s -f19) #On+Near Bait Bases / PF Bases Aligned.
totalTargetedUsableBases=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f2) #total number of usable bases. NB BQSR requires >= 100M, ideally >= 1B
meanOnTargetCoverage=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f3) #avg usable coverage
pctTargetBasesCt=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f7) #percentage panel covered with good enough data for variant detection
freemix=$(tail -n1 "$seqId"_"$sampleId"_Contamination.selfSM | cut -s -f7) #percentage DNA contamination. Should be <= 0.02
pctPfReadsAligned=$(grep ^PAIR "$seqId"_"$sampleId"_AlignmentSummaryMetrics.txt | awk '{print $7*100}') #Percentage mapped reads
atDropout=$(head -n8 "$seqId"_"$sampleId"_HsMetrics.txt | tail -n1 | cut -s -f51) #A measure of how undercovered <= 50% GC regions are relative to the mean
gcDropout=$(head -n8 "$seqId"_"$sampleId"_HsMetrics.txt | tail -n1 | cut -s -f52) #A measure of how undercovered >= 50% GC regions are relative to the mean

#gender analysis using Y chrom coverage
awk '{if ($1 == "Y") print $0}' /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed > Y.bed

if [ $(wc -l Y.bed |cut -d' ' -f1) -gt 0 ] && [ $(awk -v meanOnTargetCoverage="$meanOnTargetCoverage" 'BEGIN{printf "%3.0f", meanOnTargetCoverage}') -gt 10 ]; then

    #calc Y coverage
    /share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
    -T DepthOfCoverage \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -o "$seqId"_"$sampleId"_Y \
    --omitDepthOutputAtEachBase \
    --omitIntervalStatistics \
    --omitLocusTable \
    -L Y.bed \
    -XL Y:10000-2649520 \
    -XL Y:59034049-59363566 \
    -I "$seqId"_"$sampleId".bam \
    --countType COUNT_FRAGMENTS \
    --minMappingQuality 20 \
    -rf MappingQualityUnavailable \
    -dt NONE

    #extract Y mean coverage
    meanYCov=$(head -n2 "$seqId"_"$sampleId"_Y.sample_summary | tail -n1 | cut -s -f3)
    calcGender=$(awk -v meanOnTargetCoverage="$meanOnTargetCoverage" -v meanYCov="$meanYCov" 'BEGIN {if (meanYCov > 10 && (meanYCov / meanOnTargetCoverage) > 0.1){print "MALE"} else if (meanYCov < 10 && (meanYCov / meanOnTargetCoverage) $

    #clean up
    rm "$seqId"_"$sampleId"_Y.*

else
    calcGender="UNKNOWN"
fi

#Print QC metrics
echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tDuplicationRate\tPctSelectedBases\tPctTargetBasesCt\tMeanOnTargetCoverage\tGender\tEstimatedContamination\tMeanInsertSize\tSDInsertSize\tPercentMapped\tAtDropout\tGcDropout$
echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$duplicationRate\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage\t$calcGender\t$freemix\t$meanInsertSize\t$sdInsertSize\t$pctPfReadsAligned\t$atDropout\t$$

#print metaline for final VCF
echo \#\#SAMPLE\=\<ID\="$sampleId",Tissue\=Germline,WorklistId\="$worklistId",SeqId\="$seqId",Assay\="$panel",PipelineName\=GermlineEnrichment,PipelineVersion\="$version",RawSequenceQuality\="$rawSequenceQuality",PercentMapped\="$pctPfRe$

#Create PED file
#TSV Format: Family_ID, Individual_ID, Paternal_ID, Maternal_ID, Sex (1=male; 2=female; 0=unknown), Phenotype (Description or 1=unaffected, 2=affected, 0=missing). Missing data is 0
if [ ! -z ${familyId-} ]; then echo -ne "$familyId\t" > "$seqId"_"$sampleId"_pedigree.ped; else echo -ne "0\t" > "$seqId"_"$sampleId"_pedigree.ped; fi
echo -ne "$sampleId\t" >> "$seqId"_"$sampleId"_pedigree.ped
if [ ! -z ${paternalId-} ]; then echo -ne "$paternalId\t" >> "$seqId"_"$sampleId"_pedigree.ped; else echo -ne "0\t" >> "$seqId"_"$sampleId"_pedigree.ped; fi
if [ ! -z ${maternalId-} ]; then echo -ne "$maternalId\t" >> "$seqId"_"$sampleId"_pedigree.ped; else echo -ne "0\t" >> "$seqId"_"$sampleId"_pedigree.ped; fi
if [ ! -z ${sex-} ]; then echo -ne "$sex\t" >> "$seqId"_"$sampleId"_pedigree.ped; else echo -ne "0\t" >> "$seqId"_"$sampleId"_pedigree.ped; fi
if [ ! -z ${phenotype-} ]; then echo -e "$phenotype" >> "$seqId"_"$sampleId"_pedigree.ped; else echo -e "2" >> "$seqId"_"$sampleId"_pedigree.ped; fi

cat "$seqId"_"$sampleId"_pedigree.ped >> ../"$seqId"_pedigree.ped
