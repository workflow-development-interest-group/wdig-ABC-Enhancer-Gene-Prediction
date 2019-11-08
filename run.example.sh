# Test run ABC code
cd /seq/lincRNA/RAP/Promoters/ABC-Enhancer-Gene-Prediction;

## HiC Section
# python src/juicebox_dump.py \
# --hic_file https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined_30.hic \
# --juicebox "java -jar /seq/lincRNA/Software/juicer/GridEngine8/scripts/juicer_tools.jar" \
# --outdir example_chr22/input_data/HiC/raw/ \
# --chromosomes 22

# python src/compute_powerlaw_fit_from_hic.py \
# --hicDir example_chr22/input_data/HiC/raw/ \
# --outDir example_chr22/input_data/HiC/raw/powerlaw/ \
# --maxWindow 1000000 \
# --minWindow 5000 \
# --resolution 5000 \
# --chr 'chr22'

#Peaks
use MACS2;
macs2 callpeak \
-t example_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
-n wgEncodeUwDnaseK562AlnRep1.chr22.macs2 \
-f BAM \
-g hs \
-p .1 \
--call-summits \
--outdir example_chr22/ABC_output/Peaks/ 

unuse MACS2;
source /seq/lincRNA/thouis/VENV36/bin/activate;
use .samtools-1.8;

#Sort narrowPeak file
bedtools sort -faidx example_chr22/reference/chr22 -i example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak > example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted

python src/makeCandidateRegions.py \
--narrowPeak example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted \
--bam example_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam \
--outDir example_chr22/ABC_output/Peaks/ \
--chrom_sizes example_chr22/reference/chr22 \
--regions_blacklist reference/wgEncodeHg19ConsensusSignalArtifactRegions.bed \
--regions_whitelist example_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.TSS.500bp.chr22.bed \
--peakExtendFromSummit 250 \
--nStrongestPeaks 3000 

#Nbhds
python src/run.neighborhoods.py \
--candidate_enhancer_regions example_chr22/ABC_output/Peaks/wgEncodeUwDnaseK562AlnRep1.chr22.macs2_peaks.narrowPeak.sorted.candidateRegions.bed \
--genes example_chr22/reference/RefSeqCurated.170308.bed.CollapsedGeneBounds.chr22.bed \
--H3K27ac example_chr22/input_data/Chromatin/ENCFF384ZZM.chr22.bam \
--DHS example_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep1.chr22.bam,example_chr22/input_data/Chromatin/wgEncodeUwDnaseK562AlnRep2.chr22.bam \
--expression_table example_chr22/input_data/Expression/K562.ENCFF934YBO.TPM.txt \
--chrom_sizes example_chr22/reference/chr22 \
--ubiquitously_expressed_genes reference/UbiquitouslyExpressedGenesHG19.txt \
--cellType K562 \
--outdir example_chr22/ABC_output/Neighborhoods/ 

#Preds
python src/predict.py \
--enhancers example_chr22/ABC_output/Neighborhoods/EnhancerList.txt \
--genes example_chr22/ABC_output/Neighborhoods/GeneList.txt \
--HiCdir example_chr22/input_data/HiC/raw/ \
--hic_resolution 5000 \
--scale_hic_using_powerlaw \
--threshold .02 \
--cellType K562 \
--outdir example_chr22/ABC_output/Predictions/ \
--make_all_putative


