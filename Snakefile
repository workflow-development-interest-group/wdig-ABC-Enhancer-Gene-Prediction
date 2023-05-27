rule all:
    input:
        "{output_dir}{output_name}_peaks.narrowPeak",
        "{output_dir}{output_name}_peaks.narrowPeak.candidateRegions.bed"

rule call_peaks:
    input:
        bam_file = config['bam_file']
    output:
        "{output_dir}{output_name}_peaks.narrowPeak"
    conda:
        "macs.yml"
    shell:
        """
        macs2 callpeak \
        -t {input.bam_file} \
        -n {config[output_name]} \
        -f BAM \
        -g hs \
        -p .1 \
        --call-summits \
        --outdir {config[output_dir]}
        """

rule sort_narrowPeak:
    input:
        "{output_dir}{output_name}_peaks.narrowPeak"
    output:
        "{output_dir}{output_name}_peaks.narrowPeak.sorted"
    conda:
        "macs.yml"
    shell:
        """
        bedtools sort -faidx {config[chrom_sizes_file]} -i {input} > {output}
        """

rule make_candidate_regions:
    input:
        narrowPeak = "{output_dir}{output_name}_peaks.narrowPeak.sorted",
        bam_file = config['bam_file']
    output:
        "{output_dir}{output_name}_peaks.narrowPeak.candidateRegions.bed"
    params:
        chrom_sizes = config['chrom_sizes_file'],
        blocklist = config['regions_blocklist_file'],
        includelist = config['regions_includelist_file'],
        peakExtendFromSummit = config['peakExtendFromSummit'],
        nStrongestPeaks = config['nStrongestPeaks']
    conda:
        "abcenv.yml"
    shell:
        """
        python src/makeCandidateRegions.py \
        --narrowPeak {input.narrowPeak} \
        --bam {input.bam_file} \
        --outDir {config[output_dir]} \
        --chrom_sizes {params.chrom_sizes} \
        --regions_blocklist {params.blocklist} \
        --regions_includelist {params.includelist} \
        --peakExtendFromSummit {params.peakExtendFromSummit} \
        --nStrongestPeaks {params.nStrongestPeaks}
        """
