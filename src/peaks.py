import pandas as pd
import numpy as np
import os
import os.path
from tools import *
from neighborhoods import *

def make_candidate_regions_from_summits(macs_peaks, accessibility_file, genome_sizes, regions_whitelist, regions_blacklist, n_enhancers, peak_extend, outdir):
    ## Generate enhancer regions from MACS summits
    # 1. Count reads in DHS peaks
    # 2. Take top N regions, get summits, extend summits, merge

    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed")
    raw_counts_outfile = os.path.join(outdir, os.path.basename(macs_peaks) + "." + os.path.basename(accessibility_file) + ".Counts.bed")

    if regions_whitelist:
    	whitelist_command = "(bedtools intersect -a {regions_whitelist} -b {genome_sizes}.bed -wa | cut -f 1-3 && cat) |"
    else:
    	whitelist_command = ""

    if regions_blacklist:
    	blacklist_command = "bedtools intersect -v -wa -a stdin -b {regions_blacklist} | "
    else:
    	blacklist_command = ""

    #1. Count DHS/ATAC reads in candidate regions
    run_count_reads(accessibility_file, raw_counts_outfile, macs_peaks, genome_sizes, use_fast_count=True)

    #2. Take top N regions, get summits, extend summits, merge, remove blacklist, add whitelist, sort and merge
    #use -sorted in intersect command? Not worth it, both files are small
    command = "bedtools sort -i {raw_counts_outfile} -faidx {genome_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_enhancers} |" + \
        "bedtools intersect -b stdin -a {macs_peaks} -wa |" + \
        "awk '{{print $1 \"\\t\" $2 + $10 \"\\t\" $2 + $10}}' |" + \
        "bedtools slop -i stdin -b {peak_extend} -g {genome_sizes} |" + \
        "bedtools sort -i stdin -faidx {genome_sizes} |" + \
        "bedtools merge -i stdin | " + \
        blacklist_command + \
        "cut -f 1-3 | " + whitelist_command + \
        "bedtools sort -i stdin -faidx {genome_sizes} | bedtools merge -i stdin > {outfile}"

    command = command.format(**locals())

    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')

    return stdoutdata


def make_candidate_regions_from_peaks(macs_peaks, accessibility_file, genome_sizes, regions_whitelist, regions_blacklist, n_enhancers, peak_extend, outdir):
    ## Generate enhancer regions from MACS summits
    # 1. Count reads in DHS peaks
    # 2. Take top N regions, get summits, extend summits, merge

    outfile = os.path.join(outdir, os.path.basename(macs_peaks) + ".candidateRegions.bed")
    raw_counts_outfile = os.path.join(outdir, os.path.basename(macs_peaks) + os.path.basename(accessibility_file) + ".Counts.bed")

    if regions_whitelist:
        whitelist_command = "(bedtools intersect -a {regions_whitelist} -b {genome_sizes}.bed -wa | cut -f 1-3 && cat) |"
    else:
        whitelist_command = ""

    if regions_blacklist:
        blacklist_command = "bedtools intersect -v -wa -a stdin -b {regions_blacklist} | "
    else:
        blacklist_command = ""

    #1. Count DHS/ATAC reads in candidate regions
    run_count_reads(accessibility_file, raw_counts_outfile, macs_peaks, genome_sizes, use_fast_count=True)

    #2. Take top N regions, extend peaks (min size 500), merge, remove blacklist, add whitelist, sort and merge
    #use -sorted in intersect command? Not worth it, both files are small
    command = "bedtools sort -i {raw_counts_outfile} -faidx {genome_sizes} | bedtools merge -i stdin -c 4 -o max | sort -nr -k 4 | head -n {n_enhancers} |" + \
        "bedtools intersect -b stdin -a {macs_peaks} -wa |" + \
        "awk '{{ l=$3-$2; if (l < {peak_extend}) {{ $2 = $2 - int(({peak_extend}-l)/2); $3 = $3 + int(({peak_extend}-l)/2) }} print $1 \"\\t\" $2 \"\\t\" $3}}' |" + \
        "bedtools sort -i stdin -faidx {genome_sizes} |" + \
        "bedtools merge -i stdin | " + \
        blacklist_command + \
        "cut -f 1-3 | " + whitelist_command + \
        "bedtools sort -i stdin -faidx {genome_sizes} | bedtools merge -i stdin > {outfile}"

    command = command.format(**locals())

    p = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("Running: " + command)
    (stdoutdata, stderrdata) = p.communicate()
    err = str(stderrdata, 'utf-8')
    if not err == '':
        raise RuntimeError("Command failed.")

    return stdoutdata
