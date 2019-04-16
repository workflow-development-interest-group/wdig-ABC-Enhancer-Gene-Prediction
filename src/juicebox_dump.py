## Download and dump data from Aiden Lab juicebox, 
## matching the file formats and directory structure of the Rao et al. 2014 data from GEO
##  (so that we don't have to change any other code to read the data properly)
#
# use Java-1.8


import argparse
import subprocess

## Extracted from http://hicfiles.tc4ga.com/juicebox.properties
# hic_files = dict( 
#     THP1_Monocyte = "https://s3.amazonaws.com/hicfiles/external/phanstiel/updated_O/Snyder_O_30.hic",
#     THP1_Macrophage = "https://s3.amazonaws.com/hicfiles/external/phanstiel/A_inter_30.hic",
#     RPE1 = "https://hicfiles.s3.amazonaws.com/hiseq/rpe1/DarrowHuntley-2015/WT-combined_30.hic",
#     HCT116_ctrl ="https://s3.amazonaws.com/hicfiles/hiseq/degron/untreated/unsynchronized/combined_30.hic",
#     HCT116_cohesinKO = "https://s3.amazonaws.com/hicfiles/hiseq/degron/treated_6hr/unsynchronized/combined.hic",
#     CD34_HSPC = "https://s3.amazonaws.com/hicfiles/external/goodell/HSPC_30.hic",
#     CD3_TCell = "https://s3.amazonaws.com/hicfiles/external/goodell/tcell_30.hic",
#     ErythroidProgenitor = "https://s3.amazonaws.com/hicfiles/external/goodell/ep_30.hic",
#     H1hESC = "https://hicfiles.s3.amazonaws.com/external/dekker/4dn/h1hesc_30.hic",
#     GM12878_primary = "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/primary.hic",
#     GM12878_replicate = "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/replicate.hic",
#     mESC = "https://hicfiles.s3.amazonaws.com/hiseq/mes/primary_30.hic",
#     BJAB="/seq/lincRNA/RAP/KO-RNASeq/171108_ImmuneHiC/BJAB/mega/aligned/inter_30.hic",
#     Jurkat="/seq/lincRNA/RAP/KO-RNASeq/171108_ImmuneHiC/Jurkat/mega/aligned/inter_30.hic",
#     K562="https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic"
# )

def parseargs():
    # celltypes = list(hic_files.keys())
    # celltypes.sort()

    # readable = argparse.FileType('r')

    parser = argparse.ArgumentParser(description='Download and dump HiC data')
    parser.add_argument('--hic_file', required=True, help="Path or url to .hic file. Must be in Aiden lab format")
    parser.add_argument('--resolution', default=5000, help="Resolution of HiC to download. In units of bp.")
    parser.add_argument('--outdir', default=".")
    parser.add_argument('--juicebox', default="/seq/lincRNA/Software/juicer/GridEngine8/scripts-old/juicebox")
    parser.add_argument('--obskr', action="store_true", help="Only download the KR observed matrix (as opposed to the Raw matrix and the KR norm vector separately")
    parser.add_argument('--chromosomes', default="all", help="comma delimited list of chromosomes to download. ")

    return parser.parse_args()

def main(args):

    if args.chromosomes == "all":
        chromosomes = list(range(1,23)) + ['X']
    else:
        chromosomes = args.chromosomes.split(",")

    for chromosome in chromosomes:
        print("Starting chr" + str(chromosome) + " ... ")
        outdir = "{0}/{2}kb_resolution_intrachromosomal/chr{1}/".format(args.outdir, chromosome, int(args.resolution/1000))
        command = "mkdir -p " + outdir
        out = subprocess.getoutput(command)

        if args.obskr:
	        ## Download observed matrix with KR normalization
	        command = args.juicebox + " dump observed KR {0} {1} {1} BP {3} {2}/chr{1}_5kb.KRobserved".format(args.hic_file, chromosome, outdir, args.resolution)
	        print(command)
	        out = subprocess.getoutput(command)
        else:
	        ## Download raw observed matrix
	        command = args.juicebox + " dump observed NONE {0} {1} {1} BP {3} {2}/chr{1}_5kb.RAWobserved".format(args.hic_file, chromosome, outdir, args.resolution)
	        print(command)
	        out = subprocess.getoutput(command)
	        
	        ## Download KR norm file
	        command = args.juicebox + " dump norm KR {0} {1} {1} BP {3} {2}/chr{1}_5kb.KRnorm".format(args.hic_file, chromosome, outdir, args.resolution)
	        out = subprocess.getoutput(command)
	        print(out)

# def main(args):
#     processCellType(args)

if __name__ == '__main__':
    args = parseargs()
    main(args)
