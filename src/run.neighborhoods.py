import argparse
import os
from neighborhoods import *
from subprocess import getoutput

#TODO
#4. Look at JAVA bam file counting code
#5. Look at how peak file is loaded. Also change the name of this column

def parseargs(required_args=True):
    class formatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass

    epilog = ("")
    parser = argparse.ArgumentParser(description='Run neighborhood for a given cell type',
                                     epilog=epilog,
                                     formatter_class=formatter)
    readable = argparse.FileType('r')
    parser.add_argument('--cellType', required=required_args, help="Name of cell type")
    parser.add_argument('--params_file', required=required_args, help="File listing feature files for each cell type.")
    parser.add_argument('--outdir', required=required_args, help="Directory to write Neighborhood files to.")
    
    parser.add_argument('--genome', required=required_args, help="File listing genome annotations for each species/build")
    parser.add_argument('--gene_name_annotations', default="symbol,refseq", help="Comma delimited string of names corresponding to the gene identifiers present in the name field of the gene annotation bed file")
    parser.add_argument('--primary_gene_identifier', default="symbol", help="Primary identifier used to identify genes. Must be present in gene_name_annotations")

    parser.add_argument('--candidate_enhancer_regions', required=required_args, help="Bed file containing candidate_enhancer_regions")

    parser.add_argument('--tss_slop_for_class_assignment', default=500, type=int, help="Consider an element a promoter if it is within this many bp of a tss")
    parser.add_argument('--skip_rpkm_quantile', action="store_true", help="Do not compute RPKM and quantiles in EnhancerList")

    # replace textio wrapper returned by argparse with actual filename
    args = parser.parse_args()
    for name, val in vars(args).items():
        if hasattr(val, 'name'):
            setattr(args, name, val.name)
    print(args)
    return args

def processCellType(cellType, args):
    params = parse_params_file(cellType, args)

    params["outdir"] = args.outdir
    os.makedirs(params["outdir"], exist_ok=True)

    genome_params = pd.read_csv(args.genome, sep="\t").set_index("name").T.to_dict()
    genome = genome_params[params['genome_build']]

    #Setup Genes
    genes = load_genes(file = genome['genes'], 
                        ue_file = genome['ue_genes'], 
                        outdir = params["outdir"], 
                        expression_table_list = params["expression_table"], 
                        gene_id_names = args.gene_name_annotations, 
                        primary_id = args.primary_gene_identifier)
    genes = annotate_genes_with_features(genes = genes, 
                                            genome = genome, 
                                            **params)
    genes.to_csv(os.path.join(params["outdir"], "GeneList.txt"),
                 sep='\t', index=False, header=True, float_format="%.6f")

    #Setup Candidate Enhancers
    enhancers = load_enhancers(genes=genes, 
                                genome_sizes=genome['sizes'], 
                                candidate_peaks=args.candidate_enhancer_regions, 
                                skip_rpkm_quantile=args.skip_rpkm_quantile, 
                                cellType=cellType, 
                                tss_slop_for_class_assignment=args.tss_slop_for_class_assignment,
                                **params)
    enhancers.to_csv(os.path.join(params['outdir'], "EnhancerList.txt"),
                sep='\t', index=False, header=True, float_format="%.6f")
    enhancers[['chr', 'start', 'end', 'name']].to_csv(os.path.join(params['outdir'], "EnhancerList.bed"),
                sep='\t', index=False, header=False)

def main(args):
    processCellType(args.cellType, args)

if __name__ == '__main__':
    args = parseargs()
    main(args)