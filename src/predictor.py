import numpy as np
from proximity import HiCFetcher, DistanceModel
import json
import pandas as pd
from tools import get_gene_name
import sys
import scipy.sparse as ssp
import time



# class Predictor(object):
#     def __init__(self, loaded_enhancers, **args):

#         self.cellType = args['cellType']

#         #HiC
#         self.hic_fetcher = HiCFetcher(args['HiCdir'], 
#                                         args['hic_gamma'],  
#                                         args['hic_gamma_reference'],
#                                         scale_with_powerlaw=args['scale_hic_using_powerlaw'],
#                                         tss_hic_contribution=args['tss_hic_contribution'])

#         #Power Law
#         if args['scale_hic_using_powerlaw']:
#             model_gamma = args['hic_gamma_reference']
#         else:
#             model_gamma = args['hic_gamma']
#         self.distance_model = DistanceModel(model_gamma)

#         #HiC adjustment parameters
#         self.tss_hic_contribution = args['tss_hic_contribution']
#         self.hic_pseudocount_distance = args['hic_pseudocount_distance']
#         self.hic_cap = 100 #args['hic_cap']

#     def estimate_contact_probability_from_distance(self, distance):
#         return self.distance_model(distance)

#     def predict_from_normalized_to_enhancers(self, enhancers, gene, window, tss_slop=500):
#         midpoint = (enhancers.start.values + enhancers.end.values) / 2
#         enhancers['distance'] = abs(gene['tss'] - midpoint)
        
#         #Add gene specific annotations
#         enhancers['isSelfPromoter'] = np.logical_and.reduce((enhancers.isPromoterElement == True, enhancers.start - tss_slop < gene.tss, enhancers.end + tss_slop > gene.tss))
#         # enhancers['isSelfGenic'] = np.logical_or(np.logical_and(enhancers.start > gene.start,enhancers.start < gene.end),
#         #                                                 np.logical_and(enhancers.end > gene.start, enhancers.end < gene.end))
#         enhancers['TargetGene'] = gene['name']
#         enhancers['TargetGeneTSS'] = gene['tss']
#         if 'is_ue' in gene:
#             enhancers['TargetGeneIsUbiquitouslyExpressed'] = gene['is_ue']

#         if 'Expression' in gene.index:
#             enhancers['TargetGeneExpression'] = gene['Expression'] 
#         else: 
#             enhancers['TargetGeneExpression'] = np.nan 

#         if 'PromoterActivityQuantile' in gene.index:
#             enhancers['TargetGenePromoterActivityQuantile'] = gene['PromoterActivityQuantile'] 
#         else: 
#             enhancers['TargetGenePromoterActivityQuantile'] = np.nan

#         is_self_tss = enhancers['isSelfPromoter'].values
#         if sum(is_self_tss) == 0:
#             print("No candidate self-promoter of {} {} {}. May want to investigate!".format(gene['name'], gene['chr'], gene['tss']))
#         elif sum(is_self_tss) > 1:
#             print("Found multiple self-promoters of {} {} {}. May want to investigate - but okay if candidate regions are not merged!".format(gene['name'], gene['chr'], gene['tss']))


#         #Get Hi-C data
#         hic_vals, rowmax, self.hic_exists, hic_vals_unscaled, rowmax_unscaled  = self.hic_fetcher(gene.chr, gene.tss, midpoint, enhancers)
#         enhancers['hic.distance'] = hic_vals
#         enhancers['hic.rowmax'] = rowmax
#         enhancers['hic.distance.unscaled'] = hic_vals_unscaled
#         enhancers['hic.rowmax.unscaled'] = rowmax_unscaled

#         #add hic pseudocount
#         enhancers['hic.distance.adj'], enhancers['hic_adjustment'] = self.add_hic_pseudocount(enhancers['distance'], 
#                                                                                                     enhancers['hic.distance'], 
#                                                                                                     enhancers['hic.rowmax'], 
#                                                                                                     hic_pseudocount_distance=self.hic_pseudocount_distance)

#         #Power Law
#         enhancers['estimatedCP'], cp_rowmax = self.estimate_contact_probability_from_distance(enhancers['distance'])
#         enhancers['estimatedCP.adj'] = self.normalize_proximity_contact_probability(enhancers['estimatedCP'], cp_rowmax)

#         #Compute ABC Score and related scores
#         enhancers = compute_score(enhancers, [enhancers['activity_base'], enhancers['hic.distance.adj']], "ABC")
#         enhancers = compute_score(enhancers, [enhancers['activity_base'], enhancers['estimatedCP.adj']], "powerlaw")
#         #enhancers = compute_score(enhancers, [enhancers['activity_base_noqnorm'], enhancers['hic.distance.adj']], "ABC.noqnorm")

#     def __call__(self, *args, **kwargs):
#         return self.predict(*args, **kwargs)

#     def add_hic_pseudocount(self, dists, hic_vals, hic_rowmax, hic_pseudocount_distance=1e6):
#         powerlaw_cp, cp_rowmax = self.estimate_contact_probability_from_distance(dists)
#         cp_at_distance = self.estimate_contact_probability_from_distance(hic_pseudocount_distance)[0]
#         adjustment = np.minimum(100 * (cp_at_distance / cp_rowmax), 100 * (powerlaw_cp / cp_rowmax))
#         return np.clip(100 * (hic_vals / hic_rowmax) + adjustment, 0, self.hic_cap), adjustment

#     def normalize_proximity_contact_probability(self, contact_probability, cp_rowmax):
#         return np.clip(100 * (contact_probability / cp_rowmax), 0, self.hic_cap)

#     # def chromosomes(self):
#     #     "returns a list of which chromosomes have data for predictions."
#     #     return list(set(self.hic_fetcher.chromosomes()))

#     def get_gene_prediction_stats(self, args, nearby_enhancers):
#         stats = pd.Series( {
#             'TargetGene' : nearby_enhancers.TargetGene.values[0],
#             'TargetChr' : nearby_enhancers.chr.values[0],
#             'TargetGeneTSS' : nearby_enhancers.TargetGeneTSS.values[0],
#             'nEnhancersConsidered' : int(nearby_enhancers.shape[0]),
#             'nDistalEnhancersPredicted' : int(sum(nearby_enhancers.loc[~nearby_enhancers['isPromoterElement'], args.score_column] > args.threshold))
#             })
#         return stats

def get_hic_file(chromosome, hic_dir):
    return '/seq/lincRNA/RAP/External/Rao2014-HiC/K562/5kb_resolution_intrachromosomal/' + chromosome + '/MAPQGE30/' + chromosome + '_5kb.KRobserved'

def make_predictions(chromosome, enhancers, genes, hic_file, args):
    pred = make_pred_table(chromosome, enhancers, genes, args)
    pred = add_hic(pred, hic_file, args)

    pred = compute_score(pred, [pred['activity_base'], pred['hic_kr_pl_scaled_adj']], "ABC")
    #pred = compute_score(pred, [pred['activity_base'], pred['estimatedCP.adj']], "powerlaw")

    return(pred)

def make_pred_table(chromosome, enh, genes, args):
    print('Making putative predictions table...')
    t = time.time()

    # import pdb
    # pdb.set_trace()

    enh['temp_merge_key'] = 0
    genes['temp_merge_key'] = 0

    enh = enh.loc[enh['chr'] == chromosome, :]
    genes = genes.loc[genes['chr'] == chromosome, :]

    genes = genes[['symbol','tss','Expression','PromoterActivityQuantile','isExpressed','temp_merge_key']]
    genes.columns = ['TargetGene', 'TargetGeneTSS', 'TargetGeneExpression', 'TargetGenePromoterActivityQuantile','TargetGeneIsExpressed','temp_merge_key']
    
    #Make cartesian product and then subset to EG pairs within window. 
    #TO DO: Replace with equivalent of bedtools intersect or GRanges overlaps 
    pred = pd.merge(enh, genes, on = 'temp_merge_key')

    pred['enh_midpoint'] = (pred['start'] + pred['end'])/2
    pred['distance'] = abs(pred['enh_midpoint'] - pred['TargetGeneTSS'])
    pred = pred.loc[pred['distance'] < args.window,:]

    print('Done. There are {} putative enhancers for chromosome {}'.format(pred.shape[0], chromosome))
    print('Elapsed time: {}'.format(time.time() - t))

    return pred
    
def add_hic(pred, hic_file, args):
    print('Begin HiC')
    #t = time.time()
    HiC = load_hic(hic_file, args)

    #Add hic
    #Could also do this by indexing into the sparse matrix (instead of merge) but this seems to be slower
    t = time.time()
    pred['enh_bin'] = np.floor(pred['enh_midpoint'] / args.hic_resolution).astype(int)
    pred['tss_bin'] = np.floor(pred['TargetGeneTSS'] / args.hic_resolution).astype(int)
    pred['bin1'] = np.amin(pred[['enh_bin', 'tss_bin']], axis = 1)
    pred['bin2'] = np.amax(pred[['enh_bin', 'tss_bin']], axis = 1)
    pred = pred.merge(HiC, how = 'left', on = ['bin1','bin2'])
    pred.fillna(value={'hic_kr' : 0}, inplace=True)
    
    #Index into sparse matrix
    #pred['hic_kr'] = [HiC[i,j] for (i,j) in pred[['enh_bin','tss_bin']].values.tolist()]
    
    print('HiC added to predictions table. Elapsed time: {}'.format(time.time() - t))

    # Add powerlaw scaling
    pred = scale_with_powerlaw(pred, args)

    #Add pseudocount
    pred = add_hic_pseudocount(pred, args)

    print("HiC Complete")
    #print('Elapsed time: {}'.format(time.time() - t))

    return(pred)

def load_hic(hic_file, args):
    print("Loading HiC")
    HiC_sparse_mat = hic_to_sparse(hic_file, args.window, args.hic_resolution)
    HiC = process_hic(HiC_sparse_mat, args)

    return(HiC)

def process_hic(hic_mat, args):
    #Make doubly stochastic.
    #Juicer produces a matrix with constant row/column sums. But sum is not 1 and is variable across chromosomes
    t = time.time()

    sums = hic_mat.sum(axis = 0)
    assert(np.max(sums)/np.min(sums[sums > 0]) < 1.001)
    mean_sum = np.mean(sums[sums > 0])
    print('HiC Matrix has row sums of {}, making double stochastic...'.format(mean_sum))
    kr_vec = np.repeat(np.sqrt(mean_sum), sums.shape[1])
    norm_mat = ssp.dia_matrix((1.0 / kr_vec, [0]), (sums.shape[1], sums.shape[1]))
    hic_mat = norm_mat * hic_mat * norm_mat

    #Adjust diagonal of matrix based on neighboring bins
    hic_mat = hic_mat.tolil(copy=False)
    for ii in range(1, sums.shape[1] - 1):
        hic_mat[ii,ii] = max(hic_mat[ii,ii-1], hic_mat[ii,ii+1]) * args.tss_hic_contribution / 100

    #Turn into dataframe
    hic_mat = hic_mat.tocoo(copy=False)
    hic_df = pd.DataFrame({'bin1': hic_mat.row, 'bin2': hic_mat.col, 'hic_kr': hic_mat.data})

    print('process.hic: Elapsed time: {}'.format(time.time() - t))


    return(hic_df)

def scale_with_powerlaw(pred, args):

    #TO DO: Include Hi-C scale here?
    if ~args.scale_hic_using_powerlaw:
        pred['hic_kr_pl_scaled'] = pred['hic_kr']
    else:
        dists = pred['distance'] / args.hic_resolution
        log_dists = np.log(dists + 1)
        powerlaw_fit = -1*args.hic_gamma * log_dists
        powerlaw_fit_reference = -1*args.hic_gamma_reference * log_dists
        pred['hic_kr_pl_scaled'] = pred['hic_kr'] * np.exp(powerlaw_fit_reference - powerlaw_fit)

    return(pred)

def add_hic_pseudocount(pred, args):

    #TO DO: Include Hi-C scale here - or deal with this constant issue. The pseudocount is too big
    dists = pred['distance'] / args.hic_resolution
    log_dists = np.log(dists + 1)
    powerlaw_fit = np.exp(-1*args.hic_gamma * log_dists)
    powerlaw_fit_at_ref = np.exp(-1*args.hic_gamma * np.log(args.hic_pseudocount_distance / args.hic_resolution + 1))
    
    pseudocount = np.amin(pd.DataFrame({'a' : powerlaw_fit, 'b' : powerlaw_fit_at_ref}), axis = 1)
    pred['hic_pseudocount'] = pseudocount
    pred['hic_kr_pl_scaled_adj'] = pred['hic_kr_pl_scaled'] + pseudocount

    return(pred)

def compute_score(enhancers, product_terms, prefix):

    scores = np.column_stack(product_terms).prod(axis = 1)
    #total = sum(scores)
    #normalized_scores = scores / (total if (total > 0) else 1)

    enhancers[prefix + '.Score.Numerator'] = scores
    enhancers[prefix + '.Score'] = enhancers[prefix + '.Score.Numerator'] / enhancers.groupby('TargetGene')[prefix + '.Score.Numerator'].transform('sum')

    #enhancers[prefix + '.Score'] = normalized_scores

    # #Why is this so slow...
    # def apply_division(df):
    #     df[prefix + '.Score'] = df[prefix + '.Score.Numerator'] / df[prefix + '.Score.Numerator'].sum()
    #     return(df)
    #enhancers = enhancers.groupby('TargetGene').apply(apply_division)

    return(enhancers)

def hic_to_sparse(filename, window, resolution):
    t = time.time()
    HiC = pd.read_table(filename, names=["bin1", "bin2", "hic_kr"],
                        header=None, engine='c', memory_map=True)

    # verify our assumptions
    assert np.all(HiC.bin1 <= HiC.bin2)

    # find largest entry
    max_pos = max(HiC.bin1.max(), HiC.bin2.max())
    hic_size = max_pos // resolution + 1

    # drop NaNs from hic
    HiC = HiC.loc[~np.isnan(HiC['hic_kr']),:]
    print("HiC has {} rows after dropping NaNs".format(HiC.shape[0]))

    #Needs to happen after KR norming step
    # window distance between contacts
    # too_far = abs(HiC.bin1 - HiC.bin2) > window
    # HiC = HiC.loc[~too_far,]
    # print("HiC has {} rows after windowing to {}".format(HiC.shape[0], window))

    # convert to sparse matrix in CSR (compressed sparse row) format, chopping
    # down to HiC bin size.  note that conversion to scipy sparse matrices
    # accumulates repeated indices, so this will do the right thing.
    row = np.floor(HiC.bin1.values / resolution).astype(int)
    col = np.floor(HiC.bin2.values / resolution).astype(int)
    dat = HiC.hic_kr.values
    # we want a symmetric matrix.  Easiest to do that during creation, but have to be careful of diagonal
    mask = (row != col)  # off-diagonal
    row2 = col[mask]  # note the row/col swap
    col2 = row[mask]
    dat2 = dat[mask]

    # concat and create
    row = np.hstack((row, row2))
    col = np.hstack((col, col2))
    dat = np.hstack((dat, dat2))

    print('hic.to.sparse: Elapsed time: {}'.format(time.time() - t))

    return ssp.csr_matrix((dat, (row, col)), (hic_size, hic_size))