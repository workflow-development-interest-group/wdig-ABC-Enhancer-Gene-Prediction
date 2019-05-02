import numpy as np
from proximity import HiCFetcher, DistanceModel
import json
import pandas as pd
from tools import get_gene_name
import sys


class Predictor(object):
    def __init__(self, loaded_enhancers, **args):

        self.cellType = args['cellType']

        #HiC
        self.hic_fetcher = HiCFetcher(args['HiCdir'], 
                                        args['hic_gamma'],  
                                        args['hic_gamma_reference'],
                                        scale_with_powerlaw=args['scale_hic_using_powerlaw'],
                                        tss_hic_contribution=args['tss_hic_contribution'])

        #Power Law
        if args['scale_hic_using_powerlaw']:
            model_gamma = args['hic_gamma_reference']
        else:
            model_gamma = args['hic_gamma']
        self.distance_model = DistanceModel(model_gamma)

        #Epi
        # self.DHS_column = normalized_dhs#args['DHS_column']
        # self.H3K27ac_column = "H3K27ac.RPM"

        #HiC adjustment parameters
        self.tss_hic_contribution = args['tss_hic_contribution']
        self.hic_pseudocount_distance = args['hic_pseudocount_distance']
        self.hic_cap = args['hic_cap']

        # # Quantile normalizations
        # if args['qnorm'] != '':
        #     normalizations = json.load(open(args['qnorm'], 'r'))
        #     #count = normalizations['count']
        #     maxpercentile = normalizations['maxpercentile']

        #     self.DHS_normalizer_promoter = make_normalizer(loaded_enhancers.ranges.loc[loaded_enhancers['isPromoterElement'] == True][self.DHS_column].values,
        #                                           normalizations[self.DHS_column + '.PROMOTER'],
        #                                           maxpercentile)
        #     self.DHS_normalizer_nonpromoter = make_normalizer(loaded_enhancers.ranges.loc[loaded_enhancers['isPromoterElement'] == False][self.DHS_column].values,
        #                                             normalizations[self.DHS_column + '.NON_PROMOTER'],
        #                                             maxpercentile)
        #     self.H3K27ac_normalizer_promoter = make_normalizer(loaded_enhancers.ranges.loc[loaded_enhancers['isPromoterElement'] == True][self.H3K27ac_column].values,
        #                                               normalizations[self.H3K27ac_column + '.PROMOTER'],
        #                                               maxpercentile)
        #     self.H3K27ac_normalizer_nonpromoter = make_normalizer(loaded_enhancers.ranges.loc[loaded_enhancers['isPromoterElement'] == False][self.H3K27ac_column].values,
        #                                                 normalizations[self.H3K27ac_column + '.NON_PROMOTER'],
        #                                                 maxpercentile)
        # else:
        #     self.DHS_normalizer_promoter = self.DHS_normalizer_nonpromoter = self.H3K27ac_normalizer_promoter = self.H3K27ac_normalizer_nonpromoter = lambda x: x

    def estimate_contact_probability_from_distance(self, distance):
        return self.distance_model(distance)

    def predict_from_normalized_to_enhancers(self, enhancers, gene, window, tss_slop=500):
        midpoint = (enhancers.start.values + enhancers.end.values) / 2
        enhancers['distance'] = abs(gene['tss'] - midpoint)
        
        #Add gene specific annotations
        enhancers['isSelfPromoter'] = np.logical_and.reduce((enhancers.isPromoterElement == True, enhancers.start - tss_slop < gene.tss, enhancers.end + tss_slop > gene.tss))
        # enhancers['isSelfGenic'] = np.logical_or(np.logical_and(enhancers.start > gene.start,enhancers.start < gene.end),
        #                                                 np.logical_and(enhancers.end > gene.start, enhancers.end < gene.end))
        enhancers['TargetGene'] = gene['name']
        enhancers['TargetGeneTSS'] = gene['tss']
        if 'is_ue' in gene:
            enhancers['TargetGeneIsUbiquitouslyExpressed'] = gene['is_ue']

        if 'Expression' in gene.index:
            enhancers['TargetGeneExpression'] = gene['Expression'] 
        else: 
            enhancers['TargetGeneExpression'] = np.nan 

        if 'PromoterActivityQuantile' in gene.index:
            enhancers['TargetGenePromoterActivityQuantile'] = gene['PromoterActivityQuantile'] 
        else: 
            enhancers['TargetGenePromoterActivityQuantile'] = np.nan

        is_self_tss = enhancers['isSelfPromoter'].values
        if sum(is_self_tss) == 0:
            print("No candidate element overlapping tss of {} {} {}. May want to investigate!".format(gene['name'], gene['chr'], gene['tss']))
        elif sum(is_self_tss) > 1:
            print("Found multiple elements overlapping tss of {} {} {}. May want to investigate - but okay if candidate regions are not merged!".format(gene['name'], gene['chr'], gene['tss']))


        #Get Hi-C data
        hic_vals, rowmax, self.hic_exists, hic_vals_unscaled, rowmax_unscaled  = self.hic_fetcher(gene.chr, gene.tss, midpoint, enhancers)
        enhancers['hic.distance'] = hic_vals
        enhancers['hic.rowmax'] = rowmax
        enhancers['hic.distance.unscaled'] = hic_vals_unscaled
        enhancers['hic.rowmax.unscaled'] = rowmax_unscaled

        #add hic pseudocount
        enhancers['hic.distance.adj'], enhancers['hic_adjustment'] = self.add_hic_pseudocount(enhancers['distance'], 
                                                                                                    enhancers['hic.distance'], 
                                                                                                    enhancers['hic.rowmax'], 
                                                                                                    hic_pseudocount_distance=self.hic_pseudocount_distance)

        #Activity
        #enhancers['activity_base'] = np.sqrt(enhancers['normalized_dhs'] * enhancers['normalized_h3k27ac'])
        #enhancers['activity_base_noqnorm'] = np.sqrt(enhancers[self.DHS_column] * enhancers['H3K27ac.RPM'])

        #Power Law
        enhancers['estimatedCP'], cp_rowmax = self.estimate_contact_probability_from_distance(enhancers['distance'])
        enhancers['estimatedCP.adj'] = self.normalize_proximity_contact_probability(enhancers['estimatedCP'], cp_rowmax)

        #Compute ABC Score and related scores
        enhancers = compute_score(enhancers, [enhancers['activity_base'], enhancers['hic.distance.adj']], "ABC")
        enhancers = compute_score(enhancers, [enhancers['activity_base'], enhancers['estimatedCP.adj']], "powerlaw")
        #enhancers = compute_score(enhancers, [enhancers['activity_base_noqnorm'], enhancers['hic.distance.adj']], "ABC.noqnorm")

    def __call__(self, *args, **kwargs):
        return self.predict(*args, **kwargs)

    def add_hic_pseudocount(self, dists, hic_vals, hic_rowmax, hic_pseudocount_distance=1e6):
        powerlaw_cp, cp_rowmax = self.estimate_contact_probability_from_distance(dists)
        cp_at_distance = self.estimate_contact_probability_from_distance(hic_pseudocount_distance)[0]
        adjustment = np.minimum(100 * (cp_at_distance / cp_rowmax), 100 * (powerlaw_cp / cp_rowmax))
        return np.clip(100 * (hic_vals / hic_rowmax) + adjustment, 0, self.hic_cap), adjustment

    def normalize_proximity_contact_probability(self, contact_probability, cp_rowmax):
        return np.clip(100 * (contact_probability / cp_rowmax), 0, self.hic_cap)

    def chromosomes(self):
        "returns a list of which chromosomes have data for predictions."
        return list(set(self.hic_fetcher.chromosomes()))

    def get_gene_prediction_stats(self, args, nearby_enhancers):
        stats = pd.Series( {
            'TargetGene' : nearby_enhancers.TargetGene.values[0],
            'TargetChr' : nearby_enhancers.chr.values[0],
            'TargetGeneTSS' : nearby_enhancers.TargetGeneTSS.values[0],
            'nEnhancersConsidered' : nearby_enhancers.shape[0],
            'nDistalEnhancersPredicted' : sum(nearby_enhancers.loc[~nearby_enhancers['isPromoterElement'], args.score_column] > args.threshold),
            'Total.Score' : sum(nearby_enhancers.loc[:,'ABC.Score'].apply(lambda x: x))
            })
        return stats

    # def add_normalized_data_to_enhancers(self, enhancers):
    #     is_promoter = enhancers['isPromoterElement'] == True
    #     enhancers.ranges.loc[is_promoter, 'normalized_dhs'] = self.DHS_normalizer_promoter(enhancers.ranges.loc[is_promoter][self.DHS_column])
    #     enhancers.ranges.loc[~is_promoter, 'normalized_dhs'] = self.DHS_normalizer_nonpromoter(enhancers.ranges.loc[~is_promoter][self.DHS_column])
    #     enhancers.ranges.loc[is_promoter, 'normalized_h3k27ac'] = self.H3K27ac_normalizer_promoter(enhancers.ranges.loc[is_promoter][self.H3K27ac_column])
    #     enhancers.ranges.loc[~is_promoter, 'normalized_h3k27ac'] = self.H3K27ac_normalizer_nonpromoter(enhancers.ranges.loc[~is_promoter][self.H3K27ac_column])


def compute_score(enhancers, product_terms, prefix):

    scores = np.column_stack(product_terms).prod(axis = 1)
    total = sum(scores)
    normalized_scores = scores / (total if (total > 0) else 1)

    enhancers[prefix + '.Score.Numerator'] = scores
    enhancers[prefix + '.Score'] = normalized_scores

    return(enhancers)

# def make_normalizer(values, target_vals, maxpercentile):
#     #assert len(values) >= count, "Need at least {} source values to build normalizer".format(count)
#     #values = np.sort(values)[-count:]
#     src_vals = np.percentile(values, np.linspace(0, maxpercentile, len(target_vals)))

#     # extend linear prediction to source = 0
#     # Seto slope_0 to 0 if we have a 0/0 situation
#     slope_0 = np.nan_to_num((target_vals[1] - target_vals[0]) / (src_vals[1] - src_vals[0]))
#     target_0 = target_vals[0] - slope_0 * src_vals[0]

#     # extend linear prediction to source = max(values)
#     slope_max = (target_vals[-1] - target_vals[-2]) / (src_vals[-1] - src_vals[-2])
#     target_max = target_vals[-1] + slope_max * (max(values) - src_vals[-1])

#     src_vals = np.concatenate(([0], src_vals, [max(values)]))
#     target_vals = np.concatenate(([target_0], target_vals, [target_max]))

#     def normalizer(vals):
#         normed = np.interp(vals, src_vals, target_vals)
#         normed[normed < 0] = 0
#         return normed

#     return normalizer
