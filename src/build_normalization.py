import sys
import numpy as np
import pandas
import json


def compute_normalization(values, maxpercentile):
    #values = np.sort(values)[-count:]
    #assert len(values) == count
    # maxpercentile limits normalization to slightly trimmed values.
    # values outside the extremes will use linear fit to first or last pair
    return np.percentile(values, np.linspace(0, maxpercentile, 100)).tolist()


if __name__ == '__main__':
    enhancers = pandas.read_table(sys.argv[1])
    enhancers_tss = enhancers.loc[enhancers['isPromoterElement'] == True]
    enhancers_nontss = enhancers.loc[enhancers['isPromoterElement'] == False]

    normalizations = {}
    #normalizations['count'] = 15000
    normalizations['maxpercentile'] = 99.5
    normalizations['source'] = sys.argv[1]

    for col in enhancers.columns:
        #if 'RPM' in col or 'RPKM' in col:
        if col in ['DHS.RPM', 'ATAC.RPM','H3K27ac.RPM']:
            normalizations[col] = compute_normalization(enhancers[col].values,
                                                        normalizations['maxpercentile'])

            normalizations[col + '.PROMOTER'] = compute_normalization(enhancers_tss[col].values,
                                                        normalizations['maxpercentile'])

            normalizations[col + '.NON_PROMOTER'] = compute_normalization(enhancers_nontss[col].values,
                                                        normalizations['maxpercentile'])

    print(json.dumps(normalizations, indent=4))
