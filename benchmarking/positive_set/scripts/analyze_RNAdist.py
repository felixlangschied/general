import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import json
import scipy.stats as sts


res_path = '/home/felixl/PycharmProjects/general/def_benchmark/data/ncortho_wDistance.tsv'
fp_path = '/home/felixl/PycharmProjects/general/def_benchmark/data/pot_FP_from_blast_folds.json'


# read results
resdf = pd.read_csv(res_path, sep='\t')
resdf['mirna'] = [mirna.replace('_pre', '') for mirna in resdf['mirna']]
resdf = resdf[resdf['species'] != 'Homo_sapiens']
resdf['distance'] = pd.to_numeric(resdf['distance'])
a = np.char.array(resdf['species'].values)
b = np.char.array(resdf['mirna'].values)
resdf['specmir'] = (a + b'|' + b).astype(str)
print(resdf)
# find interquartile range of distances for cutoff
dists = resdf['distance'].to_numpy()
dists = dists[~np.isnan(dists)]
# dist_med = np.median(dists)
dist_iqr = np.percentile(dists, [25, 50, 75, 95, 98, 100])
# print(f'median: {dist_med}')
print(f'iqr: {dist_iqr}')
# plot distance distribution
sns.violinplot(x=resdf['distance'])
plt.tight_layout()
plt.show()

# read potential FPs
with open(fp_path, 'r') as fph:
    fpdf = pd.DataFrame.from_dict(json.load(fph))
    fpdf['mirna'] = [mirna.replace('Hsa-', '') for mirna in fpdf['mirna']]
    a = np.char.array(fpdf['species'].values)
    b = np.char.array(fpdf['mirna'].values)
    fpdf['specmir'] = (a + b'|' + b).astype(str)
    fpspecmir = fpdf['specmir'].to_list()

# find distances of potential FPs
fp_dist = resdf[resdf['specmir'].isin(fpspecmir)]
fp_dist = fp_dist.filter(['species', 'mirna', 'score', 'distance', 'scheme', 'scheme_ref'])
fp_dist = fp_dist[(fp_dist['distance'] > 11) & (fp_dist['distance'] < 63)]
# fp_dist.pop('Unnamed: 0')
# fp_dist.pop('mirna_coorth')
# fp_dist.pop('seq')
# fp_dist.pop('specmir')
# fp_dist.pop('mirna_fam')
fp_dist = fp_dist.sort_values(by='distance', ascending=False)
print(fp_dist)

fp_dist.to_csv('/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/structure_dist_table.tsv', sep='\t', index=False)

