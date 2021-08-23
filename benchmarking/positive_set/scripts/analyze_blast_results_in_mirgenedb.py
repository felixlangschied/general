import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# inpath = r'C:\Users\felix\PycharmProjects\general\benchmarking\positive_set\data\blast_results_in_mirgenedb.json'
inpath = '/home/felixl/project/ncOrtho/benchmark/filtered_mirbh/analysis/blast_results_in_mirgenedb.json'
outf = '/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/ncortho_res_below_90i_80c.txt'
pot_FP = '/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/fams_found_but_not_in_mirgenedb.txt'

with open(inpath, 'r') as fh:
    df_dict = json.load(fh)

# load potential FPs, miRNA families that were detected in species for which no mirgenedb entry exists
FP = pd.read_csv(pot_FP, sep='\t')
a = np.char.array(FP['species'].values)
b = np.char.array(FP['mirfam'].values)
FP['specfam'] = (a + b'|' + b).astype(str)
print(FP)

# load BLAST results
df = pd.DataFrame.from_dict(df_dict)
mirnas = list(df['miRNA'])
mirfam = ['-'.join(mirna.split('-')[1:3]) for mirna in mirnas]
df['mirfam'] = mirfam
a = np.char.array(df['species'].values)
b = np.char.array(df['mirfam'].values)
df['specfam'] = (a + b'|' + b).astype(str)
# print(df)

# find potential FPs in blast results
# bl_fp = FP.join(df.set_index('specfam'), on='specfam', how='outer', lsuffix='_fp', rsuffix='_blast')
bl_fp = df[df['specfam'] == 'Sarcophilus_harrisii|Mir-1287']
# print(bl_fp)
#



# apply threshold
below_t = df[(df['identity'] < 90) | (df['coverage'] < 80)]
print(below_t['mirfam'].unique())

# below_t.to_csv(outf, sep='\t', index=False)






# # plot
# a = np.char.array(df['identity'].values)
# b = np.char.array(df['coverage'].values)
# df['idcov'] = (a + b'_' + b).astype(str)
# df['count'] = df.groupby('idcov')['idcov'].transform('count')
# sns.set_style('darkgrid')
# sns.scatterplot(data=df, x='identity', y='coverage', size='count', sizes=(30, 300))
# plt.tight_layout()
# plt.show()

