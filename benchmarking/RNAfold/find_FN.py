import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter

nc_path = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/ncortho.json'
mirgene_path = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/mirgenedb_rnafold.json'
with open('/home/felixl/project/ncOrtho/benchmark/vertebrate_list.txt', 'r') as vh:
    vertebrates = [spec.strip().replace(' ', '_') for spec in vh]


def load_df(path, name):
    with open(path, 'r') as fh:
        df_dict = json.load(fh)
        df = pd.DataFrame.from_dict(df_dict)
        df['score'] = pd.to_numeric(df['score'])
        df['length'] = df['seq'].apply(len)
        df['type'] = name
        df['species'] = df['species'].str.replace('>', '')
        df['species'] = df['species'].str.replace('Sarcophilus_harrissii', 'Sarcophilus_harrisii')
        df['species'] = df['species'].str.replace('Canis_familiaris', 'Canis_lupus_familiaris')
        df['species'] = df['species'].str.replace('Strongylocentrotus_purpuratus_', 'Strongylocentrotus_purpuratus')
        # filter for vertebrates only
        # df = df[df['species'].isin(vertebrates)]
        return df


ncortho = load_df(nc_path, 'ncOrtho')
ncortho.pop('mirna_coorth')
mirgene = load_df(mirgene_path, 'miRGeneDB')


tmp = pd.DataFrame()
nc_df = pd.DataFrame()
a = np.char.array(ncortho['species'].values)
b = np.char.array(ncortho['mirna_fam'].values)
tmp['res'] = (a + b'|' + b).astype(str)
nc_df['res'] = tmp['res'].unique()
nc_df['in'] = True
nc_df = nc_df.set_index('res')

tmp = pd.DataFrame()
mir_df = pd.DataFrame()
a = np.char.array(mirgene['species'].values)
b = np.char.array(mirgene['mirna_fam'].values)
tmp['res'] = (a + b'|' + b).astype(str)
mir_df['res'] = tmp['res'].unique()
mir_df['in'] = True
# mir_df.set_index('res')


comp_df = nc_df.join(mir_df.set_index('res'), how='outer', lsuffix='_ncortho', rsuffix='_mirgenedb')
print(comp_df)

FN = comp_df[comp_df['in_ncortho'].isnull()]
FN_count = Counter([element.split('|')[0] for element in FN.index])
FN_df = pd.DataFrame.from_dict(FN_count, orient='index', columns=['count'])
FN_df = FN_df.sort_values(by='count')
# print(FN_df)

pot_FP = comp_df[comp_df['in_mirgenedb'].isnull()]
print(pot_FP)
FP_dict = {'species': [], 'mirfam': []}
for ind in pot_FP.index:
    species, mirfam = ind.split('|')
    FP_dict['species'].append(species)
    FP_dict['mirfam'].append(mirfam)
# FP_count = Counter([element.split('|')[0] for element in pot_FP.index])
print(FP_dict)
FP_df = pd.DataFrame.from_dict(FP_dict)
# FP_df = FP_df.sort_values(by='count', ascending=False)
print(FP_df)

fp_out = '/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/fams_found_but_not_in_mirgenedb.txt'
FP_df.to_csv(fp_out, sep='\t', index=False)