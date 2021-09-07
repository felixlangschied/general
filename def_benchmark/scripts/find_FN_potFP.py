import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter

nc_path = '/home/felixl/PycharmProjects/general/def_benchmark/data/ncortho.json'
mirgene_path = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/mirgenedb_rnafold.json'
with open('/home/felixl/project/ncOrtho/benchmark/vertebrate_list.txt', 'r') as vh:
    vertebrates = [spec.strip().replace(' ', '_') for spec in vh]
result_dir = '/home/felixl/project/ncOrtho/benchmark/ncortho'
map_path = '/home/felixl/project/ncOrtho/benchmark/mirgenedb_specs.txt'

all_specs = []
with open(map_path, 'r') as mh:
    for line in mh:
        all_specs.append(line.strip().split('\t')[1])

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

mirgene = load_df(mirgene_path, 'miRGeneDB')
# only families that are present in human can be detected and count therefore as False negatives
fam_in_hsa = mirgene['mirna_fam'][mirgene['species'] == 'Homo_sapiens'].unique()
mirgene = mirgene[mirgene['mirna_fam'].isin(fam_in_hsa)]
mirgene = mirgene[mirgene['species'].isin(all_specs)]


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
# print(comp_df)

# extract FNs
FN = comp_df[comp_df['in_ncortho'].isnull()]
fn_list = [element.split('|') for element in FN.index]


# print(fn_list)
reasons = []
with open('/home/felixl/PycharmProjects/general/def_benchmark/data/FN_list.txt', 'w') as fnh:
    for species, family in fn_list:
        if family in fam_in_hsa:
            with open(f'{result_dir}/{species}/{species}.log', 'r') as resh:
                for line in resh:
                    if family in line:
                        reasons.append(line.strip().split('\t')[1])
reason_count = Counter(reasons)
print(reason_count)
            # fnh.write(f'{species}\t{family}\n')

FN_count = Counter([element.split('|')[0] for element in FN.index])
FN_df = pd.DataFrame.from_dict(FN_count, orient='index', columns=['count'])
FN_per_spec = FN_df.sort_values(by='count', ascending=False)
#plot
sns.set_style('whitegrid')
sns.barplot(data=FN_per_spec, y=FN_per_spec.index, x='count', orient='h')
plt.tight_layout()
plt.xlabel('miRNA families not found although present in MirGeneDB')


# print(FN_per_spec)

pot_FP = comp_df[comp_df['in_mirgenedb'].isnull()]
# print(pot_FP)
FP_dict = {'species': [], 'mirfam': []}
for ind in pot_FP.index:
    species, mirfam = ind.split('|')
    FP_dict['species'].append(species)
    FP_dict['mirfam'].append(mirfam)
# FP_count = Counter([element.split('|')[0] for element in pot_FP.index])
# print(FP_dict)
FP_df = pd.DataFrame.from_dict(FP_dict)
# FP_df = FP_df.sort_values(by='count', ascending=False)

FP_count = Counter(FP_df['species'])
FP_count_df = pd.DataFrame.from_dict(FP_count, orient='index', columns=['count'])
FP_count_df = FP_count_df.sort_values(by='count', ascending=False)
plt.figure()
sns.barplot(data=FP_count_df, y=FP_count_df.index, x='count', orient='h')
plt.tight_layout()
plt.xlabel('potential FPs')
plt.show()
print(FP_count_df)

fp_out = '/home/felixl/PycharmProjects/general/def_benchmark/data/fams_found_but_not_in_mirgenedb.txt'
FP_df.to_csv(fp_out, sep='\t', index=False)
