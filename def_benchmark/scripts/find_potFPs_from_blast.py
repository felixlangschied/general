import pandas as pd
import json
import numpy as np
from collections import Counter
import tempfile
import os
import subprocess as sp
import glob

results_path = '/home/felixl/project/ncOrtho/benchmark/ncortho/analysis/PhyloProfile.long'
result_dir = '/home/felixl/project/ncOrtho/benchmark/ncortho'
map_path = '/home/felixl/project/ncOrtho/benchmark/mirgenedb_specs.txt'
blast_path = '/home/felixl/PycharmProjects/general/def_benchmark/scripts/blast_results_in_mirgenedb.json'
outf = '/home/felixl/PycharmProjects/general/def_benchmark/data/pot_FP_from_blast_folds.json'
phylo_out = '/home/felixl/PycharmProjects/general/def_benchmark/data/FP_phyloprofile.long'


def load_df(path):
    with open(path, 'r') as fh:
        df_dict = json.load(fh)
        df = pd.DataFrame.from_dict(df_dict)
        df['species'] = df['species'].str.replace('>', '')
        df['species'] = df['species'].str.replace('Sarcophilus_harrissii', 'Sarcophilus_harrisii')
        df['species'] = df['species'].str.replace('Canis_familiaris', 'Canis_lupus_familiaris')
        df['species'] = df['species'].str.replace('Strongylocentrotus_purpuratus_', 'Strongylocentrotus_purpuratus')
        # filter for vertebrates only
        # df = df[df['species'].isin(vertebrates)]
        return df


name_2_taxid = {}
taxid_2_name = {}
with open(map_path, 'r') as mh:
    for line in mh:
        taxid, name = line.strip().split()
        name_2_taxid[name] = taxid
        taxid_2_name[taxid] = name

uniq_list = []
results_dict = {'species': [], 'miRNA': [], 'combined': []}
with open(results_path, 'r') as rh:
    next(rh)
    for line in rh:
        mirid, taxon, candidate = line.strip().split()
        taxname = taxid_2_name[taxon.replace('ncbi', '')]
        combined = f'{taxname}|{mirid}'
        if combined not in uniq_list:
            uniq_list.append(combined)
            results_dict['species'].append(taxname)
            results_dict['miRNA'].append(mirid)
            results_dict['combined'].append(combined)

results_df = pd.DataFrame.from_dict(results_dict)
# print(results_df['species'].unique())
results_df = results_df[['combined']]
# print(results_df)


blast_df = load_df(blast_path)
blast_df = blast_df[(blast_df['coverage'] > 80) & (blast_df['identity'] > 90)]
a = np.char.array(blast_df['species'].values)
b = np.char.array(blast_df['miRNA'].values)
blast_df['combined'] = (a + b'|' + b).astype(str)
blast_df = blast_df[['combined']]
blast_df['blast'] = True


# join tables
comp_df = results_df.join(blast_df.set_index('combined'), on='combined', how='outer')
pot_FP = comp_df[comp_df.isnull().any(axis=1)]
pot_FP = pot_FP['combined'].str.split('|', expand=True)
pot_FP.columns = ['species', 'mirid']
int_mirids = pot_FP['mirid'].unique()
with open('/home/felixl/PycharmProjects/general/def_benchmark/data/fp_mirnas.txt', 'w') as mirh:
    for element in int_mirids:
        mirh.write(f'{element}\n')
# pot_FP['species'] = pot_FP['combined'].str.split('|')[0]

# df with all ncOrtho results that were not confirmed with blast
# print(pot_FP)

FP_dict = {}
spec_count = Counter(pot_FP['species'])
print(spec_count)
for row in pot_FP.iterrows():
    species = list(row)[1].species
    mirid = list(row)[1].mirid
    if not species in FP_dict:
        FP_dict[species] = []
    FP_dict[species].append(mirid)


# df_dict = {'species': [], 'mirna': [], 'mirna_fam': [], 'score': [], 'seq': [], 'scheme': []}
# for species in FP_dict:
#     with open(f'{result_dir}/{species}/{species}_orthologs.fa') as resh:
#         for line in resh:
#             if line.startswith('>'):
#                 header = '|'.join(line.strip().split('|')[:2])
#                 mirid = line.split('|')[1].split('_')[0]
#                 mirfam = '-'.join(mirid.strip().split('-')[1:3]).replace('_pre', '')
#                 if mirid in FP_dict[species]:
#                     seq = next(resh)
#                     with tempfile.NamedTemporaryFile(mode='w+') as fp:
#                         fp.write(f'{header}\n{seq}\n')
#                         fp.seek(0)
#                         fold_cmd = f'RNAfold -i {fp.name}'
#                         res = sp.run(fold_cmd, shell=True, capture_output=True)
#
#                         out = res.stdout.decode('utf-8').strip()
#                         folddata = out.split('\n')
#                         schemscore = folddata[-1]
#                         scheme, score = schemscore.split(' (')
#                         score = score.replace(')', '')
#                         df_dict['species'].append(species)
#                         df_dict['mirna'].append(mirid)
#                         df_dict['mirna_fam'].append(mirfam)
#                         df_dict['score'].append(score)
#                         df_dict['seq'].append(seq)
#                         df_dict['scheme'].append(scheme)
#                         header = header.replace('|', '_')
# psfiles = glob.glob(f'{os.getcwd()}/*.ps')
# for file in psfiles:
#     os.remove(file)
# with open(outf, 'w') as of:
#     json.dump(df_dict, of)


results_df = results_df['combined'].str.split('|', expand=True)
results_df.columns = ['species', 'mirid']
# print(results_df)
print(name_2_taxid)

with open(phylo_out, 'w') as pph:
    pph.write(f'geneID\tncbiID\torthoID\tis_pot_FP\n')
    for row in results_df.iterrows():
        species = list(row)[1].species
        mirid = list(row)[1].mirid

        if species in FP_dict and mirid in FP_dict[species]:
            fp = '0'
        else:
            fp = '1'
        taxid = f'ncbi{name_2_taxid[species]}'
        pph.write(f'{mirid}\t{taxid}\t{mirid}\t{fp}\n')
