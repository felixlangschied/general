import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering


# outpath = 'UTRanalysis'
taxid_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\length_distribution\taxid_name_assembly.tsv'
taxord_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\length_distribution\prim_rod_act.phyloprofile'
phylo_output = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\ingos_run.phyloprofile'
fa_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\tmp_data\{}.fa'
prot_map = r'uniprot_genename.txt'
con_target_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\list_of_connected_targetprots.txt'


con_targets = []
with open(con_target_path, 'r') as fh:
    for line in fh:
        con_targets.append(line.strip())


tax_2_name = {}
name_2_tax = {}
with open(taxid_path, 'r') as fh:
    for line in fh:
        taxid, spec_name, assembly = line.strip().split('\t')
        if len(spec_name.split('_')) > 2:
            spec_name = '_'.join(spec_name.split('_')[0:2])
        tax_2_name[taxid] = spec_name
        name_2_tax[spec_name] = taxid

tax_order = []
with open(taxord_path, 'r') as fh:
    next(fh)
    for line in fh:
        taxid = line.strip().split('\t')[1].replace('ncbi', '')
        tax_name = tax_2_name[taxid]
        if not tax_name in tax_order:
            tax_order.append(tax_name)

ens_2_name = {}
name_2_ens = {}
with open(prot_map, 'r') as fh:
    for line in fh:
        ensembl, name = line.strip().split()
        ens_2_name[ensembl] = name
        name_2_ens[name] = ensembl

prot_names = set()
tax_ids = set()
with open(phylo_output, 'r') as fh:
    next(fh)
    for line in fh:
        prot_names.add(line.strip().split('\t')[0])
        tax_ids.add(line.strip().split('\t')[1])

prot_names = list(prot_names)
prot_names = [ens_2_name[ens_name] for ens_name in prot_names]
df = pd.DataFrame(index=tax_order, columns=prot_names)


for prot in prot_names:
    ens_id = name_2_ens[prot]
    if os.path.isfile(fa_path.format(ens_id)):
        with open(fa_path.format(ens_id), 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    ncbi_id = line.strip().replace('>ncbi', '')
                    row_id = tax_2_name[ncbi_id]
                else:
                    length = len(line.strip())
                    df.loc[row_id, prot] = length
    else:
        print('# Could not find file {}'.format(fa_path.format(prot)))

print(df.head())

# filter for only targets that are connected by at least two interactions
df = df.loc[:, con_targets]



plt.clf()
# df = df.dropna(how='all', axis=0)
# df = df.dropna(how='all', axis=1)



# plt.show()

# # create length distribution phyloprofile output
# t_df = df.transpose()
# max_df = t_df.max()
# t_df = t_df.divide(max_df)
# t_df = t_df.applymap(lambda x: 'genbla#' + str(x) + '#1')
# print(t_df)
# t_df.to_csv('utr_length_phyloprofile.matrix', sep='\t', index_label='geneID')

# create heatmap

t_df = df.transpose()
mask = t_df.isna()
t_df = t_df.fillna(0)
# cluster
num_df = t_df.to_numpy()
clustering = AgglomerativeClustering(n_clusters=3).fit(num_df)
# print(clustering.labels_)
t_df['cluster'] = clustering.labels_
t_df = t_df.sort_values(by='cluster')
sns.heatmap(data=t_df, xticklabels=True, yticklabels=True, mask=mask)
# plt.xticks(rotation=20, ha='right')
plt.tight_layout()
plt.title("Length of 3'-UTR [nt]")


# create boxplot
plt.figure()
sns.boxplot(data=t_df.transpose())
plt.xticks(rotation=45, ha='right')
plt.ylabel("Length of 3'-UTR [nt]")
plt.tight_layout()
plt.show()

#
# plt.figure()
# sns.heatmap(data=df.transpose(), xticklabels=True, yticklabels=True)
# plt.show()
# for fasta in fa_files:


