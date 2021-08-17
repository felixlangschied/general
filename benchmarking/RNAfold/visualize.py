import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

nc_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\ncortho.json'
mirbh_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\mirbh_folds.json'
filtered_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\filtered_mirbh_folds.json'
mirgene_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\mirgenedb_rnafold.json'


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
        return df


def fam_count(df, spec):
    spec_df = df[df['species'] == spec]
    fam_count = len(spec_df['mirna_fam'].unique())
    return fam_count

ncortho = load_df(nc_path, 'ncOrtho')
ncortho.pop('mirna_coorth')
print(ncortho['length'])
mirbh = load_df(mirbh_path, 'mirbh')
mirbh = mirbh[mirbh['length'] > 50]
mirgene = load_df(mirgene_path, 'miRGeneDB')
mirgene_famset = mirgene[mirgene['mirna_fam'].isin(mirbh['mirna_fam'].unique())]

# concatenate dataframes
ges = pd.concat([mirgene, ncortho, mirbh])
# filter for families for which orthologs  have been detected by mirbh
fam_set = ges[ges['mirna_fam'].isin(mirbh['mirna_fam'].unique())]

# count families
c_dict = {'species': [], 'mirgenedb': [], 'mirbh': [], 'ncortho': []}
for species in mirbh['species'].unique():
    c_dict['species'].append(species)
    c_dict['mirgenedb'].append(fam_count(mirgene_famset, species))
    c_dict['mirbh'].append(fam_count(mirbh, species))
    c_dict['ncortho'].append(fam_count(ncortho, species))


count_df = pd.DataFrame.from_dict(c_dict).sort_values(by=['mirgenedb', 'mirbh', 'ncortho'], ascending=False)
print(count_df)



# plot length distribution of orthologs
plt.subplot(211)
# sns.set_style('dark')
# sns.histplot(data=fam_set, x='length', hue='type', element='poly', fill=False)
sns.kdeplot(data=fam_set, x='length', hue='type')
plt.xlabel('Length [nt]')
plt.xlim([40, 80])
# plt.figure()

# plot score distribution
plt.subplot(212)
# sns.histplot(data=fam_set, x='score', hue='type', element='poly', fill=False)
sns.kdeplot(data=fam_set, x='score', hue='type')
plt.xlabel('RNAfold minimum free energy')
plt.xlim([-40, -10])
plt.show()

# plot ortholog families per species
# is_unique = pd.unique(mirbh['mirna_fam'])
# print(len(is_unique))