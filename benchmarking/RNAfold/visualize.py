import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# nc_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\ncortho.json'
# mirbh_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\mirbh_folds.json'
# filtered_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\filtered_mirbh_folds.json'
# mirgene_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\data\mirgenedb_rnafold.json'
nc_path = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/ncortho.json'
mirbh_path = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/mirbh_folds.json'
filtered_path = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/filtered_mirbh_folds.json'
mirgene_path = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/mirgenedb_rnafold.json'
filtered = False
style = 'whitegrid'


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


def make_distplot(df, sty, fil):
    sns.set_style(sty)
    plt.subplot(211)
    # sns.histplot(data=fam_set, x='length', hue='type', element='poly', fill=False)
    sns.kdeplot(data=df, x='length', hue='type', hue_order=['miRGeneDB', 'ncOrtho', 'mirbh_length_cutoff', 'mirbh'])
    plt.tight_layout()
    plt.xlabel('Length [nt]')

    if fil:
        plt.title('Outliers filtered')
        plt.xlim([40, 80])
    else:
        plt.xlim([20, 100])
        # plt.xlim([140, 1000])
        # plt.ylim([-0.001, 0.06])

    # plot score distribution
    plt.subplot(212)
    # sns.histplot(data=fam_set, x='score', hue='type', element='poly', fill=False)
    sns.kdeplot(data=df, x='score', hue='type', hue_order=['miRGeneDB', 'ncOrtho', 'mirbh_length_cutoff', 'mirbh'])

    plt.xlabel('RNAfold minimum free energy')

    if fil:
        # plt.title('Length outside 3sig interval filtered out')
        plt.xlim([-40, -10])
    else:
        plt.xlim([-60, 0])
        # plt.ylim([-0.001, 0.05])
    plt.tight_layout()
    plt.show()


def fam_count(df, spec):
    spec_df = df[df['species'] == spec]
    fam_count = len(spec_df['mirna_fam'].unique())
    return fam_count


ncortho = load_df(nc_path, 'ncOrtho')
ncortho.pop('mirna_coorth')
mirbh = load_df(mirbh_path, 'mirbh')
filtered = load_df(filtered_path, 'mirbh_length_cutoff')
print(mirbh['mirna'].head())
print(ncortho['mirna'].head())
mirgene = load_df(mirgene_path, 'miRGeneDB')
mirgene_famset = mirgene[mirgene['mirna_fam'].isin(mirbh['mirna_fam'].unique())]
# show largest miRNAs in database
# print(mirgene_famset.sort_values(by='length', ascending=False).head(10))

# concatenate dataframes
ges = pd.concat([mirgene, ncortho, mirbh, filtered])
# filter for families for which orthologs  have been detected by mirbh
fam_set = ges[ges['mirna_fam'].isin(mirbh['mirna_fam'].unique())]

# plot distribution plot
make_distplot(fam_set, style, fil=False)

# calculate length cutoff
length = mirgene_famset['length'].to_numpy()
m_threesig = length.mean() - 3 * length.std()
p_threesig = length.mean() + 3 * length.std()
print(f'3sig: {m_threesig} - {p_threesig}')


# determine number of sequences above 3sig level for each dataset
large_mirgene = mirgene['length'][mirgene['length'] > p_threesig].size
large_mirbh = mirbh['length'][mirbh['length'] > p_threesig].size
large_ncortho = ncortho['length'][ncortho['length'] > p_threesig].size
print(large_mirgene)
print(large_ncortho)
print(large_mirbh)
# determine number of sequences below 3sig level for each dataset
small_mirgene = mirgene['length'][mirgene['length'] < m_threesig].size
small_mirbh = mirbh['length'][mirbh['length'] < m_threesig].size
small_ncortho = ncortho['length'][ncortho['length'] < m_threesig].size
print(small_mirgene)
print(small_ncortho)
print(small_mirbh)


# filter dataframes fro length cutoff
# mirbh = mirbh[(mirbh['length'] > m_threesig) & (mirbh['length'] < p_threesig)]
mirbh = mirbh[(mirbh['length'] > m_threesig)]
# ncortho = ncortho[(ncortho['length'] > m_threesig) & (ncortho['length'] < p_threesig)]

# concatenate
ges = pd.concat([mirgene, ncortho, mirbh])
# filter for families for which orthologs  have been detected by mirbh
fil_fam_set = ges[ges['mirna_fam'].isin(mirbh['mirna_fam'].unique())]

# plot length distribution of orthologs
# make_distplot(fil_fam_set, style, fil=True)

# plot ortholog families per species
# count families
# c_dict = {'species': [], 'mirgenedb': [], 'mirbh': [], 'ncortho': []}
# for species in mirbh['species'].unique():
#     c_dict['species'].append(species)
#     c_dict['mirgenedb'].append(fam_count(mirgene_famset, species))
#     c_dict['mirbh'].append(fam_count(mirbh, species))
#     c_dict['ncortho'].append(fam_count(ncortho, species))
# count_df = pd.DataFrame.from_dict(c_dict).sort_values(by=['mirgenedb', 'mirbh', 'ncortho'], ascending=False)
# print(count_df)

c_vis = {'species': [], 'count': [], 'source': []}
for species in mirbh['species'].unique():
    c_vis['species'].append(species)
    c_vis['count'].append(fam_count(mirgene_famset, species))
    c_vis['source'].append('mirgenedb')

    c_vis['species'].append(species)
    c_vis['count'].append(fam_count(mirbh, species))
    c_vis['source'].append('mirbh_length_cutoff')

    c_vis['species'].append(species)
    c_vis['count'].append(fam_count(ncortho, species))
    c_vis['source'].append('ncortho')
vis_df = pd.DataFrame.from_dict(c_vis)

plt.figure()
sns.set_style(style)
count_sorted = vis_df.copy().sort_values(by=['count'], ascending=False)
sns.barplot(data=count_sorted, y='species', x='count', hue='source', orient='h')
plt.tight_layout()

plt.figure()
sns.set_style(style)
ever_second = [True, True, True, False, False, False] * int((len(vis_df) / 3) / 2)
small_count = vis_df.copy()[ever_second].sort_values(by=['count'], ascending=False)
sns.barplot(data=small_count, y='species', x='count', hue='source', orient='h', hue_order=['mirgenedb', 'ncortho', 'mirbh_length_cutoff'])
# plt.xticks(rotation=30, ha='right')
plt.tight_layout()

plt.show()
