import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats

# https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.kstest.html

def shorten(seq):
    if len(seq) > 80:
        seq = seq[:80]
        seq = seq + '@'
    return seq


# json_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\all05_RNAfold_results.txt'
json_path = '/benchmarking/data/all05_RNAfold.json'
six_path = '/benchmarking/data/all06_RNAfold.json'
avail_path = '/benchmarking/data/present.txt'
ref_path = '/benchmarking/data/reference_RNAfold.json'
out_dir = '/benchmarking/data/rnafold_analysis'
mouse_path = '/benchmarking/data/mmu_reference_RNAfold.json'
mirbase_path = '/benchmarking/data/hsa_mirbase_RNAfold.json'

with open(json_path, 'r') as fh:
    sdict = json.load(fh)
with open(avail_path, 'r') as fh:
    available = fh.read().split('\n')
    available.pop()
with open(ref_path, 'r') as fh:
    ref_dict = json.load(fh)
with open(mouse_path, 'r') as fh:
    mmu_dict = json.load(fh)
with open(mirbase_path, 'r') as fh:
    mirbase_dict = json.load(fh)
with open(six_path, 'r') as fh:
    six_dict = json.load(fh)

# initialize dataframes
df = pd.DataFrame.from_dict(sdict)
ref_df = pd.DataFrame.from_dict(ref_dict)
ref_df['score'] = pd.to_numeric(ref_df['score'])
mmu_df = pd.DataFrame.from_dict(mmu_dict)
mmu_df['score'] = pd.to_numeric(mmu_df['score'])
mirbase_df = pd.DataFrame.from_dict(mirbase_dict)
mirbase_df['score'] = pd.to_numeric(mirbase_df['score'])
six_df = pd.DataFrame.from_dict(six_dict)
six_df['score'] = pd.to_numeric(six_df['score'])
print(df.columns)

# analyze score column
df['score'] = pd.to_numeric(df['score'])
scores = df['score'].to_numpy()
print(f'mean: {scores.mean()}')
print(f'median: {np.median(scores)}')
std = scores.std()
mean = scores.mean()
threesig = mean - 3 * std
plustwosig = mean + 2 * std
plusthreesig = mean + 3 * std
print(f'standart deviation: {std}')
print(f'2sig: {mean - 2 * std}')
print(f'plus2sig: {plustwosig}')
print(f'3sig: {mean - 3 * std}')
print(f'plus2sig: {plusthreesig}')

six_scores = six_df['score'].to_numpy()
six_thres = six_scores.mean() - 3 * six_scores.std()

# filter for human
human = df[df['species'] == 'Homo_sapiens']
# filter for 19 specs
df = df[df['species'].isin(available)]
# filter for score
filtered_five = df[df['score'] >= threesig]
filtered_six = six_df[six_df['score'] >= six_thres]

# # plotting
# sns.histplot(data=df, x='score')
# plt.xlabel('RNAfold minimum free energy')
#
# plt.figure()
# sns.histplot(data=filtered, x='score')
# plt.xlabel('RNAfold minimum free energy')
# plt.title('Outliers filtered')
#
# plt.figure()
# sns.histplot(data=human, x='score')
# plt.xlabel('RNAfold minimum free energy')
# plt.title('Human ncOrtho results')
#
# plt.figure()
# sns.histplot(data=ref_df, x='score')
# plt.xlabel('RNAfold minimum free energy')
# plt.title('Reference miRNAs')
# plt.show()

# distribution plot
sns.kdeplot(data=filtered_five, x='score')
sns.kdeplot(data=human, x='score')
sns.kdeplot(data=ref_df, x='score')
sns.kdeplot(data=mirbase_df, x='score')
sns.kdeplot(data=mmu_df, x='score')
plt.legend(['All results (mirGeneDB)', 'Human ncOrtho results (mirGeneDB)', 'Human mirGeneDB miRNAs', 'Human miRBase miRNAs', 'Mouse miRBase miRNAs'])
plt.xlabel('RNAfold minimum free energy')
plt.show()

# sns.kdeplot(data=filtered_five, x='score')
# sns.kdeplot(data=filtered_six, x='score')
# plt.legend(['0.5 CMsearch cutoff', '0.6 CMsearch cutoff'])
# plt.show()


small_score = df[df['score'] <= threesig]
print(small_score)
fams = set()
for element in small_score['mirna_fam']:
    fams.add(element)
# for scheme in small_score['scheme']:
#     print(scheme)
# for seq in small_score['seq']:
#     print(seq)
print(fams)


repfams = df[df['mirna_fam'].isin(list(fams))].sort_values(by='score')
# filter columns
repfams = repfams.filter(items=['species', 'mirna_coorth', 'score'])
# print(repfams)

exfam = df[df['mirna'] == 'Hsa-Mir-154-P4b'].sort_values(by='score')
exfam['length'] = exfam['seq'].apply(len)
exfam = exfam.filter(items=['species', 'mirna_coorth', 'score', 'length', 'seq', 'scheme'])
exfam['seq'] = exfam['seq'].apply(shorten)
exfam['scheme'] = exfam['scheme'].apply(shorten)
exfam.to_csv(f'{out_dir}/example_table.tsv', sep='\t', index=False)
highlights = pd.concat([exfam.iloc[-3:-1, :], exfam.iloc[1:3, :]])
highlights.to_csv(f'{out_dir}/highlights_table.tsv', sep='\t', index=False)

# print(exfam.filter(items=['species', 'mirna_coorth', 'score', 'length']))
# for scheme in exfam['scheme']:
#     print(scheme)
# for seq in exfam['seq']:
#     print(seq)

# highlights = pd.concat([exfam.iloc[1:3, :], exfam.iloc[-3:-1, :]])
# highlights['seq'] = highlights['seq'].apply(shorten)




# below MIRNAminer threshold
# below = df[df['score'] >= -25]
# print(below)