import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

cm_dir = '/share/project/felixl/ncOrtho/data/mirGeneDB/coreset/output/CMs'


miss_path = '/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/fam_missing_TP.tsv'
missed = []
with open(miss_path, 'r') as fh:
    next(fh)
    for line in fh:
        data = line.strip().split()
        if data[-1] == '0':
            continue
        missed.append(data[0])

# files = glob.glob(f'{cm_dir}/*.log')
# df_dict = {'mirna': [], 'fam': [], 'count': []}
# for file in files:
#     with open(file, 'r') as fh:
#         for line in fh:
#             if not line.startswith('#'):
#                 # print(line.split())
#                 data = line.strip().split()
#                 fam = '-'.join(data[1].split('-')[1:3]).replace('_pre', '')
#                 df_dict['mirna'].append(data[1])
#                 df_dict['fam'].append(fam)
#                 df_dict['count'].append(int(data[2]) - 1)
# count_df = pd.DataFrame.from_dict(df_dict)
# count_df.to_csv('/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/training_count.csv')
count_df = pd.read_csv('/home/felixl/PycharmProjects/general/benchmarking/positive_set/data/training_count.csv')
print(count_df.head())
missed_df = count_df[count_df['fam'].isin(missed)]
print(missed_df.head())

sns.set_style('dark')
# sns.histplot(count_df, x='count', discrete=True)
# plt.xlabel('Number of core-orthologs used for training')
#
# plt.figure()
# sns.histplot(missed_df, x='count', discrete=True)
# plt.xlabel('Number of core-orthologs used for training')
# plt.title('miRNA families with missed TPs')


plt.figure()
sns.histplot(data=count_df, x='count', discrete=True)
sns.histplot(data=missed_df, x='count', discrete=True, color='orange')
plt.xlabel('Number of core-orthologs used for training')
plt.legend(['All miRNAs', 'miRNA families with missed TPs'])

# plt.figure()
# sns.histplot(data=count_df, x='count', discrete=True, stat='density')
# sns.histplot(data=missed_df, x='count', discrete=True, stat='density', color='orange')
# plt.xlabel('Number of core-orthologs used for training')
# plt.legend(['All miRNAs', 'miRNA families with missed TPs'])
plt.show()
