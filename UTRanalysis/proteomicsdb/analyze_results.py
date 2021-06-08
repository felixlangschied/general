import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

res_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\proteomicsdb\proteomics_results.csv'
# res_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\proteomicsdb\cell_line_only_proteomics.csv'

# target_name_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\list_of_connected_targetprots.txt'
target_name_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\hsa_rodentmis_genenames.txt'

target_protnames = []
with open(target_name_path, 'r') as fh:
    for line in fh:
        target_protnames.append(line.strip())
print(target_protnames)

# bot_exp_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\proteomicsdb\exp_mouse_human.txt'
bot_exp_path = r'C:\Users\felix\PycharmProjects\general\UTRanalysis\curated_connected_targetprots.txt'
to_keep = []
with open(bot_exp_path, 'r') as fh:
    for line in fh:
        to_keep.append(line.strip())

# LOAD DATA FROM DISK
df = pd.read_csv(res_path)

# FILTER FOR TISSUES THAT HAVE DATA IN HUMAN AND IN MOUSE
# print(df.head())
# tiss_in_mouse = df['tissue'] if df['organism'] == 'mouse'
# tiss_in_mouse = df[df.organism == 'mouse'].tissue.unique()
# # print(tiss_in_mouse)
# keep_tiss = df['tissue'].isin(tiss_in_mouse)
# tiss_df = df[keep_tiss]
# sns.catplot(x='tissue', y='norm_intensity', hue='organism', data=tiss_df, kind='violin', legend=False)

df = df.sort_values(by='organism')
sns.catplot(x='tissue', y='norm_intensity', hue='organism', data=df, kind='box', legend=False)
plt.xticks(rotation=60, ha='right')
# plt.xticks(rotation=30, ha='right')
# ax = plt.gca()
# locs, labels = plt.xticks()
# ax.set_xticks(np.add(0.5, locs), minor=True)
# plt.grid(axis='x', which='minor')
plt.ylabel('log 10 normalized iBAQ intensity')
plt.legend(loc='lower right')
plt.tight_layout()

# # PLOT EXPRESSION PER TISSUE
# sns.catplot(x='tissue', y='norm_intensity', hue='organism', data=df, kind='box')
# plt.xticks(rotation=45, ha='right')
# plt.tight_layout()





# initialize manual filter of proteins
# to_remove = ['IGF1R', 'BMPR1A', 'RORA', 'TWIST1', 'HMGA2']
# to_remove = ['PIM1', ]
# to_keep = [entry for entry in target_protnames if entry not in to_remove]
# do filter
keeper = df['gene_name'].isin(to_keep)
df = df[keeper]



# plot expression per protein

plt.rcParams['font.size'] = '16'
sns.catplot(x='gene_name', y='norm_intensity', hue='organism', data=df, kind='violin', legend=False)
plt.ylabel('log 10 normalized iBAQ intensity')
plt.xticks(rotation=30, ha='right')
ax = plt.gca()
locs, labels = plt.xticks()
ax.set_xticks(np.add(0.5, locs), minor=True)
plt.grid(axis='x', which='minor')


plt.xlabel('')
plt.legend(loc='lower center')
plt.tight_layout()
plt.show()
