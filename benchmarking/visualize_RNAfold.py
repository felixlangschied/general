import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# json_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\all05_RNAfold_results.txt'
json_path = '/home/felixl/PycharmProjects/general/benchmarking/data/all05_RNAfold.json'

with open(json_path, 'r') as fh:
    sdict = json.load(fh)

plt.clf()
df = pd.DataFrame.from_dict(sdict)
df['score'] = pd.to_numeric(df['score'])
print(df.columns)
# for i in list(df.columns):
#     print(df[i].head(1))

# sns.histplot(data=df, x='score')
# plt.title('Unfiltered score distribution')
# plt.show()

# analyze score array
scores = df['score'].to_numpy()
print(f'mean: {scores.mean()}')
print(f'median: {np.median(scores)}')
std = scores.std()
mean = scores.mean()
threesig = mean - 3 * std

print(f'standart deviation: {std}')
print(f'2sig: {mean - 2 * std}')
print(f'3sig: {mean - 3 * std}')


plt.figure()
small_score = df[df['score'] <= threesig]
fams = set()
for element in small_score['mirna_fam']:
    fams.add(element)
print(fams)
# print(small_score['scheme'])
# print(small_score)


# below MIRNAminer threshold
below = df[df['score'] >= -25]
print(below)