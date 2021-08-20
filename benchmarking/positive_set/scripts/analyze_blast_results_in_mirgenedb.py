import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

inpath = r'C:\Users\felix\PycharmProjects\general\benchmarking\positive_set\data\blast_results_in_mirgenedb.json'

with open(inpath, 'r') as fh:
    df_dict = json.load(fh)

df = pd.DataFrame.from_dict(df_dict)
df['coverage'] = df['coverage'].apply(round)
df['identity'] = df['identity'].apply(round)
a = np.char.array(df['identity'].values)
b = np.char.array(df['coverage'].values)

df['idcov'] = (a + b'_' + b).astype(str)
df['count'] = df.groupby('idcov')['idcov'].transform('count')
# df['idcov'] = df.identity.map(str) + '_' + df.coverage
print(df)

sns.set_style('darkgrid')
sns.scatterplot(data=df, x='identity', y='coverage', size='count', sizes=(30, 300))
plt.tight_layout()
plt.show()
