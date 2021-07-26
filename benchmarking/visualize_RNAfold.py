import json
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# json_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\all05_RNAfold_results.txt'
json_path = '/home/felixl/PycharmProjects/general/benchmarking/data/all05_RNAfold.json'

with open(json_path, 'r') as fh:
    sdict = json.load(fh)

df = pd.DataFrame.from_dict(sdict)
df['score'] = pd.to_numeric(df['score'])
print(df.columns)
for i in list(df.columns):
    print(df[i].head(1))

sns.histplot(data=df, x='score')
plt.title('Unfiltered score distribution')
plt.show()

plt.figure()