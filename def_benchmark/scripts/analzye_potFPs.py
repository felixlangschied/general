import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# initialize data
fp_path = '/home/felixl/PycharmProjects/general/def_benchmark/data/pot_FP_from_blast_folds.json'
with open(fp_path) as fh:
    fp_dict = json.load(fh)
dffp = pd.DataFrame.from_dict(fp_dict)
dffp = dffp.filter(['species', 'mirna', 'score', 'scheme'])
dffp['score'] = pd.to_numeric(dffp['score'])
print(dffp)

# # plot RNAfold score distribution
# sns.histplot(dffp, x='score')
# plt.show()
# filrow = dffp[(dffp['species'] == 'Dasypus_novemcinctus') & (dffp['mirna'] == 'Hsa-Mir-3126')]
# filrow = dffp[(dffp['species'] == 'Mus_musculus') & (dffp['mirna'] == 'Hsa-Mir-2355')]

# print(filrow['scheme'].to_string())
