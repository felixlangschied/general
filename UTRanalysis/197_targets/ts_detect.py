import json
import subprocess as sp
import os
import pandas as pd

human_path = '/home/felixl/PycharmProjects/general/UTRanalysis/197_targets/target_utrs/GCF_000001405.39_GRCh38.p13_UTRs.json'
mouse_path = '/home/felixl/PycharmProjects/general/UTRanalysis/197_targets/target_utrs/GCF_000001635.27_GRCm39_UTRs.json'

outdir = '/home/felixl/PycharmProjects/general/UTRanalysis/197_targets/data'
# target_path = '/home/felixl/PycharmProjects/general/UTRanalysis/197_targets/potential_targets.txt'
target_path = '/home/felixl/PycharmProjects/general/UTRanalysis/197_targets/data/matthias_targets.txt'

pot_ts = []
with open(target_path, 'r') as th:
    for line in th:
        pot_ts.append(line.strip())

if not os.path.isdir(outdir):
    os.mkdir(outdir)


def string_analysis(utr_dict, strings):
    # df_dict = {'gene': [], 'score': [], 'formula': []}
    df_dict = {'gene': [], 'score': []}
    for protname in utr_dict:
        # print(f'### human: {protname} ###')
        counter = 0
        for iso_count, seq in enumerate(utr_dict[protname].values(), 1):
            # print('# {}'.format(iso_id))
            for motif in strings:
                if motif in seq:
                   counter += 1
        if iso_count == 0:
            score = 'NaN'
            formula = 'NaN'
            continue
        else:
            formula = '{} / {}'.format(counter, iso_count)
            score = counter / iso_count
        df_dict['gene'].append(protname)
        df_dict['score'].append(score)
        # df_dict['formula'].append(formula)
    df = pd.DataFrame.from_dict(df_dict)
    return df


list_197 = ['TGGTGAA', 'TGGTGAG']
list_769 = ['AGGTCTC', 'AGGTTTC']
# str_list = ['TGGTGA']
with open(human_path, 'r') as hh:
    raw_d = json.load(hh)
    human_d = {k: v for k, v in raw_d.items() if k in pot_ts}
with open(mouse_path, 'r') as mh:
    raw_d = json.load(mh)
    mouse_d = {k.upper(): v for k, v in raw_d.items() if k.upper() in pot_ts}

hu_df = string_analysis(human_d, list_197)
mu_df = string_analysis(mouse_d, list_197)


df_197 = hu_df.join(mu_df.set_index('gene'), on="gene", lsuffix='_hsa', rsuffix='_mmu', how='inner')
hu_df = string_analysis(human_d, list_769)
mu_df = string_analysis(mouse_d, list_769)
df_769 = hu_df.join(mu_df.set_index('gene'), on="gene", lsuffix='_hsa', rsuffix='_mmu', how='inner')

print(df_197.head())
print(df_769.head())
df = df_197.join(df_769.set_index('gene'), on="gene", lsuffix='_197', rsuffix='_769', how='outer')
print(df)
# df = df.dropna(how='any')
# print(df.size)
# df = df.sort_values(by=['score_mmu', 'score_hsa'], ascending=False)
# # print(df.head(30))
# df.to_csv(f'{outdir}/result.tsv', sep='\t')
df.to_csv('/home/felixl/PycharmProjects/general/UTRanalysis/matthias_targets_scores.txt', sep='\t', index=False)


# above_1 = df['score_hsa'] >= 1
# fil_df = df[above_1].sort_values(by=['score_mmu', 'score_hsa'], ascending=False)
# print(fil_df)
# fil_df.to_csv(f'{outdir}/filtered_string_result.tsv', sep='\t')
# # print(f'{outdir}/string_result.tsv')






