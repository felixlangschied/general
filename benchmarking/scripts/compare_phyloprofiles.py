import pandas as pd
import json

# ortho_path = '/benchmarking/phylo_input/family_results_19specs.long'
# mirgene_path = '/benchmarking/phylo_input/family_mirgenedb_19specs.long'
# map_path = r'/benchmarking/data/available_specs.txt'
# # ortho_path = '/home/felixl/PycharmProjects/general/benchmarking/phylo_output/family_results_19specs.txt'
# out_dir = '/benchmarking/data'

ortho_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\phylo_input\family_results_19specs.long'
mirgene_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\phylo_input\family_mirgenedb_19specs.long'
map_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\available_specs.txt'
out_dir = r'C:\Users\felix\PycharmProjects\general\benchmarking\positive_set\data'

# load mapping
taxid_2_name = {}
name_2_taxid = {}
with open(map_path, 'r') as fh:
    for line in fh:
        data = line.strip().split('|')
        taxid = data[0]
        name = data[1].replace(' ', '_')
        taxid_2_name[taxid] = name
        name_2_taxid[name] = taxid

# initialize dataframe
spec_names = set()
fam_names = set()
with open(ortho_path, 'r') as fh:
    #skip header
    header = next(fh)
    for line in fh:
        data = line.strip().split()
        taxid = data[1].replace('ncbi', '')
        spec_names.add(taxid_2_name[taxid])
        fam_names.add(data[0])
spec_names = sorted(list(spec_names))
fam_names = sorted(list(fam_names))
raw_df = pd.DataFrame(columns=spec_names, index=fam_names)
# print(df)

def fill_df(path, df, with05=True):
    ret_df = df.copy()
    with open(path, 'r') as fh:
        # skip header
        header = next(fh)
        for line in fh:
            data = line.strip().split()
            if not with05 and len(data) == 4:
                if data[3] == '0':
                    continue
            taxid = data[1].replace('ncbi', '')
            spec_name = taxid_2_name[taxid]
            fam_name = data[0]
            ret_df.loc[fam_name, spec_name] = 1
    ret_df = ret_df.fillna(0)
    return ret_df


def find_missing(df1, df2):
    # subtracting the true positive set -> negative values are bad
    df = df1.subtract(df2)
    s_df = df[df.isin([-1.0]).any(1)]
    bo = df.isin([-1.0])
    # print(s_df)

    miss_spec = bo.sum().sort_values(ascending=False)
    # miss_spec = miss_spec[miss_spec != 0]
    miss_fam = bo.sum(1).sort_values(ascending=False)
    # miss_fam = miss_fam[miss_fam != 0]
    return miss_spec, miss_fam, s_df

# cmsearch cutoff 0.5
ortho_df = fill_df(ortho_path, raw_df, with05=True)
mirgene_df = fill_df(mirgene_path, raw_df, with05=True)
five_miss_spec, five_miss_fam, five_summary = find_missing(ortho_df, mirgene_df)
five_summary.to_csv(f'{out_dir}/summary_missed_TPs.tsv', sep='\t')

sum_dict = {}
for col in five_summary.columns:
    is_mis = five_summary[col][five_summary[col] == -1.0]
    sum_dict[col] = sorted(list(is_mis.index))
print(sum_dict)
with open(f'{out_dir}/missing_fams_per_species.json', 'w') as of:
    json.dump(sum_dict, of)

# # cmsearch cutoff 0.6
# ortho_df = fill_df(ortho_path, raw_df, with05=False)
# mirgene_df = fill_df(mirgene_path, raw_df, with05=False)
# six_miss_spec, six_miss_fam = find_missing(ortho_df, mirgene_df)
#
# # combine
# spec_df = pd.concat([six_miss_spec, five_miss_spec], axis=1)
# spec_df.columns = ['06cutoff', '05cutoff']
# fam_df = pd.concat([six_miss_fam, five_miss_fam], axis=1)
# fam_df.columns = ['06cutoff', '05cutoff']
#
# # filter and sort
# fam_df = fam_df[~fam_df.isin([0]).all(1)]
# fam_df = fam_df.sort_values(by=['06cutoff', '05cutoff'], ascending=False)
# spec_df = spec_df[~spec_df.isin([0]).all(1)]
# spec_df = spec_df.sort_values(by=['06cutoff', '05cutoff'], ascending=False)
#
# # print results
# print(fam_df)
# print(spec_df)
#
# # save results
# spec_df.to_csv(f'{out_dir}/spec_missing_TP.tsv', sep='\t')
# fam_df.to_csv(f'{out_dir}/fam_missing_TP.tsv', sep='\t')


