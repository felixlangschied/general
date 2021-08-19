import pandas as pd
import sys
import numpy as np

mirbh_path = '/home/felixl/project/ncOrtho/benchmark/filtered_mirbh/analysis/phyloprofile/pp.long'
ncortho_path = '/home/felixl/project/ncOrtho/benchmark/ncortho/analysis/phyloprofile/pp.long'
outf = '/home/felixl/project/ncOrtho/benchmark/mirna_ncortho_vs_mirbh.long'


def read_pp(path):
    out_dict = {'id': []}
    with open(path, 'r') as fh:
        next(fh)
        for line in fh:
            linedata = line.strip().split()
            idstr = '_'.join(linedata[:2])
            out_dict['id'].append(idstr)
    df = pd.DataFrame.from_dict(out_dict)
    df['found'] = 1
    df = df.set_index('id')
    return df


mirbh_df = read_pp(mirbh_path)
ncortho_df = read_pp(ncortho_path)

comp_df = mirbh_df.join(ncortho_df, how='outer', lsuffix='_mirbh', rsuffix='_ncortho')
# print(df[df.isnull().any(axis=1)])
# print(df)

with open(outf, 'w') as of:
    of.write('geneID\tncbiID\torthoID\tmirbh_only\n')
    for row in comp_df.iterrows():
        # print(list(row)[1])
        entry = list(row)[1]
        # print(entry['found_mirbh'])
        # print(entry['found_ncortho'])
        # if entry['found_mirbh'] and entry['found_ncortho']:
        #     qualifier = '0.5'
        # elif entry['found_mirbh']:
        #     qualifier = '1'
        # elif entry['found_ncortho']:
        #     qualifier = '0'
        # print(type(entry['found_ncortho']))
        if np.isnan(entry['found_ncortho']):
            mirbh_only = '1'
        elif np.isnan(entry['found_mirbh']):
            mirbh_only = '0'
        else:
            mirbh_only = '0.5'

        mirna, taxid = str(entry.name).split('_')
        fam = '-'.join(mirna.split('-')[1:3])

        print(fam, taxid, mirna, mirbh_only)
        # outstr = f'{fam}\t{taxid}\t{mirna}\t{mirbh_only}\n'
        outstr = f'{mirna}\t{taxid}\t{mirna}\t{mirbh_only}\n'
        of.write(outstr)
