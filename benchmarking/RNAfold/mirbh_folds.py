import subprocess as sp
import glob
import os
import tempfile
import json

# data_dir = r'C:\Users\felix\PycharmProjects\general\benchmarking\RNAfold\mirbh_data'
data_dir = '/home/felixl/project/ncOrtho/benchmark/mirbh'
outf = '/home/felixl/project/ncOrtho/benchmark/mirbh/mirbh_folds.json'

files = glob.glob(f'{data_dir}/*.fa')

df_dict = {'species': [], 'mirna_coorth': [], 'mirna': [], 'mirna_fam': [], 'score': [], 'seq': [], 'scheme': []}
for file in files:
    species = file.split(os.sep)[-1].replace('.fa', '')
    # print(species)
    with open(file, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                mirna = '-'.join(line.strip().split('-')[1:])
                mirfam = '-'.join(line.strip().split('-')[1:3]).replace('_pre', '')
            else:
                seq = line.strip()
                # do RNAfold
                with tempfile.NamedTemporaryFile(mode='w+') as fp:
                    fp.write(f'>{mirna}\n{seq}\n')
                    fp.seek(0)
                    fold_cmd = f'RNAfold -i {fp.name}'
                    res = sp.run(fold_cmd, shell=True, capture_output=True)

                    out = res.stdout.decode('utf-8').strip()
                    folddata = out.split('\n')
                    schemscore = folddata[-1]
                    scheme, score = schemscore.split(' (')
                    score = score.replace(')', '')
                    df_dict['species'].append(species)
                    df_dict['mirna'].append(mirna)
                    df_dict['mirna_fam'].append(mirfam)
                    df_dict['score'].append(score)
                    df_dict['seq'] = seq
                    df_dict['scheme'].append(scheme)

                    if os.path.isfile(f'{mirna}_ss.ps'):
                        os.remove(f'{mirna}_ss.ps')
with open(outf, 'w') as of:
    json.dump(df_dict, of)