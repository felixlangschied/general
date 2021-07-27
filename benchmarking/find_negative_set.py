import os
import subprocess as sp
import json

# over05_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\05overview.txt'
over_path = '/home/felixl/PycharmProjects/general/benchmarking/data/06overview.txt'
outfile = '/home/felixl/PycharmProjects/general/benchmarking/data/all06_RNAfold.json'


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

df_dict = {'species': [], 'mirna_coorth': [], 'mirna': [], 'mirna_fam': [], 'score': [], 'seq': [], 'scheme': []}
# out_stdout = f'{outpath}/all05_RNAfold_results.txt'
# with open(over05_path, 'r') as fh, open(out_stdout, 'w') as of:
with open(over_path, 'r') as fh:
    for c, line in enumerate(fh):
        print(f'{c} / 15880')
        ipath = line.split()[0].replace('mirGeneDB', 'mirgenedb')
        specname = ipath.split('/')[7]
        # run RNAfold
        fold_cmd = f'RNAfold -i {ipath}'
        res = sp.run(fold_cmd, shell=True, capture_output=True)
        out = res.stdout.decode('utf-8').strip()
        data = out.split('\n')
        data = chunks(data, 3)
        for chunk in list(data):
            # print(chunk)
            header, seq, schemscore = chunk
            mirname = header.split('|')[0].replace('>', '')
            mirna = mirname.split('_')[0]
            mirfam = '-'.join(mirna.split('-')[0:3])
            # print(schemscore.split(' ('))
            scheme, score = schemscore.split(' (')
            score = score.replace(')', '')
            # fill dict
            df_dict['species'].append(specname)
            df_dict['mirna_coorth'].append(mirname)
            df_dict['mirna'].append(mirna)
            df_dict['mirna_fam'].append(mirfam)
            df_dict['score'].append(score)
            df_dict['seq'].append(seq)
            df_dict['scheme'].append(scheme)
            # remove tmp file
            cur_dir = os.getcwd()
            to_rmv = cur_dir + '/' + header.replace('>', '').replace('|', '_') + '_ss.ps'
            if os.path.isfile(to_rmv):
                os.remove(to_rmv)
        # if c == 100:
        #     break
# print(df_dict)

with open(outfile, 'w') as of:
    json.dump(df_dict, of)
