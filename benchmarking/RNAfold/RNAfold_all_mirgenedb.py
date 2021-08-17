import tempfile
import subprocess as sp
import os
import json

all_file = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb/ALL-pre.fas'
spec_trash = '/home/felixl/PycharmProjects/general/benchmarking/data/mirgenedb_specs_abb.txt'
pres_path = '/home/felixl/PycharmProjects/general/benchmarking/data/present.txt'
# out_file = '/home/felixl/PycharmProjects/general/benchmarking/data/all19_mirgenedb_rnafold'
out_file = '/home/felixl/PycharmProjects/general/benchmarking/RNAfold/data/mirgenedb_rnafold.json'


with open(pres_path, 'r') as fh:
    present = fh.read().split('\n')
    present.pop()

spec_map = {}
with open(spec_trash, 'r') as fh:
    for line in fh:
        data = line.strip().split('\t')
        if data[0]:
            name = data[1].replace(' ', '_')
            abb = data[2].capitalize()
            spec_map[abb] = name

df_dict = {'species': [], 'mirna': [], 'mirna_fam': [], 'score': [], 'seq': [], 'scheme': []}
with open(all_file, 'r') as fh:
    # mirna =
    # specname =
    # seq =
    for line in fh:
        if line.startswith('>'):
            abb = line.split('-')[0].replace('>', '')
            specname = spec_map[abb]
            # if specname not in present:
            #     continue
            mirna = '-'.join(line.strip().split('-')[1:])
            mirfam = '-'.join(line.strip().split('-')[1:3]).replace('_pre', '')
        else:
            # if specname not in present:
            #     continue
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
                df_dict['species'].append(specname)
                df_dict['mirna'].append(mirna)
                df_dict['mirna_fam'].append(mirfam)
                df_dict['score'].append(score)
                df_dict['seq'].append(seq)
                df_dict['scheme'].append(scheme)

                if os.path.isfile(f'{mirna}_ss.ps'):
                    os.remove(f'{mirna}_ss.ps')
        # print(line.strip())

with open(out_file, 'w') as of:
    json.dump(df_dict, of)