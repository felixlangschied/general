import subprocess as sp
import tempfile
import json
import os

# ref_path = '/share/project/felixl/ncOrtho/data/mirGeneDB/data/mirgenedb.tsv'
# ref_path = '/share/project/felixl/ncOrtho/data/mmu_ref_core/mirnas/mmu_highc_felix.tsv'
ref_path = '/share/project/felixl/ncOrtho/data/hsa_ref/mirnas/hsa_highc_505.tsv'
out_file = '/home/felixl/PycharmProjects/general/benchmarking/data/hsa_mirbase_RNAfold.json'

df_dict = {'species': [], 'mirna': [], 'mirna_fam': [], 'score': [], 'seq': [], 'scheme': []}
with open(ref_path, 'r') as fh:
    for line in fh:
        mirna, chrom, start, end, strand, seq, mat = line.strip().split('\t')
        specname = 'human_ref'
        with tempfile.NamedTemporaryFile(mode='w+') as fp:
            fp.write(f'>{mirna}\n{seq}\n')
            fp.seek(0)
            fold_cmd = f'RNAfold -i {fp.name}'
            res = sp.run(fold_cmd, shell=True, capture_output=True)
        out = res.stdout.decode('utf-8').strip()
        data = out.split('\n')
        schemscore = data[-1]
        mirfam = '-'.join(mirna.split('-')[0:3])
        scheme, score = schemscore.split(' (')
        score = score.replace(')', '')
        df_dict['species'].append(specname)
        df_dict['mirna'].append(mirna)
        df_dict['mirna_fam'].append(mirfam)
        df_dict['score'].append(score)
        df_dict['seq'] = seq
        df_dict['scheme'].append(scheme)

        if os.path.isfile(f'{mirna}_ss.ps'):
            os.remove(f'{mirna}_ss.ps')


with open(out_file, 'w') as of:
    json.dump(df_dict, of)