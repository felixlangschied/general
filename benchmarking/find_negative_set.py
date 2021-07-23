import os
import subprocess as sp

# over05_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\05overview.txt'
over05_path = '/home/felixl/PycharmProjects/general/benchmarking/data/05overview.txt'
outpath = '/home/felixl/PycharmProjects/general/benchmarking/data/RNAfold'
if not os.path.isdir(outpath):
    os.mkdir(outpath)

with open(over05_path, 'r') as fh:
    for line in fh:
        ipath = line.split()[0].replace('mirGeneDB', 'mirgenedb')
        print(ipath)
        fname = '_'.join(ipath.split('/')[-2:]).replace('_pre_orthologs.fa', '')
        # print(f'{outpath}/{fname}')

        fold_cmd = f'RNAfold -i {ipath} -o {outpath}/{fname}'
        sp.run(fold_cmd, shell=True)


        break
