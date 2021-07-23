import os
import subprocess as sp

# over05_path = r'C:\Users\felix\PycharmProjects\general\benchmarking\data\05overview.txt'
over05_path = '/home/felixl/PycharmProjects/general/benchmarking/data/05overview.txt'
outpath = '/home/felixl/PycharmProjects/general/benchmarking/data/RNAfold'
if not os.path.isdir(outpath):
    os.mkdir(outpath)


outdict = {}
allscores = []
with open(over05_path, 'r') as fh:
    for line in fh:
        ipath = line.split()[0].replace('mirGeneDB', 'mirgenedb')
        specname = ipath.split('/')[7]
        mirname = ipath.split('/')[-1].replace('_pre_orthologs.fa', '')

        # run RNAfold
        fold_cmd = f'RNAfold -i {ipath}'
        res = sp.run(fold_cmd, shell=True, capture_output=True)
        out = res.stdout.decode('utf-8')
        score = out.split(' (')[1].replace(')', '')
        # fill dict
        if specname not in outdict:
            outdict[specname] = {}
        outdict[specname][mirname] = score
        allscores.append(score)
        print(out)
        break
print(outdict)